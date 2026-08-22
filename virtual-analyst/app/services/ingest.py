"""mzML -> ClickHouse (points) + SQLite (channel metadata).

It writes EXACTLY the way the production system writes:
  * `chromatogram_points` - the ready-made chromatograms shipped in the mzML
    (TIC, BasePeak). Each becomes its own channel, with `ion` filled in and
    scan_filter_id NULL.
  * `scans` + `spectra_points` - the run centroid by centroid, grouped per
    channel, where a channel is (sample, scan filter). This is what the pipeline
    reads to build an XIC for any m/z.

`scan_index` is counted PER CHANNEL, starting at 1 - it is not the mzML's global
index. That is what makes the scans x spectra_points LEFT JOIN line up.
"""

from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path

from pyteomics import mzml

from app import config
from app.db.clickhouse import insert_retry, lake_ch
from app.db.sqlite import CompatConnection

_CONTROL_PREFIXES = ("CQCL_REINJ", "CQCL", "CQN", "CQHP_REINJ", "CQHP")

SCANS_COLS = ["year", "sample_id", "channel_id", "scan_index", "rt", "ms_level"]
POINTS_COLS = ["year", "sample_id", "channel_id", "scan_index", "mz", "intensity"]
CHROM_COLS = ["year", "sample_id", "channel_id", "rt", "intensity"]


def is_control_sample(sample_name: str) -> bool:
    """Guess whether a sample is a control from its name prefix.

    This only pre-ticks a box on screen - who decides which sample is the CQCL,
    the CQCL_reinj and the CQN is the user, in the task form.
    """
    name = sample_name.strip().upper()
    return any(name.startswith(prefix) for prefix in _CONTROL_PREFIXES)


def year_from_batch(batch: str) -> int:
    """Extract the year from the batch name (`20260812_0442_26` -> 2026).

    With no recognisable prefix it falls back to the current year. The year
    matters because the ClickHouse tables are partitioned by it - getting it
    wrong corrupts nothing, it only costs partition pruning.
    """
    if re.match(r"^(19|20)\d{6}", batch.strip()):
        return int(batch[:4])
    return datetime.now().year


# ---------------------------------------------------------------------------
# Metadata (lake SQLite)
# ---------------------------------------------------------------------------

def get_or_create_scan_filter_id(conn: CompatConnection, scan_filter: str, cache: dict) -> int:
    """Resolve (or create) the id of the scan filter text in the dimension."""
    cached = cache.get(scan_filter)
    if cached is not None:
        return cached
    cursor = conn.cursor()
    cursor.execute("INSERT OR IGNORE INTO scan_filters (scan_filter) VALUES (%s)", (scan_filter,))
    cursor.execute("SELECT id FROM scan_filters WHERE scan_filter = %s", (scan_filter,))
    filter_id = int(cursor.fetchone()[0])
    cursor.close()
    cache[scan_filter] = filter_id
    return filter_id


def insert_channel(
    conn: CompatConnection,
    *,
    sample_id: int,
    scan_filter: str | None,
    ion: str | None,
    year: int,
    screening: str | None,
    cache: dict,
) -> int:
    scan_filter_id = get_or_create_scan_filter_id(conn, scan_filter, cache) if scan_filter else None
    cursor = conn.cursor()
    cursor.execute(
        """
        INSERT INTO sample_channels (sample_id, scan_filter_id, ion, year, screening)
        VALUES (%s, %s, %s, %s, %s)
        """,
        (sample_id, scan_filter_id, ion, int(year), screening),
    )
    channel_id = int(cursor.lastrowid)
    cursor.close()
    return channel_id


# ---------------------------------------------------------------------------
# Reading the mzML
# ---------------------------------------------------------------------------

def detect_mzml_contents(mzml_path: Path) -> tuple[bool, bool]:
    """(has ready-made chromatograms, has spectra)."""
    has_chrom = False
    has_spec = False

    with mzml.MzML(str(mzml_path)) as reader:
        for _ in reader.iterfind("chromatogram"):
            has_chrom = True
            break

    with mzml.read(str(mzml_path)) as reader:
        for spectrum in reader:
            mzs = spectrum.get("m/z array")
            intensities = spectrum.get("intensity array")
            if mzs is not None and intensities is not None and len(mzs) > 0:
                has_spec = True
                break

    return has_chrom, has_spec


def ingest_chromatograms(ch, conn, mzml_path: Path, *, sample_id: int, year: int, screening: str | None, cache: dict) -> int:
    channels: dict[str, int] = {}
    total = 0

    with mzml.MzML(str(mzml_path)) as reader:
        for chromatogram in reader.iterfind("chromatogram"):
            chrom_id = chromatogram.get("id", "")
            if not chrom_id:
                continue

            if chrom_id not in channels:
                channels[chrom_id] = insert_channel(
                    conn, sample_id=sample_id, scan_filter=None, ion=chrom_id,
                    year=year, screening=screening, cache=cache,
                )

            rt_array = chromatogram.get("time array")
            intensity_array = chromatogram.get("intensity array")
            if rt_array is None or intensity_array is None:
                continue

            channel_id = channels[chrom_id]
            rows = []
            for rt, intensity in zip(rt_array, intensity_array):
                rows.append([year, sample_id, channel_id, float(rt), float(intensity)])
                if len(rows) >= config.BATCH_SCANS:
                    insert_retry(ch, "chromatogram_points", rows, CHROM_COLS)
                    total += len(rows)
                    rows.clear()

            if rows:
                insert_retry(ch, "chromatogram_points", rows, CHROM_COLS)
                total += len(rows)

    return total


def ingest_scans_and_spectra(ch, conn, mzml_path: Path, *, sample_id: int, year: int, screening: str | None, cache: dict) -> tuple[int, int]:
    channels: dict[str, int] = {}
    counters: dict[str, int] = {}
    scan_rows: list[list] = []
    point_rows: list[list] = []
    n_scans = 0
    n_points = 0

    with mzml.read(str(mzml_path)) as reader:
        for spectrum in reader:
            ms_level = spectrum.get("ms level")
            if ms_level is None:
                continue

            scan = spectrum.get("scanList", {}).get("scan", [{}])[0]
            rt = scan.get("scan start time")
            if rt is None:
                continue

            # With no `filter string` (an mzML from another vendor) every
            # spectrum lands in a single channel. The method then needs
            # NO_FILTER configured as its scan filter to find that channel.
            scan_filter = scan.get("filter string") or "NO_FILTER"
            mzs = spectrum.get("m/z array")
            intensities = spectrum.get("intensity array")
            if mzs is None or intensities is None:
                continue

            if scan_filter not in channels:
                channels[scan_filter] = insert_channel(
                    conn, sample_id=sample_id, scan_filter=str(scan_filter), ion=None,
                    year=year, screening=screening, cache=cache,
                )
                counters[scan_filter] = 0

            counters[scan_filter] += 1
            scan_index = counters[scan_filter]
            channel_id = channels[scan_filter]
            scan_rows.append([year, sample_id, channel_id, scan_index, float(rt), int(ms_level)])

            for mz, intensity in zip(mzs, intensities):
                point_rows.append([year, sample_id, channel_id, scan_index, float(mz), float(intensity)])

            # The cut is on the number of POINTS: one full-scan spectrum holds
            # thousands of centroids, so counting scans would not bound memory.
            if len(point_rows) >= config.BATCH_POINTS:
                insert_retry(ch, "scans", scan_rows, SCANS_COLS)
                n_scans += len(scan_rows)
                scan_rows.clear()
                insert_retry(ch, "spectra_points", point_rows, POINTS_COLS)
                n_points += len(point_rows)
                point_rows.clear()

    if scan_rows:
        insert_retry(ch, "scans", scan_rows, SCANS_COLS)
        n_scans += len(scan_rows)
    if point_rows:
        insert_retry(ch, "spectra_points", point_rows, POINTS_COLS)
        n_points += len(point_rows)

    return n_scans, n_points


# ---------------------------------------------------------------------------
# Cleanup (re-ingestion)
# ---------------------------------------------------------------------------

def purge_sample_points(ch, sample_id: int) -> None:
    """Delete a sample's points before ingesting it again.

    Mutations are asynchronous in ClickHouse; `mutations_sync=2` waits for them
    to finish so re-ingestion cannot end up duplicating points.
    """
    for table in ("spectra_points", "scans", "chromatogram_points"):
        ch.command(
            f"ALTER TABLE {table} DELETE WHERE sample_id = {int(sample_id)}",
            settings={"mutations_sync": 2},
        )


def purge_sample_channels(conn: CompatConnection, sample_id: int) -> None:
    cursor = conn.cursor()
    cursor.execute("DELETE FROM sample_channels WHERE sample_id = %s", (sample_id,))
    cursor.close()
    conn.commit()


# ---------------------------------------------------------------------------
# Per-sample orchestration
# ---------------------------------------------------------------------------

def ingest_mzml(conn: CompatConnection, mzml_path: Path, *, sample_id: int, year: int, screening: str | None) -> dict:
    """Read the whole mzML and store it. Returns a summary for the on-screen log."""
    ch = lake_ch()
    cache: dict[str, int] = {}
    try:
        purge_sample_points(ch, sample_id)
        purge_sample_channels(conn, sample_id)

        has_chrom, has_spec = detect_mzml_contents(mzml_path)
        summary = {
            "chromatogram_points": 0,
            "scans": 0,
            "spectra_points": 0,
            "ingested_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }

        if has_chrom:
            summary["chromatogram_points"] = ingest_chromatograms(
                ch, conn, mzml_path, sample_id=sample_id, year=year,
                screening=screening, cache=cache,
            )
        if has_spec:
            n_scans, n_points = ingest_scans_and_spectra(
                ch, conn, mzml_path, sample_id=sample_id, year=year,
                screening=screening, cache=cache,
            )
            summary["scans"] = n_scans
            summary["spectra_points"] = n_points

        conn.commit()
        if not has_spec:
            # No spectra means no XIC, and with no XIC the pipeline has nothing
            # to analyse. Better to fail here, during ingestion, than to return
            # INVALID for every substance later.
            raise RuntimeError(
                f"{mzml_path.name} has no spectra - only ready-made chromatograms. "
                f"Without spectra there is no way to extract an XIC."
            )
        return summary
    finally:
        ch.close()


def scan_filters_of_sample(conn: CompatConnection, sample_id: int) -> list[str]:
    """The scan filters observed in a sample that has already been ingested.

    The method screen uses this to offer the EXACT text: the match against
    `method_substances.scan_filter` is string for string, and one extra space
    makes the substance come out INVALID without saying why.
    """
    cursor = conn.cursor()
    cursor.execute(
        """
        SELECT sf.scan_filter
        FROM sample_channels sc
        JOIN scan_filters sf ON sf.id = sc.scan_filter_id
        WHERE sc.sample_id = %s
        ORDER BY sf.scan_filter
        """,
        (sample_id,),
    )
    filters = [row[0] for row in cursor.fetchall()]
    cursor.close()
    return filters
