"""
data_provider.py - DbDataProvider, the "virtual MSFileReader".

FAITHFUL PORT of the module of the same name in the production system. The only
difference is where the metadata comes from: there `sample_channels` lives in
MariaDB, here in SQLite (through the compatibility layer in app/db/sqlite.py,
which accepts the same `%s` dialect). The calibration against MSFileReader - the
RT edge tolerance, the smoothing sigma, the LEFT JOIN that turns a scan with no
point into a zero - was kept character for character: it is what guarantees an
XIC identical to Xcalibur's.

WHY IT EXISTS:
    The original chemistry pipeline (verify_two_peaks,
    evaluate_against_reference_fast, and friends) was written assuming every
    sample/control is an MSFileReader object opened on a Thermo .raw file. Every
    chromatogram extraction goes through `raw.GetChroData(...)`.

    To reuse ALL of that pipeline while reading from the database (instead of
    the .raw), we need a class that LOOKS exactly like MSFileReader but fetches
    its points from ClickHouse. That is what this class does.

DROP-IN REPLACEMENT - it implements the minimal API the pipeline uses:
    raw.GetChroData(startTime, endTime, massRange1, scanFilter, smoothingType, smoothingValue)
        -> queries scans + spectra_points in ClickHouse, aggregates per scan,
           applies smoothing, returns [[rt_array, intensity_array]] (same shape).
    raw.SetCurrentController(...)   # no-op (only relevant for multi-detector .raw)
    raw.SetMassTolerance(...)       # adjusts `mass_tolerance_ppm` at runtime
    raw.Close()                     # no-op

CALIBRATION against MSFileReader:
    Thermo's `GetChroData` does three things ClickHouse does not do by default:
      1. It includes scans at the UPPER EDGE of the RT window because of Float32
         rounding -> handled by RT_EDGE_TOLERANCE.
      2. It applies gaussian smoothing with a specific sigma
         -> calibrated in apply_thermo_smoothing (sigma = (N-1)/5).
      3. It applies the mass tolerance via `SetMassTolerance(units=1)`
         -> controlled by `mass_tolerance_ppm`. For the TVII method leave it
            None (the mass_range already has the 6 ppm baked in).

    Validated against MSFileReader over 9 batches x ~2200 pairs: median
    RMSE/peak below 0.4% on the sample chromatograms.
"""

from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np
from scipy.ndimage import gaussian_filter1d


# Edge tolerance on the RT axis (in minutes) when building the XIC.
# It mirrors MSFileReader's behaviour, which includes scans at the "edge" of the
# [startTime, endTime] window (absorbing Float32 rounding and Thermo's own
# inclusion convention).
RT_EDGE_TOLERANCE = 0.001  # 60 ms - just enough to absorb Float32 rounding error


def parse_mass_range(mass_range: str) -> Tuple[float, float]:
    """Parse the `mass_range` field configured in method_substances.

    Expected format: 'min-max' (e.g. '279.14750510492-279.15085489508').
    A single value is also accepted (legacy) - it becomes the degenerate
    interval [v, v].
    """
    if mass_range is None:
        return (0.0, 0.0)
    s = str(mass_range).strip()
    if "-" not in s:
        v = float(s)
        return (v, v)
    lo, hi = s.split("-", 1)
    return (float(lo.strip()), float(hi.strip()))


def apply_thermo_smoothing(intensity: np.ndarray, smoothing_type: int, smoothing_value: int) -> np.ndarray:
    """Replicate the smoothing applied by `MSFileReader.GetChroData`.

    Thermo accepts these 3 modes through `smoothingType`:
      0 = none (passthrough)
      1 = Boxcar - uniform moving average over a `smoothing_value`-point window
      2 = Gaussian - centred gaussian kernel, sigma derived from `smoothing_value`

    For Gaussian (type=2), `sigma` was calibrated empirically: we ran the
    pipeline against dozens of (sample, .raw) pairs with smoothing_value=7 and
    tuned sigma until `scipy.ndimage.gaussian_filter1d` reproduced the
    MSFileReader output with a residual error below 1% at the peaks. The
    resulting formula:

        sigma = (smoothing_value - 1) / 5.0

    `truncate=4.0` sets how far the kernel is evaluated (4 sigmas - the classic
    gaussian cut-off). `mode="nearest"` avoids pulling the XIC edges to zero.

    Conceptual reference: http://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
    """
    # Trivial cases: not enough data, or smoothing turned off -> passthrough.
    if intensity.size < 3 or not smoothing_value or smoothing_value <= 0:
        return intensity

    n = int(smoothing_value)

    if smoothing_type == 1:
        # Boxcar - uniform kernel [1/n, 1/n, ...]
        kernel = np.ones(n, dtype=float) / float(n)
        return np.convolve(intensity, kernel, mode="same")

    if smoothing_type == 2:
        # Gaussian - sigma calibrated against MSFileReader (see the docstring).
        if n <= 1:
            return intensity
        sigma = max((n - 1) / 5.0, 1e-3)
        return gaussian_filter1d(intensity, sigma=sigma, truncate=4.0, mode="nearest")

    # Unknown type - return it unsmoothed.
    return intensity


class DbDataProvider:
    """Stand-in for `MSFileReader` that pulls chromatograms from the database.

    One instance represents ONE sample (a row in `samples`). The original
    pipeline has 4 of these per processed sample:
        raw_sample, raw_cqcl, raw_cqcl_reinj, raw_cqn

    Every `GetChroData(scanFilter, ...)` call resolves which `channel_id` in
    `sample_channels` corresponds to that scan_filter for that sample (and
    caches it).
    """

    def __init__(
        self,
        mariadb_conn,
        clickhouse_client,
        sample_id: int,
        sample_label: str = "",
        mass_tolerance_ppm: Optional[float] = None,
        year: Optional[int] = None,
    ):
        """
        :param mariadb_conn: lake SQLite connection (reads sample_channels).
        :param clickhouse_client: clickhouse_connect client for the `lake`
            database (reads scans + spectra_points).
        :param sample_id: id of the `samples` row this instance represents.
        :param sample_label: human-readable label (the sample name) - logs only.
        :param mass_tolerance_ppm: tolerance in ppm applied on top of the
            `massRange1` received in GetChroData. It works like MSFileReader's
            `SetMassTolerance(units=1)`.

            WARNING: do NOT use it for the TVII method - the `mass_range`
            configured in `method_substances` already has the 6 ppm tolerance
            baked in (e.g. m/z 223.119 -> "223.1176-223.1203" = +/-6 ppm).
            Applying ppm here DOUBLES the tolerance and picks up extra m/z
            points that MSFileReader.GetChroData would never see.

            Calibrated validation: leaving it `None` matches MSFileReader within
            1% RMSE across 21k+ tested pairs.
        """
        self.mdb = mariadb_conn
        self.ch = clickhouse_client
        self.sample_id = int(sample_id)
        self.sample_label = sample_label or f"sample_id={sample_id}"
        self.mass_tolerance_ppm = mass_tolerance_ppm
        # The sample's year. When set, it becomes a `year = X` filter in the
        # ClickHouse query and enables partition pruning (scans, spectra_points
        # and chromatogram_points are partitioned by year in this schema).
        # Without it the query scans every partition - it works, but slower.
        self.year: Optional[int] = int(year) if year is not None else None
        # Cache (scan_filter -> channel_id). The pipeline asks for the same
        # scan_filter many times (once per substance) - caching avoids N
        # redundant metadata lookups for the same sample.
        self._channel_cache: dict[str, Optional[int]] = {}

    # ------------------------------------------------------------------------
    # The minimal MSFileReader API the pipeline uses.
    # These are no-ops here because their equivalents make no sense against a
    # database - but they have to EXIST to keep the signature compatible.
    # ------------------------------------------------------------------------

    def SetCurrentController(self, *args, **kwargs):
        # In a .raw an MS run can have multiple "detectors" (controllers). The
        # lake only stores the MS channel that matters, so there is nothing to
        # select.
        return None

    def SetMassTolerance(self, *args, **kwargs):
        """Allow the tolerance to be adjusted at runtime (units=1 = ppm).

        The original pipeline calls this right after opening the .raw:
            raw.SetMassTolerance(userDefined=True, massTolerance=6.0, units=1)
        Here we accept the call and store the value - but the .env guidance is to
        leave it None for the TVII method (see the __init__ docstring).
        """
        ppm = kwargs.get("massTolerance")
        if ppm is None and len(args) >= 2:
            ppm = args[1]
        if ppm is not None:
            self.mass_tolerance_ppm = float(ppm)
        return None

    def Close(self, *args, **kwargs):
        # On a .raw this releases the file handle. Here there is no resource to
        # release (the ClickHouse/SQLite connection is shared by every
        # DbDataProvider in the batch - the orchestrator closes it).
        return None

    # ------------------------------------------------------------------------
    # Channel resolution: given a scan_filter, find which `channel_id` in
    # ClickHouse holds the matching spectral points.
    # ------------------------------------------------------------------------

    def _resolve_channel_id(self, scan_filter: str) -> Optional[int]:
        """Every (sample, scan_filter) combination becomes a row in
        `sample_channels` with its own `id`. The points in ClickHouse are stored
        referencing that `channel_id` - it is the link between metadata (SQLite)
        and data (ClickHouse).

        Cached in memory so we do not hit the metadata database every time.
        Returns None when the sample has no data for that scan_filter - in which
        case GetChroData returns an empty XIC.
        """
        if scan_filter in self._channel_cache:
            return self._channel_cache[scan_filter]
        cur = self.mdb.cursor()
        # The scan_filter text is normalised into the `scan_filters` dimension
        # and `sample_channels` references it by `scan_filter_id` (INT), so we
        # resolve it with a JOIN through the dimension.
        # The JOIN also excludes the direct chromatogram channels (BasePeak and
        # friends), which have a NULL scan_filter_id and are never used by the
        # pipeline.
        cur.execute(
            """
            SELECT sc.id
            FROM sample_channels sc
            JOIN scan_filters sf ON sf.id = sc.scan_filter_id
            WHERE sc.sample_id = %s AND sf.scan_filter = %s
            LIMIT 1
            """,
            (self.sample_id, scan_filter),
        )
        row = cur.fetchone()
        cur.close()
        cid = int(row[0]) if row else None
        self._channel_cache[scan_filter] = cid
        return cid

    # ------------------------------------------------------------------------
    # The main API - the equivalent of MSFileReader.GetChroData
    # ------------------------------------------------------------------------

    def GetChroData(
        self,
        startTime: float,
        endTime: float,
        massRange1: str,
        scanFilter: str,
        smoothingType: int = 0,
        smoothingValue: int = 0,
        **_ignored,
    ):
        """Extract a chromatogram (XIC) from the database in the SAME format as
        MSFileReader.

        Output: [[rt_array, intensity_array]] - a list of two lists, exactly what
        `MSFileReader.GetChroData` returns. The rest of the pipeline assumes that
        shape and would break if it changed.

        Internal steps:
          1. Resolve `channel_id` for the (sample, scanFilter) pair. If there is
             none -> return an empty XIC (a substance with no signal in that
             channel - common, not an error).
          2. Parse `massRange1` ('mz_min-mz_max').
          3. If `mass_tolerance_ppm` is set, widen the interval:
                 mz_min' = mz_min - mz_min*ppm/1e6
                 mz_max' = mz_max + mz_max*ppm/1e6
          4. Run the LEFT JOIN query in ClickHouse:
               - `scans` lists every scan of the channel in the RT range (X axis).
               - `spectra_points` filters m/z to the configured range (Y axis).
               - The LEFT JOIN guarantees that scans WITHOUT valid m/z points
                 become zero instead of a hole - that keeps the XIC continuous,
                 the same way Thermo does it.
          5. Apply the smoothing (gaussian/boxcar) calibrated against MSFileReader.

        RT edge tolerance: we extend ONLY the upper bound, by RT_EDGE_TOLERANCE.
        That reproduces MSFileReader's behaviour, which includes the last scan
        even when its RT is micro-rounded just above `endTime` (Float32). We do
        not extend the lower bound because `startTime` tends to be a round number
        in the method, and extending it would pull in a scan outside the window.
        """
        channel_id = self._resolve_channel_id(scanFilter)
        if channel_id is None:
            # No channel registered -> empty chromatogram. The pipeline treats
            # that as "substance with no valid data" and marks it
            # INVALID/NEGATIVE.
            return [[[], []]]

        mz_min, mz_max = parse_mass_range(massRange1)
        if mz_max < mz_min:
            mz_min, mz_max = mz_max, mz_min

        # Optional ppm widening (mirrors SetMassTolerance(units=1) on a .raw).
        # For the TVII method ppm stays None - the tolerance is already in
        # the mass_range.
        if self.mass_tolerance_ppm:
            ppm = self.mass_tolerance_ppm / 1e6
            mz_min = mz_min - mz_min * ppm
            mz_max = mz_max + mz_max * ppm

        # Asymmetric tolerance on the RT edges - upper bound only.
        rt_lo = float(startTime)
        rt_hi = float(endTime) + RT_EDGE_TOLERANCE

        # ClickHouse query:
        #   - GROUP BY scan: sums the intensity of every m/z point inside the
        #     m/z window (= the classic XIC).
        #   - LEFT JOIN: scans with no m/z point in the window become
        #     intensity=0 (keeping the RT grid continuous).
        #   - ORDER BY rt: guarantees the chromatogram comes out chronological.
        #   - The `year = X` filter on both tables enables partition pruning
        #     (the tables are partitioned by year in the lake schema). When the
        #     year is unavailable the clause is omitted and the query scans
        #     every partition.
        year_clause_s = f"AND s.year = {int(self.year)}" if self.year is not None else ""
        year_clause_p = f"AND p.year = {int(self.year)}" if self.year is not None else ""
        sql = f"""
            SELECT s.rt AS rt, sum(coalesce(p.intensity, 0)) AS intensity
            FROM scans s
            LEFT JOIN spectra_points p
                ON p.channel_id = s.channel_id
                AND p.scan_index = s.scan_index
                {year_clause_p}
                AND p.mz BETWEEN {float(mz_min)} AND {float(mz_max)}
            WHERE s.channel_id = {int(channel_id)}
              {year_clause_s}
              AND s.rt BETWEEN {rt_lo} AND {rt_hi}
            GROUP BY s.rt
            ORDER BY s.rt
        """
        result = self.ch.query(sql)
        rows = result.result_rows or []

        if not rows:
            return [[[], []]]

        # Convert to numpy (exactly what the pipeline expects to receive).
        rt = np.array([float(r[0]) for r in rows], dtype=np.float64)
        intensity = np.array([float(r[1]) for r in rows], dtype=np.float64)

        # Calibrated smoothing - the last step before handing back to the
        # pipeline.
        # NOTE: there is no "isolated spike" filter here any more. The spurious
        # centroids (lock mass, AGC overflow, internal calibration) are already
        # removed AT THE SOURCE by the ingestion, which uses ThermoRawFileParser
        # with `--excludeExceptionData`. The data in ClickHouse is already clean
        # - any isolated point left is legitimate sample signal.
        intensity = apply_thermo_smoothing(intensity, smoothingType, smoothingValue)

        # The exact format MSFileReader.GetChroData returns.
        return [[rt.tolist(), intensity.tolist()]]
