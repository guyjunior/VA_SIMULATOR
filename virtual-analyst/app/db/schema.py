"""Creates the whole schema: two SQLite files plus two ClickHouse databases.

    py -m app.db.schema    # creates everything and seeds the 4 analysis rules (idempotent)

The ClickHouse tables are a FAITHFUL COPY of the production ones (same engine,
same ORDER BY, same CODECs). The SQLite ones mirror the MariaDB tables, with the
types translated and the foreign keys kept.
"""

from __future__ import annotations

import sys

from app import config
from app.db.clickhouse import server_ch
from app.db.sqlite import analysis_conn, lake_conn

# ---------------------------------------------------------------------------
# ClickHouse
# ---------------------------------------------------------------------------

CH_LAKE_TABLES = [
    # Ready-made chromatograms that ship inside the mzML (TIC, BasePeak). The
    # analytical pipeline does not use them - they are there to eyeball the run.
    """
    CREATE TABLE IF NOT EXISTS {db}.chromatogram_points
    (
        year       UInt16,
        sample_id  UInt64,
        channel_id UInt64,
        rt         Float32 CODEC(ZSTD(3)),
        intensity  Float32 CODEC(ZSTD(3))
    )
    ENGINE = MergeTree
    PARTITION BY year
    ORDER BY (channel_id, rt)
    SETTINGS index_granularity = 8192
    """,
    # The run's RT grid: one row per scan. This is the X axis of the XIC - the
    # LEFT JOIN against it is what turns a scan with no point in the m/z window
    # into a zero instead of a hole.
    """
    CREATE TABLE IF NOT EXISTS {db}.scans
    (
        year       UInt16,
        sample_id  UInt64,
        channel_id UInt64,
        scan_index UInt32,
        rt         Float32 CODEC(ZSTD(3)),
        ms_level   UInt8,
        INDEX idx_rt rt TYPE minmax GRANULARITY 4
    )
    ENGINE = MergeTree
    PARTITION BY year
    ORDER BY (channel_id, scan_index)
    SETTINGS index_granularity = 8192
    """,
    # The centroids. This is the big table (tens of millions of rows per
    # sample). ORDER BY (channel_id, mz, scan_index) is what makes the XIC
    # cheap: the m/z window becomes a contiguous range scan.
    """
    CREATE TABLE IF NOT EXISTS {db}.spectra_points
    (
        year       UInt16,
        sample_id  UInt64,
        channel_id UInt64,
        scan_index UInt32,
        mz         Float32 CODEC(ZSTD(3)),
        intensity  Float32 CODEC(ZSTD(3)),
        INDEX idx_mz mz TYPE minmax GRANULARITY 1,
        INDEX idx_scan_index scan_index TYPE minmax GRANULARITY 4
    )
    ENGINE = MergeTree
    PARTITION BY year
    ORDER BY (channel_id, mz, scan_index)
    SETTINGS index_granularity = 8192
    """,
]

CH_ANALYSIS_TABLES = [
    # The 4 traces (sample/cqcl/cqcl_reinj/cqcn) for every substance x sample,
    # exactly as the pipeline saw them. This is what the results screen plots.
    """
    CREATE TABLE IF NOT EXISTS {db}.processing_sampleparameter
    (
        id UUID DEFAULT generateUUIDv4(),
        sample_result_id UInt64,
        sample_type String,
        retention_time Float64,
        intensity Float64,
        timestamp DateTime64(3, 'UTC')
    )
    ENGINE = MergeTree
    ORDER BY (sample_result_id, retention_time)
    PRIMARY KEY (sample_result_id, retention_time)
    """,
    """
    ALTER TABLE {db}.processing_sampleparameter
    ADD INDEX IF NOT EXISTS idx_sample_type sample_type TYPE set(100) GRANULARITY 1
    """,
]


def create_clickhouse() -> None:
    ch = server_ch()
    try:
        for db, ddls in (
            (config.CH_LAKE_DB, CH_LAKE_TABLES),
            (config.CH_ANALYSIS_DB, CH_ANALYSIS_TABLES),
        ):
            ch.command(f"CREATE DATABASE IF NOT EXISTS {db}")
            print(f"  OK database {db}")
            for ddl in ddls:
                ch.command(ddl.format(db=db))
            print(f"  OK tables in {db}")
    finally:
        ch.close()


# ---------------------------------------------------------------------------
# SQLite - lake
# ---------------------------------------------------------------------------

LAKE_TABLES = [
    # A "folder" here is the BATCH the user uploaded: the set of .raw files sent
    # in one go. In the production system it would be the instrument's folder.
    """
    CREATE TABLE IF NOT EXISTS folders (
        id                        INTEGER PRIMARY KEY AUTOINCREMENT,
        path                      TEXT    NOT NULL,
        screening                 TEXT,
        computer                  TEXT    NOT NULL DEFAULT 'upload',
        chromatographic_technique TEXT    NOT NULL DEFAULT 'LC-MS',
        created_at                TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP,
        deleted_at                TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS samples (
        id                       INTEGER PRIMARY KEY AUTOINCREMENT,
        folder_id                INTEGER NOT NULL REFERENCES folders(id),
        name                     TEXT    NOT NULL,
        batch                    TEXT    NOT NULL,
        year                     INTEGER NOT NULL,
        screening                TEXT,
        computer                 TEXT    NOT NULL DEFAULT 'upload',
        chromatography_technique TEXT    NOT NULL DEFAULT 'LC-MS',
        gender                   TEXT,
        competition              INTEGER,
        beta_block               INTEGER,
        modality                 TEXT,
        sport                    TEXT,
        is_control               INTEGER NOT NULL DEFAULT 0,
        -- Conversion/ingestion state: pending | converting | ingesting | ready | error
        status                   TEXT    NOT NULL DEFAULT 'pending',
        converter                TEXT,
        raw_path                 TEXT,
        mzml_path                TEXT,
        error                    TEXT,
        created_at               TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP
    )
    """,
    # The same per-ACQUISITION uniqueness the production system uses: (name,
    # batch) alone does not identify a sample - the same batch runs under more
    # than one screening, with repeated control names (CQN in TVII and TXVI).
    """
    CREATE UNIQUE INDEX IF NOT EXISTS uniq_sample_acquisition
        ON samples (name, batch, IFNULL(screening, ''), computer)
    """,
    "CREATE INDEX IF NOT EXISTS idx_samples_folder ON samples (folder_id)",
    "CREATE INDEX IF NOT EXISTS idx_samples_batch ON samples (batch)",
    # Dimension holding the scan filter text (~87 chars, repeated across
    # millions of points). Normalising it cuts the footprint and turns the
    # search's LIKE into a query over a small table.
    """
    CREATE TABLE IF NOT EXISTS scan_filters (
        id          INTEGER PRIMARY KEY AUTOINCREMENT,
        scan_filter TEXT NOT NULL UNIQUE
    )
    """,
    # The link between metadata (SQLite) and data (ClickHouse): the points in
    # ClickHouse reference channel_id.
    """
    CREATE TABLE IF NOT EXISTS sample_channels (
        id             INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_id      INTEGER NOT NULL REFERENCES samples(id),
        scan_filter_id INTEGER REFERENCES scan_filters(id),
        ion            TEXT,
        year           INTEGER,
        screening      TEXT,
        created_at     TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_channels_sample ON sample_channels (sample_id)",
    "CREATE INDEX IF NOT EXISTS idx_channels_lookup ON sample_channels (sample_id, scan_filter_id)",
]


# ---------------------------------------------------------------------------
# SQLite - analysis
# ---------------------------------------------------------------------------

ANALYSIS_TABLES = [
    """
    CREATE TABLE IF NOT EXISTS methods (
        id         INTEGER PRIMARY KEY AUTOINCREMENT,
        name       TEXT NOT NULL UNIQUE,
        version    TEXT NOT NULL DEFAULT '1',
        created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
        updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
        deleted_at TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS analysis_rules (
        id          INTEGER PRIMARY KEY AUTOINCREMENT,
        code        TEXT NOT NULL UNIQUE,
        name        TEXT NOT NULL,
        description TEXT,
        created_at  TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
        updated_at  TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
        deleted_at  TEXT
    )
    """,
    # CAREFUL: namesake of the lake's `scan_filters`, with different semantics.
    # Here it is the filter CONFIGURED in the method (column `name`); there it
    # is the text observed in the mzML (column `scan_filter`). The two must
    # match string for string, or the substance finds no channel and comes out
    # INVALID.
    """
    CREATE TABLE IF NOT EXISTS scan_filters (
        id         INTEGER PRIMARY KEY AUTOINCREMENT,
        name       TEXT NOT NULL UNIQUE,
        created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
        updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
        deleted_at TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS method_substances (
        id                      INTEGER PRIMARY KEY AUTOINCREMENT,
        name                    TEXT    NOT NULL,
        version                 INTEGER NOT NULL DEFAULT 1,
        method_id               INTEGER NOT NULL REFERENCES methods(id),
        analysis_rule_id        INTEGER REFERENCES analysis_rules(id),
        -- 'substance' | 'internal_standard' | 'substance_group'
        --
        -- COLLATE NOCASE is not decoration. In MariaDB the default collation is
        -- case-insensitive, so the pipeline's `WHERE type = 'substance'` matches
        -- the 'SUBSTANCE' stored by the production system. SQLite compares TEXT
        -- case-sensitively, so without NOCASE an imported uppercase value would
        -- match NOTHING and every task would come back with zero results - with
        -- no error to explain it.
        --
        -- 'substance_group' marks a substance that belongs to a group: the
        -- individual-substance query filters it out, and it is reached only
        -- through substance_group_memberships. Storing it as 'substance' would
        -- process it TWICE.
        type                    TEXT    NOT NULL COLLATE NOCASE,
        competition             INTEGER NOT NULL DEFAULT 0,
        beta_blocker            INTEGER NOT NULL DEFAULT 0,
        start_time              REAL    NOT NULL,
        end_time                REAL    NOT NULL,
        mz                      REAL    NOT NULL,
        mass_error_ppm          REAL    NOT NULL DEFAULT 6.0,
        scan_filter_id          INTEGER REFERENCES scan_filters(id),
        fortified_concentration REAL,
        cut_off                 REAL,
        r2_threshold            REAL,
        dtw_limit               REAL,
        any_points_ratio_limit  REAL,
        smoothing_type          INTEGER NOT NULL DEFAULT 2,
        smoothing_value         INTEGER NOT NULL DEFAULT 7,
        created_at              TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP,
        updated_at              TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP,
        deleted_at              TEXT
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_ms_method ON method_substances (method_id, type, deleted_at)",
    # Fragments of the same molecule. Every member gets its own row in
    # sample_results; an UPDATE then consolidates `result_group`.
    """
    CREATE TABLE IF NOT EXISTS substance_groups (
        id         INTEGER PRIMARY KEY AUTOINCREMENT,
        method_id  INTEGER NOT NULL REFERENCES methods(id),
        name       TEXT,
        created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
        updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS substance_group_memberships (
        id                  INTEGER PRIMARY KEY AUTOINCREMENT,
        substance_group_id  INTEGER NOT NULL REFERENCES substance_groups(id),
        method_substance_id INTEGER NOT NULL REFERENCES method_substances(id),
        created_at          TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
        updated_at          TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
    )
    """,
    # Which internal standard is the method's GATE and the cut-off reference.
    """
    CREATE TABLE IF NOT EXISTS cutoff_processing_groups (
        id                INTEGER PRIMARY KEY AUTOINCREMENT,
        method_id         INTEGER NOT NULL REFERENCES methods(id),
        internal_standard INTEGER NOT NULL REFERENCES method_substances(id),
        created_at        TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
        updated_at        TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
    )
    """,
    # One run of the pipeline over a set of samples plus the 3 controls.
    """
    CREATE TABLE IF NOT EXISTS sample_processing_tasks (
        id                       INTEGER PRIMARY KEY AUTOINCREMENT,
        samples                  TEXT    NOT NULL,
        method_id                INTEGER NOT NULL,
        batch                    TEXT    NOT NULL,
        computer                 TEXT,
        screening                TEXT,
        -- Frozen copy of the method at run time. Without it, editing a
        -- substance tomorrow would rewrite the meaning of yesterday's result.
        method_snapshot          TEXT    NOT NULL,
        internal_standards_audit TEXT,
        status                   TEXT    NOT NULL,
        result                   TEXT,
        cqcl_id                  INTEGER,
        cqcl_reinj_id            INTEGER,
        cqn_id                   INTEGER,
        with_controls            INTEGER NOT NULL DEFAULT 1,
        validate_is              INTEGER NOT NULL DEFAULT 1,
        is_auto_processed        INTEGER NOT NULL DEFAULT 0,
        feedback_batch           TEXT,
        reprocessed              INTEGER NOT NULL DEFAULT 0,
        created_at               TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP,
        updated_at               TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS sample_results (
        id                        INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_processing_task_id INTEGER NOT NULL REFERENCES sample_processing_tasks(id),
        sample_name               TEXT    NOT NULL,
        substance                 TEXT    NOT NULL,
        result                    TEXT    NOT NULL,
        result_validation         TEXT,
        concentration_value       REAL,
        result_group              TEXT,
        created_at                TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP,
        updated_at                TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_results_task ON sample_results (sample_processing_task_id)",
    "CREATE INDEX IF NOT EXISTS idx_results_sample ON sample_results (sample_processing_task_id, sample_name)",
]


# The 4 rules. The `code` is what the pipeline reads - changing the text here
# changes which function runs. See app/services/virtual_analyst.py.
SEED_RULES = [
    (
        "standard",
        "Standard (R2 + DTW)",
        "Gate over the CQCL peak window: SUSPECT when R2 > r2_threshold AND DTW < dtw_limit. "
        "If it comes out NEGATIVE against the CQCL, it tries again against the CQCL_reinj.",
    ),
    (
        "two_peaks",
        "Two peaks",
        "For substances eluting in two peaks (isomers/metabolites). Evaluates both CQCL peaks "
        "separately; one passing is enough.",
    ),
    (
        "any_points",
        "Height ratio",
        "No R2/DTW: compares the sample height against the control height over a fixed window "
        "centred on the CQCL peak, against any_points_ratio_limit.",
    ),
    (
        "no_analyze",
        "Do not interpret",
        "Visualisation-only substance: the XIC is extracted and stored, but never classified. "
        "Comes out as NOT_INTERPRETED.",
    ),
]


def create_sqlite() -> None:
    for label, conn_factory, ddls in (
        ("lake", lake_conn, LAKE_TABLES),
        ("analysis", analysis_conn, ANALYSIS_TABLES),
    ):
        conn = conn_factory()
        try:
            cursor = conn.cursor()
            for ddl in ddls:
                cursor.execute(ddl)
            conn.commit()
            print(f"  OK SQLite {label}")
        finally:
            conn.close()


def seed_rules() -> None:
    conn = analysis_conn()
    try:
        cursor = conn.cursor()
        for code, name, description in SEED_RULES:
            cursor.execute(
                "INSERT OR IGNORE INTO analysis_rules (code, name, description) VALUES (%s, %s, %s)",
                (code, name, description),
            )
        conn.commit()
        print(f"  OK {len(SEED_RULES)} analysis rules")
    finally:
        conn.close()


def main() -> int:
    print("\n=== SQLite ===")
    create_sqlite()
    seed_rules()
    print("\n=== ClickHouse ===")
    create_clickhouse()
    print("\nDone.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
