"""Tells you what is missing before you find out halfway through a batch.

    py -m tools.check_env

It checks, in this order: the converters, ClickHouse, the schema, SQLite and the
current state of the databases. Exits with code 1 if anything essential is
missing.
"""

from __future__ import annotations

import sys
from pathlib import Path

from app import config

OK = "  ok  "
MISSING = "MISSING"
WARNING = "  warn"


def line(mark: str, text: str) -> None:
    print(f"[{mark}] {text}")


def check_converters() -> int:
    problems = 0
    trfp = Path(config.THERMO_RAW_PARSER_PATH)
    if trfp.exists():
        line(OK, f"ThermoRawFileParser: {trfp}")
    else:
        problems += 1
        line(MISSING, f"ThermoRawFileParser is not at {trfp}")
        print("        It ships with the repository. Without it a Thermo .raw cannot")
        print("        be converted, and msconvert will not do: it keeps the exception")
        print("        data, which turns into false positives.")

    msconvert = Path(config.MSCONVERT_PATH)
    if msconvert.exists():
        line(OK, f"msconvert: {msconvert}")
    else:
        line(WARNING, f"msconvert is not at {msconvert}")
        print("        Only needed for non-Thermo vendors (Waters, Agilent, Sciex...).")
    return problems


def check_clickhouse() -> int:
    from app.db.clickhouse import server_ch

    try:
        ch = server_ch()
        version = ch.query("SELECT version()").result_rows[0][0]
        line(OK, f"ClickHouse {version} at {config.CH_HOST}:{config.CH_PORT}")
    except Exception as exc:  # noqa: BLE001
        line(MISSING, f"ClickHouse unreachable at {config.CH_HOST}:{config.CH_PORT} - {exc}")
        print("        Start it with: docker compose up -d")
        return 1

    missing = []
    for database, tables in (
        (config.CH_LAKE_DB, ("scans", "spectra_points", "chromatogram_points")),
        (config.CH_ANALYSIS_DB, ("processing_sampleparameter",)),
    ):
        existing = {
            r[0] for r in ch.query(
                "SELECT name FROM system.tables WHERE database = %(db)s",
                parameters={"db": database},
            ).result_rows
        }
        for table in tables:
            if table not in existing:
                missing.append(f"{database}.{table}")

    if missing:
        line(MISSING, f"missing tables: {', '.join(missing)}")
        print("        Create them with: py -m app.db.schema")
        ch.close()
        return 1

    line(OK, f"tables in {config.CH_LAKE_DB} and {config.CH_ANALYSIS_DB} created")

    # Disk usage - the spectra table grows fast and the Docker volume usually
    # lives on the system drive.
    rows = ch.query(
        """
        SELECT database, table, formatReadableSize(sum(bytes_on_disk)), sum(rows)
        FROM system.parts WHERE active AND database IN (%(a)s, %(b)s)
        GROUP BY database, table ORDER BY sum(bytes_on_disk) DESC
        """,
        parameters={"a": config.CH_LAKE_DB, "b": config.CH_ANALYSIS_DB},
    ).result_rows
    for database, table, size, n in rows:
        line(OK, f"{database}.{table}: {size} across {n:,} rows")
    ch.close()
    return 0


def check_sqlite() -> int:
    from app.db.sqlite import analysis_conn, lake_conn

    for label, file, factory, queries in (
        ("lake", config.LAKE_DB, lake_conn, {
            "samples": "SELECT COUNT(*) FROM samples",
            "ready": "SELECT COUNT(*) FROM samples WHERE status = 'ready'",
        }),
        ("analysis", config.ANALYSIS_DB, analysis_conn, {
            "methods": "SELECT COUNT(*) FROM methods WHERE deleted_at IS NULL",
            "substances": "SELECT COUNT(*) FROM method_substances WHERE deleted_at IS NULL",
            "tasks": "SELECT COUNT(*) FROM sample_processing_tasks",
            "results": "SELECT COUNT(*) FROM sample_results",
        }),
    ):
        if not Path(file).exists():
            line(MISSING, f"the {label} SQLite file does not exist ({file}) - run py -m app.db.schema")
            return 1
        conn = factory()
        try:
            cursor = conn.cursor()
            parts = []
            for name, sql in queries.items():
                cursor.execute(sql)
                parts.append(f"{name}={cursor.fetchone()[0]}")
            cursor.close()
            line(OK, f"SQLite {label}: " + ", ".join(parts))
        except Exception as exc:  # noqa: BLE001
            line(MISSING, f"the {label} SQLite file is unreadable: {exc}")
            return 1
        finally:
            conn.close()
    return 0


def main() -> int:
    print(f"\nProject: {config.ROOT}\n")
    problems = check_converters()
    print()
    problems += check_clickhouse()
    print()
    problems += check_sqlite()
    print()
    if problems:
        print(f"{problems} problem(s) - fix them before uploading a sample.\n")
        return 1
    print("Environment ready.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
