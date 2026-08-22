r"""Bring a fresh installation up with a working example, in one command.

    venv\Scripts\python -m tools.seed_example

It creates the schema, loads the example method and ingests the four example
samples, leaving the app ready for the user to create the first analysis. It does
NOT run the analysis - that is the part worth doing by hand the first time.

Everything is idempotent: run it again and it re-does only what is missing. Use
--force to rebuild the method and re-ingest the samples from scratch.

The four .raw files ship with the repository, in `examples/raw/`. They are all
CONTROL injections - no athlete sample is published. `Sample.raw` is a further
CQCL injection standing in for the sample under analysis; its internal header
still says CQCL, and the README explains why that is the honest choice.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from app import config
from app.db.sqlite import analysis_conn, lake_conn

METHOD_NAME = "Method Example"
METHOD_JSON = "data/method_article.json"
BATCH = "EXAMPLE"
RAW_DIR = "examples/raw"

# The example set: one sample and the three controls it is judged against.
EXPECTED = ("Sample.raw", "CQCL.raw", "CQCL_reinj.raw", "CQN.raw")


def step(number: int, text: str) -> None:
    print(f"\n[{number}/4] {text}")


def create_schema() -> None:
    from app.db.schema import create_clickhouse, create_sqlite, seed_rules

    create_sqlite()
    seed_rules()
    create_clickhouse()


def load_method(force: bool) -> bool:
    """Returns True when the method is in place."""
    conn = analysis_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute("SELECT id FROM methods WHERE name = %s", (METHOD_NAME,))
        existing = cursor.fetchone()
        cursor.close()
    finally:
        conn.close()

    if existing and not force:
        print(f"  '{METHOD_NAME}' is already loaded (id {existing['id']}) - skipping")
        return True

    json_path = Path(METHOD_JSON)
    if not json_path.exists():
        print(f"  {json_path} not found. It is the exported method; see the README.")
        return False

    argv = [sys.argv[0], METHOD_JSON, "--name", METHOD_NAME]
    if existing:
        argv.append("--replace")

    from tools import import_method

    saved, sys.argv = sys.argv, argv
    try:
        return import_method.main() == 0
    finally:
        sys.argv = saved


def register_samples(raw_dir: Path, batch: str, force: bool) -> list[int]:
    """Create the sample rows, pointing at the .raw files where they already are.

    They are referenced in place rather than copied into data/uploads/: the files
    are already inside the project, and a second 324 MB copy would buy nothing.
    """
    from app.services.ingest import is_control_sample, year_from_batch

    missing = [name for name in EXPECTED if not (raw_dir / name).exists()]
    if missing:
        print(f"  missing from {raw_dir}: {', '.join(missing)}")
        return []

    year = year_from_batch(batch)
    conn = lake_conn()
    ids: list[int] = []
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute("SELECT id, name, status FROM samples WHERE batch = %s", (batch,))
        present = {row["name"]: row for row in cursor.fetchall()}

        if present and not force:
            ready = [r for r in present.values() if r["status"] == "ready"]
            print(f"  batch '{batch}' already registered "
                  f"({len(ready)}/{len(present)} ingested) - skipping")
            return [int(r["id"]) for r in present.values() if r["status"] != "ready"]

        if present and force:
            from app.db.clickhouse import lake_ch
            from app.services.ingest import purge_sample_points

            ch = lake_ch()
            try:
                for row in present.values():
                    purge_sample_points(ch, int(row["id"]))
            finally:
                ch.close()
            writer = conn.cursor()
            for row in present.values():
                writer.execute("DELETE FROM sample_channels WHERE sample_id = %s", (row["id"],))
                writer.execute("DELETE FROM samples WHERE id = %s", (row["id"],))
            writer.close()
            conn.commit()
            print(f"  removed the previous '{batch}' ({len(present)} samples)")

        writer = conn.cursor()
        writer.execute(
            "INSERT INTO folders (path, screening, computer, chromatographic_technique) "
            "VALUES (%s, %s, 'example', 'LC-MS')",
            (str(raw_dir.resolve()), "TVII"),
        )
        folder_id = int(writer.lastrowid)

        for name in EXPECTED:
            stem = Path(name).stem
            writer.execute(
                """
                INSERT INTO samples
                    (folder_id, name, batch, year, screening, computer,
                     chromatography_technique, is_control, status, raw_path)
                VALUES (%s, %s, %s, %s, 'TVII', 'example', 'LC-MS', %s, 'pending', %s)
                """,
                (folder_id, stem, batch, year,
                 1 if is_control_sample(stem) else 0,
                 str((raw_dir / name).resolve())),
            )
            ids.append(int(writer.lastrowid))
        writer.close()
        conn.commit()
        print(f"  registered {len(ids)} samples in batch '{batch}' (year {year})")
    finally:
        conn.close()
    return ids


def ingest(sample_ids: list[int]) -> None:
    """Convert and ingest, in the foreground.

    The web app does this on a background worker so the browser is not left
    hanging; here running in the foreground is the point - the operator watches
    it happen and sees any failure immediately.
    """
    from app.services.runner import process_sample_upload

    for position, sample_id in enumerate(sample_ids, start=1):
        print(f"  ({position}/{len(sample_ids)})", end=" ", flush=True)
        process_sample_upload(sample_id, "auto")


def report() -> None:
    lake = lake_conn()
    analysis = analysis_conn()
    try:
        cursor = lake.cursor(dictionary=True)
        cursor.execute("SELECT name, status, error FROM samples WHERE batch = %s ORDER BY id", (BATCH,))
        samples = cursor.fetchall()
        cursor.close()

        cursor = analysis.cursor(dictionary=True)
        cursor.execute(
            """
            SELECT m.name,
                   (SELECT COUNT(*) FROM method_substances ms
                     WHERE ms.method_id = m.id AND ms.deleted_at IS NULL
                       AND ms.type <> 'internal_standard') AS n_substances
            FROM methods m WHERE m.name = %s
            """,
            (METHOD_NAME,),
        )
        method = cursor.fetchone()
        cursor.close()
    finally:
        lake.close()
        analysis.close()

    print("\n" + "=" * 64)
    if method:
        print(f"method   {method['name']} - {method['n_substances']} substances")
    for sample in samples:
        mark = "ok " if sample["status"] == "ready" else "!! "
        detail = f"  {sample['error']}" if sample["error"] else ""
        print(f"{mark}      {sample['name']:14} {sample['status']}{detail}")
    print("=" * 64)

    if all(s["status"] == "ready" for s in samples) and method:
        print("\nReady. Start the app and create the first analysis:")
        print("    venv\\Scripts\\python run.py")
        print("\nOn the New analysis screen, pick 'Method Example', tick 'Sample'")
        print("as the sample to analyse, and set CQCL / CQCL_reinj / CQN as the")
        print("controls (the 'Guess from names' button fills them in).")
        print("\nAll four files are control injections - 'Sample' is a CQCL standing")
        print("in for a sample, so expect findings rather than a clean screen.")
    else:
        print("\nSomething is missing above. `py -m tools.check_env` usually says why.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default=RAW_DIR)
    parser.add_argument("--batch", default=BATCH)
    parser.add_argument("--force", action="store_true",
                        help="rebuild the method and re-ingest the samples")
    args = parser.parse_args()

    started = time.time()
    print(f"Seeding the example in {config.ROOT}")

    step(1, "schema (SQLite + ClickHouse)")
    try:
        create_schema()
    except Exception as exc:  # noqa: BLE001 - the message is the product here
        print(f"  failed: {exc}")
        print("  Is ClickHouse up? docker compose up -d")
        return 1

    step(2, f"method '{METHOD_NAME}'")
    if not load_method(args.force):
        return 1

    step(3, f"example samples from {args.raw_dir}")
    pending = register_samples(Path(args.raw_dir), args.batch, args.force)

    step(4, "conversion and ingestion" if pending else "conversion and ingestion (nothing to do)")
    if pending:
        ingest(pending)

    report()
    print(f"\ntotal: {(time.time() - started) / 60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
