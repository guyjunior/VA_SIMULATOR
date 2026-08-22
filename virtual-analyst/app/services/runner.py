"""Heavy work, outside the request.

Converting a .raw and ingesting its spectra takes minutes; processing a batch
takes longer. None of it can happen inside a Flask request - the browser would
time out and the user would see no progress at all.

DESIGN: a single queue with ONE worker.

Parallelism here would be counterproductive: the converter already uses the
whole CPU, and two simultaneous ingests fight over ClickHouse memory (which
answers with MEMORY_LIMIT_EXCEEDED, not with slowness). Serial is faster overall
and far more predictable.

STATE LIVES IN THE DATABASE, not in the worker's memory: `samples.status` and
`sample_processing_tasks.status`. The screen reads from there. If Flask dies
mid-run, whatever was running stays marked as it was - visible, not a ghost.
"""

from __future__ import annotations

import json
import queue
import threading
import time
import traceback
from datetime import datetime
from pathlib import Path

from app import config
from app.db.clickhouse import analysis_ch, lake_ch
from app.db.sqlite import analysis_conn, lake_conn
from app.services import convert, ingest

_queue: "queue.Queue[tuple]" = queue.Queue()
_worker: threading.Thread | None = None
_lock = threading.Lock()

# The last log lines of each job, so the screen can show them without opening a
# file.
_LOG_MAX = 400
_log_lines: list[str] = []
_log_lock = threading.Lock()


def log(message: str) -> None:
    line = f"{datetime.now():%H:%M:%S} {message}"
    with _log_lock:
        _log_lines.append(line)
        del _log_lines[:-_LOG_MAX]
    print(line, flush=True)
    try:
        with open(config.LOG_DIR / "runner.log", "a", encoding="utf-8") as handle:
            handle.write(line + "\n")
    except OSError:
        pass


def recent_log(n: int = 60) -> list[str]:
    with _log_lock:
        return _log_lines[-n:]


def queue_size() -> int:
    return _queue.qsize()


def _run_worker() -> None:
    while True:
        fn, args, kwargs = _queue.get()
        try:
            fn(*args, **kwargs)
        except Exception as exc:  # noqa: BLE001 - the worker must never die
            log(f"Unhandled error in the worker: {exc}")
            log(traceback.format_exc())
        finally:
            _queue.task_done()


def submit(fn, *args, **kwargs) -> None:
    """Enqueue. Starts the worker on the first call."""
    global _worker
    with _lock:
        if _worker is None or not _worker.is_alive():
            _worker = threading.Thread(target=_run_worker, name="va-worker", daemon=True)
            _worker.start()
    _queue.put((fn, args, kwargs))


# ---------------------------------------------------------------------------
# Job 1 - convert and ingest one sample
# ---------------------------------------------------------------------------

def _set_sample(conn, sample_id: int, **fields) -> None:
    assignments = ", ".join(f"{key} = %s" for key in fields)
    cursor = conn.cursor()
    cursor.execute(f"UPDATE samples SET {assignments} WHERE id = %s", (*fields.values(), sample_id))
    cursor.close()
    conn.commit()


def process_sample_upload(sample_id: int, converter_choice: str) -> None:
    """Convert the sample's .raw file and ingest the mzML.

    Every sample is its own job: if one fails (file being acquired, empty .raw),
    the rest of the batch carries on. The reason for the failure is kept in
    `samples.error` and shows up on screen - it does not vanish into a log.
    """
    conn = lake_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute("SELECT * FROM samples WHERE id = %s", (sample_id,))
        sample = cursor.fetchone()
        cursor.close()
        if not sample:
            log(f"sample {sample_id} disappeared before processing")
            return

        name = sample["name"]
        raw_path = Path(sample["raw_path"])
        started = time.time()

        # --- conversion ---
        _set_sample(conn, sample_id, status="converting", error=None)
        vendor = convert.detect_vendor(raw_path)
        converter = convert.pick_converter(vendor, converter_choice)
        warning = convert.converter_warning(vendor, converter)
        if warning:
            log(f"[{name}] WARNING: {warning}")
        log(f"[{name}] converting ({vendor} -> {converter})")

        try:
            mzml_path, converter_used = convert.convert_to_mzml(
                raw_path, config.MZML_DIR, converter_choice
            )
        except convert.ConversionFailure as failure:
            reason = convert.SKIP_REASON.get(failure.kind)
            if reason:
                # Not a system error: the file simply has nothing to convert.
                # Marking it 'skipped' keeps the screen from turning red over
                # something that will never work.
                _set_sample(conn, sample_id, status="skipped", error=reason)
                log(f"[{name}] skipped: {reason}")
                return
            raise

        _set_sample(conn, sample_id, status="ingesting", mzml_path=str(mzml_path), converter=converter_used)
        log(f"[{name}] mzML ready in {time.time() - started:.0f}s - ingesting")

        # --- ingestion ---
        summary = ingest.ingest_mzml(
            conn, mzml_path,
            sample_id=sample_id,
            year=int(sample["year"]),
            screening=sample["screening"],
        )
        _set_sample(conn, sample_id, status="ready", error=None)
        log(
            f"[{name}] ready in {time.time() - started:.0f}s - "
            f"{summary['scans']} scans, {summary['spectra_points']} points"
        )
    except Exception as exc:  # noqa: BLE001 - the error has to reach the screen
        log(f"[sample {sample_id}] FAILED: {exc}")
        log(traceback.format_exc())
        try:
            _set_sample(conn, sample_id, status="error", error=str(exc)[:2000])
        except Exception:
            pass
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Job 2 - process one task (N samples against the 3 controls)
# ---------------------------------------------------------------------------

def run_task(task_id: int) -> None:
    """Run the analytical pipeline over every sample in the task.

    It mirrors the production orchestrator, including in what looks like detail:

      * ONE analysis connection and ONE ClickHouse client for the whole batch,
        handed down to each sample. Opening them per sample exhausts ephemeral
        ports on Windows (WinError 10048/10055) on a large batch.
      * A failure in ONE sample does not take down the batch: it becomes
        `error: ...` in that sample's result and the loop carries on.
      * The status returned per sample is stored as-is ("processed",
        "blocked_by_is_gate", "method_not_found"). Treating everything as
        "processed" has already hidden a sample aborted by the IS gate: task
        'done', zero results, and the field claiming success.
    """
    from app.services.virtual_analyst import process_sample_from_db

    lake_db = lake_conn()
    lake_click = lake_ch()
    analysis_db = analysis_conn()
    analysis_click = analysis_ch()

    try:
        cursor = analysis_db.cursor(dictionary=True)
        cursor.execute("SELECT * FROM sample_processing_tasks WHERE id = %s", (task_id,))
        task = cursor.fetchone()
        cursor.close()
        if not task:
            log(f"task {task_id} not found")
            return

        sample_ids = json.loads(task["samples"])
        method_id = int(task["method_id"])
        batch = task["batch"]

        cursor = analysis_db.cursor()
        cursor.execute(
            "UPDATE sample_processing_tasks SET status = 'processing', updated_at = NOW() WHERE id = %s",
            (task_id,),
        )
        cursor.close()
        analysis_db.commit()

        # Sample metadata (name, year, flags) comes from the lake.
        placeholders = ", ".join("%s" for _ in sample_ids)
        cursor = lake_db.cursor(dictionary=True)
        cursor.execute(f"SELECT * FROM samples WHERE id IN ({placeholders})", tuple(sample_ids))
        by_id = {int(row["id"]): row for row in cursor.fetchall()}
        cursor.close()

        log(f"task {task_id}: batch {batch}, {len(sample_ids)} sample(s), method {method_id}")
        started = time.time()
        result: dict[str, str] = {}

        for sample_id in sample_ids:
            sample = by_id.get(int(sample_id))
            if not sample:
                result[str(sample_id)] = "error: sample not found"
                continue
            name = sample["name"]
            try:
                # NULL -> True on both flags. Safe default: when in doubt,
                # analyse EVERYTHING; only an explicit False narrows the scope.
                competition = sample.get("competition")
                in_competition = True if competition is None else bool(competition)
                beta_block = sample.get("beta_block")
                is_beta_blocker = True if beta_block is None else bool(beta_block)

                status = process_sample_from_db(
                    src_mdb_conn=lake_db,
                    src_ch_client=lake_click,
                    sample_id=int(sample_id),
                    sample_name=name,
                    cqcl_id=int(task["cqcl_id"]),
                    cqcl_reinj_id=int(task["cqcl_reinj_id"]),
                    cqn_id=int(task["cqn_id"]),
                    batch_name=batch,
                    method_id=method_id,
                    task_id=task_id,
                    mass_tolerance_ppm=config.MASS_TOLERANCE_PPM,
                    in_competition=in_competition,
                    is_beta_blocker=is_beta_blocker,
                    validate_internal_standards=bool(task["validate_is"]),
                    year=int(sample["year"]),
                    analysis_mdb=analysis_db,
                    analysis_ch_client=analysis_click,
                )
                result[name] = status or "processed"
                log(f"task {task_id}: {name} -> {result[name]}")
            except Exception as exc:  # noqa: BLE001 - one sample must not kill the batch
                log(f"task {task_id}: error on {name}: {exc}")
                log(traceback.format_exc())
                result[name] = f"error: {exc}"

        cursor = analysis_db.cursor()
        cursor.execute(
            "UPDATE sample_processing_tasks SET status = 'done', result = %s, updated_at = NOW() WHERE id = %s",
            (json.dumps(result, ensure_ascii=False), task_id),
        )
        cursor.close()
        analysis_db.commit()
        log(f"task {task_id} finished in {(time.time() - started) / 60:.1f} min")

    except Exception as exc:  # noqa: BLE001
        log(f"task {task_id} FAILED: {exc}")
        log(traceback.format_exc())
        try:
            cursor = analysis_db.cursor()
            cursor.execute(
                "UPDATE sample_processing_tasks SET status = 'error', result = %s, updated_at = NOW() WHERE id = %s",
                (json.dumps({"error": str(exc)}), task_id),
            )
            cursor.close()
            analysis_db.commit()
        except Exception:
            pass
    finally:
        for resource in (lake_click, analysis_click):
            try:
                resource.close()
            except Exception:
                pass
        for resource in (lake_db, analysis_db):
            try:
                resource.close()
            except Exception:
                pass
