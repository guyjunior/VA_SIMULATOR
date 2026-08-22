"""Uploading .raw files, converting and ingesting them.

The upload creates the `samples` rows with status 'pending' and returns the page
immediately; converting and ingesting happen on the background worker. A batch
of 20 samples takes tens of minutes - holding the request until the end would
time out the browser and show no progress at all.
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

from flask import Blueprint, flash, jsonify, redirect, render_template, request, url_for
from werkzeug.utils import secure_filename

from app import config
from app.db.clickhouse import lake_ch
from app.db.sqlite import lake_conn
from app.services import convert, ingest, runner

bp = Blueprint("samples", __name__, url_prefix="/samples")

# Extensions the pipeline accepts. `.raw` is the normal case; mzML goes straight
# through (skipping conversion) and helps compare converters outside the app.
EXTENSIONS = {".raw", ".mzml", ".mzxml", ".mgf", ".wiff", ".wiff2", ".lcd"}


def _slug_batch(value: str) -> str:
    """A batch name safe to use as a folder name and as a label."""
    value = (value or "").strip()
    value = re.sub(r"[^A-Za-z0-9_.\-]+", "_", value)
    return value or f"BATCH_{datetime.now():%Y%m%d_%H%M%S}"


@bp.route("/")
def index():
    conn = lake_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute(
            """
            SELECT f.id, f.path, f.screening, f.created_at,
                   COUNT(s.id)                                         AS total,
                   SUM(CASE WHEN s.status = 'ready' THEN 1 ELSE 0 END) AS ready,
                   SUM(CASE WHEN s.status = 'error' THEN 1 ELSE 0 END) AS failed
            FROM folders f
            LEFT JOIN samples s ON s.folder_id = f.id
            WHERE f.deleted_at IS NULL
            GROUP BY f.id
            ORDER BY f.id DESC
            """
        )
        batches = cursor.fetchall()

        cursor.execute(
            """
            SELECT s.*, f.path AS folder_path
            FROM samples s
            JOIN folders f ON f.id = s.folder_id
            WHERE f.deleted_at IS NULL
            ORDER BY s.folder_id DESC, s.name
            """
        )
        samples = cursor.fetchall()
        cursor.close()
    finally:
        conn.close()

    return render_template(
        "samples/index.html",
        batches=batches,
        samples=samples,
        converters=convert.CONVERTERS,
        queued=runner.queue_size(),
    )


@bp.route("/upload", methods=["POST"])
def upload():
    batch = _slug_batch(request.form.get("batch", ""))
    screening = (request.form.get("screening") or "").strip() or None
    choice = request.form.get("converter", convert.AUTO)
    if choice not in convert.CONVERTERS:
        choice = convert.AUTO

    files = [f for f in request.files.getlist("files") if f and f.filename]
    if not files:
        flash("No file selected.", "error")
        return redirect(url_for("samples.index"))

    target_dir = config.UPLOAD_DIR / batch
    target_dir.mkdir(parents=True, exist_ok=True)
    year = ingest.year_from_batch(batch)

    conn = lake_conn()
    created, ignored = [], []
    try:
        cursor = conn.cursor()
        cursor.execute(
            "INSERT INTO folders (path, screening, computer, chromatographic_technique) VALUES (%s, %s, %s, %s)",
            (str(target_dir), screening, "upload", "LC-MS"),
        )
        folder_id = int(cursor.lastrowid)

        for uploaded in files:
            filename = secure_filename(Path(uploaded.filename).name)
            if Path(filename).suffix.lower() not in EXTENSIONS:
                ignored.append(filename)
                continue

            path = target_dir / filename
            uploaded.save(path)
            name = Path(filename).stem

            cursor.execute(
                """
                INSERT INTO samples
                    (folder_id, name, batch, year, screening, computer,
                     chromatography_technique, is_control, status, raw_path)
                VALUES (%s, %s, %s, %s, %s, 'upload', 'LC-MS', %s, 'pending', %s)
                """,
                (folder_id, name, batch, year, screening,
                 1 if ingest.is_control_sample(name) else 0, str(path)),
            )
            created.append(int(cursor.lastrowid))

        cursor.close()
        conn.commit()
    finally:
        conn.close()

    for sample_id in created:
        runner.submit(runner.process_sample_upload, sample_id, choice)

    if ignored:
        flash(f"Ignored (unsupported extension): {', '.join(ignored)}", "warning")
    flash(f"{len(created)} file(s) queued for batch {batch}.", "ok")
    return redirect(url_for("samples.index"))


@bp.route("/<int:sample_id>/reprocess", methods=["POST"])
def reprocess(sample_id: int):
    """Convert and ingest again, possibly with a different converter.

    Ingestion deletes the sample's old points before writing (see
    `ingest.purge_sample_points`), so switching converters and reprocessing
    replaces the data - it does not add to it.
    """
    choice = request.form.get("converter", convert.AUTO)
    if choice not in convert.CONVERTERS:
        choice = convert.AUTO
    runner.submit(runner.process_sample_upload, sample_id, choice)
    flash("Sample queued again.", "ok")
    return redirect(url_for("samples.index"))


@bp.route("/<int:sample_id>/delete", methods=["POST"])
def delete(sample_id: int):
    conn = lake_conn()
    ch = lake_ch()
    try:
        ingest.purge_sample_points(ch, sample_id)
        cursor = conn.cursor()
        cursor.execute("DELETE FROM sample_channels WHERE sample_id = %s", (sample_id,))
        cursor.execute("DELETE FROM samples WHERE id = %s", (sample_id,))
        cursor.close()
        conn.commit()
    finally:
        ch.close()
        conn.close()
    flash("Sample removed (metadata and points).", "ok")
    return redirect(url_for("samples.index"))


@bp.route("/api/status")
def api_status():
    """Screen polling: the status of each sample plus the worker's recent log."""
    conn = lake_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute("SELECT id, name, status, converter, error FROM samples ORDER BY id")
        samples = cursor.fetchall()
        cursor.close()
    finally:
        conn.close()
    return jsonify({
        "samples": samples,
        "queued": runner.queue_size(),
        "log": runner.recent_log(30),
    })


@bp.route("/<int:sample_id>/scan-filters")
def scan_filters(sample_id: int):
    """The scan filters observed in this sample - used by the method screen."""
    conn = lake_conn()
    try:
        return jsonify(ingest.scan_filters_of_sample(conn, sample_id))
    finally:
        conn.close()
