"""The main screen: pick N samples plus the 3 controls and run the pipeline.

Every sample is analysed against all three controls at once:

    CQCL        reference - where and how the substance shows up
    CQCL_reinj  second reference, used when the CQCL comes out NEGATIVE
    CQN         baseline, subtracted from the others before comparing

That is why all three are required: with no CQN there is no matrix correction,
and with no CQCL there is nothing to compare against (the substance would come
out INVALID, which does not even produce a result row).
"""

from __future__ import annotations

import json

from flask import Blueprint, flash, jsonify, redirect, render_template, request, url_for

from app.db.sqlite import analysis_conn, lake_conn
from app.services import runner
from app.services.build_method_snapshot import build_method_snapshot

bp = Blueprint("tasks", __name__, url_prefix="/processing")


@bp.route("/")
def index():
    conn = analysis_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute(
            """
            SELECT t.*,
                   (SELECT COUNT(*) FROM sample_results r
                     WHERE r.sample_processing_task_id = t.id) AS n_results,
                   (SELECT COUNT(*) FROM sample_results r
                     WHERE r.sample_processing_task_id = t.id
                       AND r.result LIKE 'SUSPECT%') AS n_suspects
            FROM sample_processing_tasks t
            ORDER BY t.id DESC
            """
        )
        tasks = cursor.fetchall()
        cursor.close()
    finally:
        conn.close()
    return render_template("tasks/index.html", tasks=tasks)


@bp.route("/new")
def new():
    lake = lake_conn()
    analysis = analysis_conn()
    try:
        lake_cursor = lake.cursor(dictionary=True)
        lake_cursor.execute(
            """
            SELECT s.*, f.path AS folder_path
            FROM samples s
            JOIN folders f ON f.id = s.folder_id
            WHERE f.deleted_at IS NULL
            ORDER BY s.batch DESC, s.is_control DESC, s.name
            """
        )
        samples = lake_cursor.fetchall()
        lake_cursor.close()

        analysis_cursor = analysis.cursor(dictionary=True)
        analysis_cursor.execute(
            """
            SELECT m.*,
                   -- Same total as the method screen: individual substances plus
                   -- the ones analysed through a group.
                   (SELECT COUNT(*) FROM method_substances ms
                     WHERE ms.method_id = m.id AND ms.deleted_at IS NULL
                       AND ms.type <> 'internal_standard') AS n_substances
            FROM methods m
            WHERE m.deleted_at IS NULL
            ORDER BY m.name
            """
        )
        methods = analysis_cursor.fetchall()
        analysis_cursor.close()
    finally:
        lake.close()
        analysis.close()

    # Grouped by batch: the user thinks in batches, not in a flat list.
    batches: dict[str, list] = {}
    for sample in samples:
        batches.setdefault(sample["batch"], []).append(sample)

    return render_template("tasks/new.html", batches=batches, methods=methods)


@bp.route("/new", methods=["POST"])
def create():
    method_id = request.form.get("method_id", type=int)
    sample_ids = request.form.getlist("samples", type=int)
    cqcl_id = request.form.get("cqcl_id", type=int)
    cqcl_reinj_id = request.form.get("cqcl_reinj_id", type=int)
    cqn_id = request.form.get("cqn_id", type=int)
    validate_is = 1 if request.form.get("validate_is") else 0
    batch_label = (request.form.get("batch") or "").strip()

    errors = []
    if not method_id:
        errors.append("pick a method")
    if not sample_ids:
        errors.append("select at least one sample")
    if not (cqcl_id and cqcl_reinj_id and cqn_id):
        errors.append("all three controls are required (CQCL, CQCL_reinj and CQN)")

    controls = {cqcl_id, cqcl_reinj_id, cqn_id}
    if cqn_id and cqn_id in {cqcl_id, cqcl_reinj_id}:
        errors.append("the CQN cannot be the same file as a CQCL - it is the baseline")
    if controls & set(sample_ids):
        errors.append("a selected sample is also marked as a control")

    if errors:
        for error in errors:
            flash(error, "error")
        return redirect(url_for("tasks.new"))

    # Everyone involved has to be ingested: with no points in ClickHouse there
    # is no XIC, and the pipeline would return INVALID for everything.
    lake = lake_conn()
    try:
        everyone = list(sample_ids) + [cqcl_id, cqcl_reinj_id, cqn_id]
        placeholders = ", ".join("%s" for _ in everyone)
        cursor = lake.cursor(dictionary=True)
        cursor.execute(f"SELECT id, name, batch, status FROM samples WHERE id IN ({placeholders})", tuple(everyone))
        found = {int(row["id"]): row for row in cursor.fetchall()}
        cursor.close()
    finally:
        lake.close()

    not_ready = [row["name"] for row in found.values() if row["status"] != "ready"]
    if not_ready:
        flash(f"Not ingested yet: {', '.join(not_ready)}", "error")
        return redirect(url_for("tasks.new"))

    if cqcl_reinj_id == cqcl_id:
        flash(
            "CQCL_reinj is the same as the CQCL: the fallback will re-evaluate exactly "
            "the same data, so there is no real second chance.",
            "warning",
        )

    batch = batch_label or found.get(sample_ids[0], {}).get("batch") or "ADHOC"

    analysis = analysis_conn()
    try:
        # Snapshot: a frozen copy of the method. Processing reads the method
        # live; the snapshot is there so the result stays interpretable after
        # someone edits a threshold.
        snapshot = build_method_snapshot(method_id)

        cursor = analysis.cursor()
        cursor.execute(
            """
            INSERT INTO sample_processing_tasks
                (samples, method_id, batch, computer, screening, method_snapshot,
                 status, cqcl_id, cqcl_reinj_id, cqn_id, with_controls,
                 validate_is, is_auto_processed)
            VALUES (%s, %s, %s, 'upload', %s, %s, 'queued', %s, %s, %s, 1, %s, 0)
            """,
            (
                json.dumps(sample_ids),
                method_id,
                batch,
                None,
                json.dumps(snapshot, ensure_ascii=False),
                cqcl_id,
                cqcl_reinj_id,
                cqn_id,
                validate_is,
            ),
        )
        task_id = int(cursor.lastrowid)
        cursor.close()
        analysis.commit()
    finally:
        analysis.close()

    runner.submit(runner.run_task, task_id)
    flash(f"Task {task_id} queued: {len(sample_ids)} sample(s).", "ok")
    return redirect(url_for("tasks.detail", task_id=task_id))


@bp.route("/<int:task_id>")
def detail(task_id: int):
    analysis = analysis_conn()
    lake = lake_conn()
    try:
        cursor = analysis.cursor(dictionary=True)
        cursor.execute("SELECT * FROM sample_processing_tasks WHERE id = %s", (task_id,))
        task = cursor.fetchone()
        if not task:
            flash("Task not found.", "error")
            return redirect(url_for("tasks.index"))

        cursor.execute("SELECT name FROM methods WHERE id = %s", (task["method_id"],))
        row = cursor.fetchone()
        method_name = row["name"] if row else f"#{task['method_id']}"
        cursor.close()

        ids = json.loads(task["samples"]) + [task["cqcl_id"], task["cqcl_reinj_id"], task["cqn_id"]]
        placeholders = ", ".join("%s" for _ in ids)
        lake_cursor = lake.cursor(dictionary=True)
        lake_cursor.execute(f"SELECT id, name, batch FROM samples WHERE id IN ({placeholders})", tuple(ids))
        names = {int(row["id"]): row["name"] for row in lake_cursor.fetchall()}
        lake_cursor.close()
    finally:
        analysis.close()
        lake.close()

    return render_template(
        "tasks/detail.html",
        task=task,
        method_name=method_name,
        names=names,
        samples=json.loads(task["samples"]),
        outcome=json.loads(task["result"]) if task["result"] else None,
        audit=json.loads(task["internal_standards_audit"]) if task["internal_standards_audit"] else None,
    )


@bp.route("/<int:task_id>/api")
def api(task_id: int):
    analysis = analysis_conn()
    try:
        cursor = analysis.cursor(dictionary=True)
        cursor.execute("SELECT id, status, result FROM sample_processing_tasks WHERE id = %s", (task_id,))
        task = cursor.fetchone()
        cursor.execute(
            "SELECT COUNT(*) AS n FROM sample_results WHERE sample_processing_task_id = %s",
            (task_id,),
        )
        n_results = cursor.fetchone()["n"]
        cursor.close()
    finally:
        analysis.close()
    return jsonify({
        "status": task["status"] if task else "?",
        "result": json.loads(task["result"]) if task and task["result"] else None,
        "n_results": n_results,
        "queued": runner.queue_size(),
        "log": runner.recent_log(30),
    })
