"""Results: the per-sample table and the chromatogram behind each row.

ClickHouse stores the four traces RAW, exactly as the pipeline extracted them.
The CQN correction is applied HERE, at display time - the same decision the
production system makes, and for a good reason: the raw data stays auditable,
and the correction can change without reprocessing anything.
"""

from __future__ import annotations

import json

import numpy as np
from flask import Blueprint, flash, jsonify, redirect, render_template, request, url_for

from app.db.clickhouse import analysis_ch
from app.db.sqlite import analysis_conn

bp = Blueprint("results", __name__, url_prefix="/results")

# The four traces the pipeline stores. `cqcn` is the name the pipeline uses for
# the negative control in the points table.
TRACE_TYPES = ("sample", "cqcl", "cqcl_reinj", "cqcn")

VALIDATIONS = ("SUSPECT", "NEGATIVE")


@bp.route("/task/<int:task_id>")
def task(task_id: int):
    suspects_only = request.args.get("suspects") == "1"
    conn = analysis_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute("SELECT * FROM sample_processing_tasks WHERE id = %s", (task_id,))
        task_row = cursor.fetchone()
        if not task_row:
            flash("Task not found.", "error")
            return redirect(url_for("tasks.index"))

        clause = "AND result LIKE 'SUSPECT%'" if suspects_only else ""
        cursor.execute(
            f"""
            SELECT * FROM sample_results
            WHERE sample_processing_task_id = %s {clause}
            ORDER BY sample_name, substance
            """,
            (task_id,),
        )
        rows = cursor.fetchall()

        # Breakdown by classification - it gives you the batch at a glance,
        # before opening sample by sample.
        cursor.execute(
            """
            SELECT result, COUNT(*) AS n
            FROM sample_results
            WHERE sample_processing_task_id = %s
            GROUP BY result
            ORDER BY n DESC
            """,
            (task_id,),
        )
        breakdown = cursor.fetchall()
        cursor.close()
    finally:
        conn.close()

    by_sample: dict[str, list] = {}
    for row in rows:
        by_sample.setdefault(row["sample_name"], []).append(row)

    return render_template(
        "results/task.html",
        task=task_row,
        by_sample=by_sample,
        breakdown=breakdown,
        suspects_only=suspects_only,
        total=len(rows),
        validations=VALIDATIONS,
    )


def _points(ch, sample_result_id: int) -> dict[str, dict]:
    """Read the four traces of one result."""
    sql = """
        SELECT sample_type, retention_time, intensity
        FROM processing_sampleparameter
        WHERE sample_result_id = %(id)s
        ORDER BY sample_type, retention_time
    """
    rows = ch.query(sql, parameters={"id": int(sample_result_id)}).result_rows

    traces: dict[str, dict] = {t: {"rt": [], "intensity": []} for t in TRACE_TYPES}
    for trace_type, rt, intensity in rows:
        if trace_type in traces:
            traces[trace_type]["rt"].append(float(rt))
            traces[trace_type]["intensity"].append(float(intensity))
    return traces


def _subtract_cqn(target: dict, cqn: dict) -> dict:
    """Subtract the CQN point by point, matching on the nearest RT.

    In practice both traces come off the same RT grid, but there is no
    guarantee: different runs have slightly different scan counts. Matching on
    the nearest neighbour (via binary search) is what the production system
    does - the only difference here is that numpy does it in O(n log n) instead
    of O(n^2).

    A negative result becomes zero: negative intensity does not exist, and
    letting it through would only ruin the chart's scale.
    """
    if not target["rt"] or not cqn["rt"]:
        return target

    target_rt = np.asarray(target["rt"], dtype=float)
    target_intensity = np.asarray(target["intensity"], dtype=float)
    cqn_rt = np.asarray(cqn["rt"], dtype=float)
    cqn_intensity = np.asarray(cqn["intensity"], dtype=float)

    order = np.argsort(cqn_rt)
    cqn_rt, cqn_intensity = cqn_rt[order], cqn_intensity[order]

    position = np.searchsorted(cqn_rt, target_rt)
    position = np.clip(position, 1, len(cqn_rt) - 1)
    left, right = cqn_rt[position - 1], cqn_rt[position]
    take_left = np.abs(target_rt - left) <= np.abs(right - target_rt)
    index = np.where(take_left, position - 1, position)

    corrected = np.maximum(target_intensity - cqn_intensity[index], 0.0)
    return {"rt": target_rt.tolist(), "intensity": corrected.tolist()}


@bp.route("/xic/<int:sample_result_id>")
def xic(sample_result_id: int):
    """The traces, ready to plot.

    `?raw=1` returns the data exactly as it sits in ClickHouse, with no CQN
    correction - that is what you look at when a result does not match what you
    expected.
    """
    raw = request.args.get("raw") == "1"
    ch = analysis_ch()
    try:
        traces = _points(ch, sample_result_id)
    finally:
        ch.close()

    if raw:
        return jsonify({"corrected": False, "traces": traces})

    cqn = traces["cqcn"]
    return jsonify({
        "corrected": True,
        "traces": {
            "sample": _subtract_cqn(traces["sample"], cqn),
            "cqcl": _subtract_cqn(traces["cqcl"], cqn),
            "cqcl_reinj": _subtract_cqn(traces["cqcl_reinj"], cqn),
            "cqcn": cqn,
        },
    })


@bp.route("/<int:result_id>/validate", methods=["POST"])
def validate(result_id: int):
    """Record the analyst's reading without erasing the pipeline's.

    `result` is what the machine decided; `result_validation` is what the person
    decided. Keeping both side by side is what makes it possible to measure
    where the method gets it wrong.
    """
    value = request.form.get("result_validation")
    if value not in VALIDATIONS:
        value = None

    conn = analysis_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute("SELECT sample_processing_task_id FROM sample_results WHERE id = %s", (result_id,))
        row = cursor.fetchone()
        cursor.execute(
            "UPDATE sample_results SET result_validation = %s, updated_at = NOW() WHERE id = %s",
            (value, result_id),
        )
        cursor.close()
        conn.commit()
    finally:
        conn.close()

    target = url_for("results.task", task_id=row["sample_processing_task_id"]) if row else url_for("tasks.index")
    return redirect(target + f"#r{result_id}")


@bp.route("/task/<int:task_id>/json")
def export(task_id: int):
    """Export the task results. Useful for diffing against another system."""
    conn = analysis_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute("SELECT * FROM sample_processing_tasks WHERE id = %s", (task_id,))
        task_row = cursor.fetchone()
        cursor.execute(
            "SELECT * FROM sample_results WHERE sample_processing_task_id = %s ORDER BY sample_name, substance",
            (task_id,),
        )
        rows = cursor.fetchall()
        cursor.close()
    finally:
        conn.close()

    return jsonify({
        "task": {
            "id": task_row["id"],
            "batch": task_row["batch"],
            "status": task_row["status"],
            "method_snapshot": json.loads(task_row["method_snapshot"]),
            "result": json.loads(task_row["result"]) if task_row["result"] else None,
        },
        "results": rows,
    })
