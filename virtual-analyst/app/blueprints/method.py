"""The method: substances, internal standards, groups and the IS gate.

The method is what the pipeline reads for every sample. Three fields are usually
the cause when "nothing came out":

  * `scan_filter` - it matches the text seen in the mzML STRING FOR STRING. One
    extra space and the substance finds no channel and comes out INVALID (which
    does not even produce a row). That is why the screen offers the filters
    actually observed in the samples already ingested, instead of a free-text
    field.
  * `start_time`/`end_time` - the RT window, in minutes. Outside it there is no
    XIC.
  * `analysis_rule` - with no rule, the substance falls back to the standard one.
"""

from __future__ import annotations

from flask import Blueprint, flash, jsonify, redirect, render_template, request, url_for

from app.db.sqlite import analysis_conn, lake_conn

bp = Blueprint("method", __name__, url_prefix="/method")

TYPES = ("substance", "internal_standard")


def _float(name: str, *, required: bool = False) -> float | None:
    """Read a numeric field from the form. Empty becomes NULL, not zero.

    The difference matters: a NULL `cut_off` means "substance without a cut-off"
    (the pipeline skips that step); zero would mean a cut-off of zero, and
    everything would pass.
    """
    value = (request.form.get(name) or "").strip().replace(",", ".")
    if value == "":
        if required:
            raise ValueError(f"required field: {name}")
        return None
    return float(value)


def _int(name: str, default: int = 0) -> int:
    value = (request.form.get(name) or "").strip()
    return int(value) if value else default


def _scan_filter_id(conn, text: str | None) -> int | None:
    """Resolve (or create) the scan filter in the ANALYSIS database by its text."""
    text = (text or "").strip()
    if not text:
        return None
    cursor = conn.cursor()
    cursor.execute("INSERT OR IGNORE INTO scan_filters (name) VALUES (%s)", (text,))
    cursor.execute("SELECT id FROM scan_filters WHERE name = %s", (text,))
    filter_id = int(cursor.fetchone()[0])
    cursor.close()
    return filter_id


@bp.route("/")
def index():
    conn = analysis_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute(
            """
            SELECT m.*,
                   -- Everything the method looks for, which is NOT the same as
                   -- what the pipeline's individual-substance query returns: a
                   -- substance typed 'substance_group' is analysed through its
                   -- group, so counting only 'substance' hides it and the total
                   -- on screen comes out short.
                   (SELECT COUNT(*) FROM method_substances ms
                     WHERE ms.method_id = m.id AND ms.deleted_at IS NULL
                       AND ms.type <> 'internal_standard') AS n_substances,
                   (SELECT COUNT(*) FROM method_substances ms
                     WHERE ms.method_id = m.id AND ms.deleted_at IS NULL
                       AND ms.type = 'substance')          AS n_individual,
                   (SELECT COUNT(*) FROM method_substances ms
                     WHERE ms.method_id = m.id AND ms.deleted_at IS NULL
                       AND ms.type = 'substance_group')    AS n_grouped,
                   (SELECT COUNT(*) FROM substance_groups sg
                     WHERE sg.method_id = m.id)            AS n_groups,
                   (SELECT COUNT(*) FROM method_substances ms
                     WHERE ms.method_id = m.id AND ms.deleted_at IS NULL
                       AND ms.type = 'internal_standard')  AS n_internal_standards
            FROM methods m
            WHERE m.deleted_at IS NULL
            ORDER BY m.name
            """
        )
        methods = cursor.fetchall()
        cursor.close()
    finally:
        conn.close()
    return render_template("method/index.html", methods=methods)


@bp.route("/new", methods=["POST"])
def new():
    name = (request.form.get("name") or "").strip()
    version = (request.form.get("version") or "1").strip()
    if not name:
        flash("The method name is required.", "error")
        return redirect(url_for("method.index"))

    conn = analysis_conn()
    try:
        cursor = conn.cursor()
        cursor.execute("INSERT INTO methods (name, version) VALUES (%s, %s)", (name, version))
        method_id = int(cursor.lastrowid)
        cursor.close()
        conn.commit()
    except Exception as exc:  # the name is unique
        flash(f"Not created: {exc}", "error")
        return redirect(url_for("method.index"))
    finally:
        conn.close()
    return redirect(url_for("method.detail", method_id=method_id))


@bp.route("/<int:method_id>")
def detail(method_id: int):
    conn = analysis_conn()
    lake = lake_conn()
    try:
        cursor = conn.cursor(dictionary=True)
        cursor.execute("SELECT * FROM methods WHERE id = %s", (method_id,))
        method = cursor.fetchone()
        if not method:
            flash("Method not found.", "error")
            return redirect(url_for("method.index"))

        cursor.execute(
            """
            SELECT ms.*, ar.code AS analysis_rule_code, sf.name AS scan_filter
            FROM method_substances ms
            LEFT JOIN analysis_rules ar ON ar.id = ms.analysis_rule_id
            LEFT JOIN scan_filters sf   ON sf.id = ms.scan_filter_id
            WHERE ms.method_id = %s AND ms.deleted_at IS NULL
            ORDER BY ms.type DESC, ms.name
            """,
            (method_id,),
        )
        substances = cursor.fetchall()

        cursor.execute("SELECT * FROM analysis_rules WHERE deleted_at IS NULL ORDER BY id")
        rules = cursor.fetchall()

        # The IS gate: the internal standard that, if it fails, aborts the whole
        # sample.
        cursor.execute(
            "SELECT internal_standard FROM cutoff_processing_groups WHERE method_id = %s LIMIT 1",
            (method_id,),
        )
        row = cursor.fetchone()
        gate_id = row["internal_standard"] if row else None

        cursor.execute(
            """
            SELECT sg.id, sg.name,
                   GROUP_CONCAT(ms.name, ' + ') AS members,
                   COUNT(sgm.id)                AS n
            FROM substance_groups sg
            LEFT JOIN substance_group_memberships sgm ON sgm.substance_group_id = sg.id
            LEFT JOIN method_substances ms            ON ms.id = sgm.method_substance_id
            WHERE sg.method_id = %s
            GROUP BY sg.id
            ORDER BY sg.id
            """,
            (method_id,),
        )
        groups = cursor.fetchall()
        cursor.close()

        # Filters observed in the samples already ingested - the source of truth
        # for the scan_filter field.
        lake_cursor = lake.cursor()
        lake_cursor.execute("SELECT scan_filter FROM scan_filters ORDER BY scan_filter")
        observed_filters = [row[0] for row in lake_cursor.fetchall()]
        lake_cursor.close()
    finally:
        conn.close()
        lake.close()

    return render_template(
        "method/detail.html",
        method=method,
        substances=substances,
        rules=rules,
        gate_id=gate_id,
        groups=groups,
        observed_filters=observed_filters,
        types=TYPES,
    )


@bp.route("/<int:method_id>/substance", methods=["POST"])
def save_substance(method_id: int):
    substance_id = request.form.get("id", type=int)
    conn = analysis_conn()
    try:
        fields = {
            "name": (request.form.get("name") or "").strip(),
            "type": request.form.get("type", "substance"),
            "analysis_rule_id": request.form.get("analysis_rule_id", type=int),
            "competition": 1 if request.form.get("competition") else 0,
            "beta_blocker": 1 if request.form.get("beta_blocker") else 0,
            "start_time": _float("start_time", required=True),
            "end_time": _float("end_time", required=True),
            "mz": _float("mz", required=True),
            "mass_error_ppm": _float("mass_error_ppm") or 6.0,
            "scan_filter_id": _scan_filter_id(conn, request.form.get("scan_filter")),
            "fortified_concentration": _float("fortified_concentration"),
            "cut_off": _float("cut_off"),
            "r2_threshold": _float("r2_threshold"),
            "dtw_limit": _float("dtw_limit"),
            "any_points_ratio_limit": _float("any_points_ratio_limit"),
            "smoothing_type": _int("smoothing_type", 2),
            "smoothing_value": _int("smoothing_value", 7),
        }
        if not fields["name"]:
            raise ValueError("the substance name is required")
        if fields["type"] not in TYPES:
            raise ValueError("invalid type")

        cursor = conn.cursor()
        if substance_id:
            assignments = ", ".join(f"{key} = %s" for key in fields)
            cursor.execute(
                f"UPDATE method_substances SET {assignments}, updated_at = NOW() WHERE id = %s",
                (*fields.values(), substance_id),
            )
            flash("Substance updated.", "ok")
        else:
            columns = ", ".join(fields)
            placeholders = ", ".join("%s" for _ in fields)
            cursor.execute(
                f"INSERT INTO method_substances (method_id, {columns}) VALUES (%s, {placeholders})",
                (method_id, *fields.values()),
            )
            flash("Substance created.", "ok")
        cursor.close()
        conn.commit()
    except Exception as exc:  # noqa: BLE001 - becomes a message on screen
        flash(f"Not saved: {exc}", "error")
    finally:
        conn.close()
    return redirect(url_for("method.detail", method_id=method_id))


@bp.route("/<int:method_id>/substance/<int:substance_id>/delete", methods=["POST"])
def delete_substance(method_id: int, substance_id: int):
    """Soft delete - the pipeline filters on `deleted_at IS NULL`.

    Deliberately not a real delete: old results reference the substance by name,
    and the task snapshot holds the parameters that applied. Dropping the row
    would not erase the past, only its context.
    """
    conn = analysis_conn()
    try:
        cursor = conn.cursor()
        cursor.execute("UPDATE method_substances SET deleted_at = NOW() WHERE id = %s", (substance_id,))
        cursor.close()
        conn.commit()
    finally:
        conn.close()
    flash("Substance removed from the method.", "ok")
    return redirect(url_for("method.detail", method_id=method_id))


@bp.route("/<int:method_id>/gate", methods=["POST"])
def set_gate(method_id: int):
    """Define which internal standard is the method's gate.

    One per method: if it fails on a sample, the pipeline aborts that whole
    sample (`blocked_by_is_gate`) - there is no point looking for substances in
    a run that did not work. The other internal standards only raise alerts.
    """
    internal_standard_id = request.form.get("internal_standard", type=int)
    conn = analysis_conn()
    try:
        cursor = conn.cursor()
        cursor.execute("DELETE FROM cutoff_processing_groups WHERE method_id = %s", (method_id,))
        if internal_standard_id:
            cursor.execute(
                "INSERT INTO cutoff_processing_groups (method_id, internal_standard) VALUES (%s, %s)",
                (method_id, internal_standard_id),
            )
        cursor.close()
        conn.commit()
    finally:
        conn.close()
    flash("IS gate updated.", "ok")
    return redirect(url_for("method.detail", method_id=method_id))


@bp.route("/<int:method_id>/group", methods=["POST"])
def save_group(method_id: int):
    """Create a substance group (fragments of the same molecule).

    The group's consolidated result depends on the sample: in competition, it is
    SUSPECT if ANY member is; out of competition, only if ALL of them are.
    """
    members = request.form.getlist("members", type=int)
    name = (request.form.get("name") or "").strip() or None
    if len(members) < 2:
        flash("A group needs at least two substances.", "error")
        return redirect(url_for("method.detail", method_id=method_id))

    conn = analysis_conn()
    try:
        cursor = conn.cursor()
        cursor.execute("INSERT INTO substance_groups (method_id, name) VALUES (%s, %s)", (method_id, name))
        group_id = int(cursor.lastrowid)
        for substance_id in members:
            cursor.execute(
                "INSERT INTO substance_group_memberships (substance_group_id, method_substance_id) VALUES (%s, %s)",
                (group_id, substance_id),
            )
        cursor.close()
        conn.commit()
    finally:
        conn.close()
    flash("Group created.", "ok")
    return redirect(url_for("method.detail", method_id=method_id))


@bp.route("/<int:method_id>/group/<int:group_id>/delete", methods=["POST"])
def delete_group(method_id: int, group_id: int):
    conn = analysis_conn()
    try:
        cursor = conn.cursor()
        cursor.execute("DELETE FROM substance_group_memberships WHERE substance_group_id = %s", (group_id,))
        cursor.execute("DELETE FROM substance_groups WHERE id = %s", (group_id,))
        cursor.close()
        conn.commit()
    finally:
        conn.close()
    flash("Group removed.", "ok")
    return redirect(url_for("method.detail", method_id=method_id))


@bp.route("/api/scan-filters")
def api_scan_filters():
    conn = lake_conn()
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT scan_filter FROM scan_filters ORDER BY scan_filter")
        return jsonify([row[0] for row in cursor.fetchall()])
    finally:
        conn.close()
