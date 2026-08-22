"""
build_method_snapshot.py - freezes the method configuration into JSON.

WHY IT EXISTS:
    When a `sample_processing_task` is created, we store alongside it a complete
    copy of ALL the method's substances and parameters (RT range, m/z range,
    smoothing, R2/DTW thresholds, cut-off, analysis rules, groups, and so on).
    That lives in `sample_processing_tasks.method_snapshot` as JSON.

    The reason: if someone edits the method later (changes a threshold, adds a
    substance), the old results stay interpretable and auditable - they were
    produced with the configuration frozen right there.

    The processing itself does NOT use the snapshot - it reads the method "live"
    from the database for every sample. The snapshot is for auditing and the UI.

RETURNED STRUCTURE:
    {
      "name":  "<method name>",
      "version": "<version>",
      "substances":       [ {... params ...}, ... ],   # substances NOT in a group
      "group_substances": [ {"group_id": N, "substances": [...] } ],
    }

    "substances" + "group_substances" cover the method's complete substance set
    without duplication - a substance in a group appears only there.

    The format is identical to the production system's - whatever comes out here
    can be diffed line by line against the snapshot from there.
"""

import os
from datetime import datetime
from decimal import Decimal
from typing import Any, Dict, List

from dotenv import load_dotenv

from app.db.sqlite import analysis_conn

load_dotenv()

def db_conn():
    """The ANALYSIS database - where the method is configured."""
    return analysis_conn()


def clean_for_json(data):
    """Recursively convert non-serialisable types into JSON-safe equivalents.

    - Decimal (coming from the database) -> float
    - datetime -> ISO string

    Without it, `json.dumps` blows up with a TypeError on fields like `cut_off`
    or the creation timestamps.
    """
    if isinstance(data, dict):
        return {k: clean_for_json(v) for k, v in data.items()}
    if isinstance(data, list):
        return [clean_for_json(v) for v in data]
    if isinstance(data, Decimal):
        return float(data)
    if isinstance(data, datetime):
        return data.isoformat()
    return data


def build_method_snapshot(method_id: int) -> Dict[str, Any]:
    """Read the whole method from the database and return a dict ready to serialise.

    It reads 3 groups of information:
      1. The method header (name, version).
      2. Every substance in the method (`method_substances`).
      3. Which substances are grouped in `substance_groups` (groups share a
         result - used for fragments of the same molecule).

    Substances that appear in any group are REMOVED from the "substances" list
    and go only into "group_substances", avoiding duplication.
    """
    if not method_id:
        raise ValueError("Valid method_id must be provided.")

    conn = db_conn()
    cursor = conn.cursor(dictionary=True)

    # ---- 1. Method header ----
    cursor.execute("SELECT id, name, version FROM methods WHERE id = %s", (method_id,))
    method = cursor.fetchone()
    if not method:
        cursor.close()
        conn.close()
        raise ValueError(f"Method with ID {method_id} does not exist.")

    # ---- 2. Raw list of ALL the method's substances ----
    # (the "is it in a group or not" filter happens in Python afterwards)
    # The snapshot stores mz + mass_error_ppm (not the mass_range string) -
    # anyone who needs the "min-max" window rebuilds it at runtime with
    # compose_mass_range.
    cursor.execute(
        """
        SELECT
            ms.id, ms.name, ms.type,
            ms.competition, ms.beta_blocker,
            ms.start_time, ms.end_time,
            ms.mz, ms.mass_error_ppm, sf.name AS scan_filter,
            ms.fortified_concentration, ms.cut_off,
            ms.r2_threshold, ms.dtw_limit, ms.any_points_ratio_limit,
            ms.smoothing_type, ms.smoothing_value,
            ms.analysis_rule_id, ar.code AS analysis_rule_code
        FROM method_substances ms
        LEFT JOIN analysis_rules ar ON ar.id = ms.analysis_rule_id
        LEFT JOIN scan_filters sf ON sf.id = ms.scan_filter_id
        WHERE ms.method_id = %s
          AND ms.deleted_at IS NULL
        """,
        (method_id,),
    )
    all_substances = cursor.fetchall()

    # ---- 3. Grouped substances (substance_groups + substance_group_memberships) ----
    # Each row is a substance tied to a group. The JOIN guarantees we get the
    # group_id plus every parameter of each member substance.
    cursor.execute(
        """
        SELECT
            sg.id AS group_id,
            sgm.method_substance_id,
            ms.name AS substance_name, ms.type,
            ms.competition, ms.beta_blocker,
            ms.start_time, ms.end_time,
            ms.mz, ms.mass_error_ppm, sf.name AS scan_filter,
            ms.fortified_concentration, ms.cut_off,
            ms.r2_threshold, ms.dtw_limit, ms.any_points_ratio_limit,
            ms.smoothing_type, ms.smoothing_value,
            ms.analysis_rule_id, ar.code AS analysis_rule_code
        FROM substance_groups sg
        JOIN substance_group_memberships sgm ON sgm.substance_group_id = sg.id
        JOIN method_substances ms ON ms.id = sgm.method_substance_id
        LEFT JOIN analysis_rules ar ON ar.id = ms.analysis_rule_id
        LEFT JOIN scan_filters sf ON sf.id = ms.scan_filter_id
        WHERE sg.method_id = %s
          AND ms.deleted_at IS NULL
        ORDER BY sg.id
        """,
        (method_id,),
    )
    group_rows = cursor.fetchall()

    # Ids of substances that are ALREADY in some group - used to exclude them
    # from the flat list below (so they do not appear in both snapshot fields).
    grouped_substance_ids: set = set()
    group_substances_data: List[Dict[str, Any]] = []

    if group_rows:
        # Regroup the rows into {group_id: [substances]}.
        groups_dict: Dict[int, Dict[str, Any]] = {}
        for row in group_rows:
            gid = row["group_id"]
            groups_dict.setdefault(gid, {"group_id": gid, "substances": []})
            grouped_substance_ids.add(row["method_substance_id"])
            groups_dict[gid]["substances"].append({
                "id": row["method_substance_id"],
                "name": row["substance_name"],
                "type": row["type"],
                "competition": bool(row["competition"]),
                "beta_blocker": bool(row["beta_blocker"]),
                "start_time": row["start_time"],
                "end_time": row["end_time"],
                "mz": row["mz"],
                "mass_error_ppm": row["mass_error_ppm"],
                "scan_filter": row["scan_filter"],
                "fortified_concentration": row["fortified_concentration"],
                "cut_off": row["cut_off"],
                "r2_threshold": row["r2_threshold"],
                "dtw_limit": row["dtw_limit"],
                "any_points_ratio_limit": row["any_points_ratio_limit"],
                "smoothing_type": row["smoothing_type"],
                "smoothing_value": row["smoothing_value"],
                "analysis_rule_id": row["analysis_rule_id"],
                "analysis_rule_code": row["analysis_rule_code"],
            })
        group_substances_data = list(groups_dict.values())

    # ---- Flat list - only substances that are NOT in a group ----
    substances_data = []
    for s in all_substances:
        if s["id"] in grouped_substance_ids:
            continue
        substances_data.append({
            "id": s["id"],
            "name": s["name"], "type": s["type"],
            "competition": bool(s["competition"]),
            "beta_blocker": bool(s["beta_blocker"]),
            "start_time": s["start_time"], "end_time": s["end_time"],
            "mz": s["mz"], "mass_error_ppm": s["mass_error_ppm"],
            "scan_filter": s["scan_filter"],
            "fortified_concentration": s["fortified_concentration"],
            "cut_off": s["cut_off"],
            "r2_threshold": s["r2_threshold"],
            "dtw_limit": s["dtw_limit"],
            "any_points_ratio_limit": s["any_points_ratio_limit"],
            "smoothing_type": s["smoothing_type"],
            "smoothing_value": s["smoothing_value"],
            "analysis_rule_id": s["analysis_rule_id"],
            "analysis_rule_code": s["analysis_rule_code"],
        })

    cursor.close()
    conn.close()

    # `clean_for_json` is the last step: it converts Decimal/datetime into
    # JSON-safe types before the dict reaches `json.dumps` in the caller.
    return clean_for_json({
        "name": method["name"],
        "version": method.get("version"),
        "substances": substances_data,
        "group_substances": group_substances_data,
    })
