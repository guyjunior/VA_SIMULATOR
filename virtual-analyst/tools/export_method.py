r"""Export a method from the production analysis database into a JSON file.

    <source-venv>\python tools\export_method.py --env <path-to>\.env --out data\method.json

WHY IT IS A SEPARATE SCRIPT
This app never holds the credentials of the system the method comes from, and
does not depend on a MySQL driver. This script does: it is meant to run with an
interpreter that already has both - the source system's own virtual environment
and its own `.env`. The JSON it produces is a portable, auditable artefact;
`tools/import_method.py` turns it into a method here.

The export deliberately records only WHAT the method is (name, version, id
window), never WHERE it came from: the artefact is meant to be publishable.

READ ONLY. It issues SELECTs and nothing else.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal
from datetime import date, datetime
from pathlib import Path

import mysql.connector
from dotenv import dotenv_values

def jsonable(value):
    if isinstance(value, Decimal):
        return float(value)
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    return value


def clean(row: dict) -> dict:
    return {k: jsonable(v) for k, v in row.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env", required=True,
                        help="path to the .env holding ANALYSIS_DB_HOST/USER/PASSWORD/NAME")
    parser.add_argument("--method-id", type=int, default=1)
    # An explicit id window keeps the export reproducible: the same command
    # always yields the same substance set, and the window is recorded in the
    # file so the artefact says what it contains.
    parser.add_argument("--id-min", type=int, default=None)
    parser.add_argument("--id-max", type=int, default=None)
    parser.add_argument("--out", default="data/method_export.json")
    args = parser.parse_args()

    env = dotenv_values(args.env)
    missing = [k for k in ("ANALYSIS_DB_HOST", "ANALYSIS_DB_USER", "ANALYSIS_DB_NAME") if not env.get(k)]
    if missing:
        print(f"{args.env} does not define: {', '.join(missing)}")
        return 1

    conn = mysql.connector.connect(
        host=env["ANALYSIS_DB_HOST"],
        port=int(env.get("ANALYSIS_DB_PORT") or 3306),
        user=env["ANALYSIS_DB_USER"],
        password=env.get("ANALYSIS_DB_PASSWORD") or "",
        database=env["ANALYSIS_DB_NAME"],
        connection_timeout=20,
    )
    cursor = conn.cursor(dictionary=True)
    method_id = args.method_id

    cursor.execute("SELECT id, name, version FROM methods WHERE id = %s", (method_id,))
    method = cursor.fetchone()
    if not method:
        print(f"method {method_id} not found in {env['ANALYSIS_DB_NAME']}")
        return 1

    # Only what is alive: the pipeline itself filters on deleted_at IS NULL, so
    # exporting soft-deleted rows would mean importing substances the source
    # system no longer analyses.
    window = ""
    params: list = [method_id]
    if args.id_min is not None:
        window += " AND ms.id >= %s"
        params.append(args.id_min)
    if args.id_max is not None:
        window += " AND ms.id <= %s"
        params.append(args.id_max)

    cursor.execute(
        f"""
        SELECT ms.*, ar.code AS analysis_rule_code, sf.name AS scan_filter
        FROM method_substances ms
        LEFT JOIN analysis_rules ar ON ar.id = ms.analysis_rule_id
        LEFT JOIN scan_filters sf   ON sf.id = ms.scan_filter_id
        WHERE ms.method_id = %s AND ms.deleted_at IS NULL{window}
        ORDER BY ms.id
        """,
        tuple(params),
    )
    substances = [clean(r) for r in cursor.fetchall()]

    # Groups are exported as lists of substance ids. The importer remaps them to
    # its own ids - the source ids mean nothing here.
    cursor.execute(
        """
        SELECT sg.id AS group_id, sgm.method_substance_id
        FROM substance_groups sg
        JOIN substance_group_memberships sgm ON sgm.substance_group_id = sg.id
        JOIN method_substances ms ON ms.id = sgm.method_substance_id
        WHERE sg.method_id = %s AND ms.deleted_at IS NULL
        ORDER BY sg.id, sgm.method_substance_id
        """,
        (method_id,),
    )
    groups: dict[int, list[int]] = {}
    for row in cursor.fetchall():
        groups.setdefault(row["group_id"], []).append(row["method_substance_id"])

    # A group whose members are not all inside the window would be exported
    # half-formed, quietly changing how that compound is judged. Split them out
    # explicitly so the importer can decide, and the operator can see it.
    exported_ids = {s["id"] for s in substances}
    whole, partial = {}, {}
    for gid, members in groups.items():
        inside = [m for m in members if m in exported_ids]
        if not inside:
            continue
        (whole if len(inside) == len(members) else partial)[gid] = {
            "member_ids": inside,
            "excluded_ids": [m for m in members if m not in exported_ids],
        }
    groups = whole

    # The IS gate. Without it the imported method would never block a sample,
    # and a failed run would be reported as a clean negative.
    cursor.execute(
        """
        SELECT cpg.internal_standard
        FROM cutoff_processing_groups cpg
        JOIN method_substances ms ON ms.id = cpg.internal_standard
        WHERE cpg.method_id = %s AND ms.deleted_at IS NULL
        LIMIT 1
        """,
        (method_id,),
    )
    gate_row = cursor.fetchone()

    cursor.close()
    conn.close()

    payload = {
        # No host and no database name: this file is meant to be publishable,
        # and naming the internal server in it would be a leak that travels.
        "source": {
            "method_id": method_id,
            "method_name": method["name"],
            "method_version": method["version"],
            "id_min": args.id_min,
            "id_max": args.id_max,
        },
        "substances": substances,
        "groups": [
            {"source_group_id": gid, "member_ids": data["member_ids"]}
            for gid, data in sorted(groups.items())
        ],
        # Groups the id window cut in half. Their surviving members would never
        # be analysed if left typed as group members - the importer promotes
        # them and says so.
        "partial_groups": [
            {"source_group_id": gid, "member_ids": data["member_ids"],
             "excluded_ids": data["excluded_ids"]}
            for gid, data in sorted(partial.items())
        ],
        "is_gate_substance_id": gate_row["internal_standard"] if gate_row else None,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    from collections import Counter
    by_type = Counter(s["type"] for s in substances)
    print(f"exported {len(substances)} substances from {method['name']} -> {out}")
    print(f"  by type: {dict(by_type)}")
    print(f"  groups: {len(payload['groups'])}  |  IS gate: {payload['is_gate_substance_id']}")
    for partial_group in payload["partial_groups"]:
        print(f"  PARTIAL group {partial_group['source_group_id']}: kept "
              f"{partial_group['member_ids']}, outside the window "
              f"{partial_group['excluded_ids']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
