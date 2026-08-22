r"""Build a local method from the JSON produced by tools/export_method.py.

    venv\Scripts\python tools\import_method.py data\method_tvii.json --name Method

It creates the method, its substances, the scan filters they reference, the
substance groups and the IS gate. Source ids are never reused - everything is
remapped, so the same JSON can be imported into any database.

Three mappings carry the weight, and each one silently breaks the analysis if it
is wrong:

  * `type` is lowercased. The source stores 'SUBSTANCE'; the pipeline queries
    `type = 'substance'`. MariaDB does not care about case, SQLite does - the
    column is COLLATE NOCASE for that reason, and this is the second belt.
  * `analysis_rule` is matched by CODE, never by id. Ids line up today by
    coincidence; the code is what the pipeline actually dispatches on.
  * `scan_filter` is copied as TEXT, character for character. It has to equal
    the filter string observed in the mzML, or the substance finds no channel
    and comes out INVALID - which produces no row at all.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from app.db.sqlite import analysis_conn, lake_conn

# Columns copied straight from the source row into method_substances.
COPIED = (
    "name", "version", "competition", "beta_blocker",
    "start_time", "end_time", "mz", "mass_error_ppm",
    "fortified_concentration", "cut_off",
    "r2_threshold", "dtw_limit", "any_points_ratio_limit",
    "smoothing_type", "smoothing_value",
)


def get_or_create_scan_filter(cursor, name: str | None, cache: dict) -> int | None:
    if not name:
        return None
    if name in cache:
        return cache[name]
    cursor.execute("INSERT OR IGNORE INTO scan_filters (name) VALUES (%s)", (name,))
    cursor.execute("SELECT id FROM scan_filters WHERE name = %s", (name,))
    filter_id = int(cursor.fetchone()[0])
    cache[name] = filter_id
    return filter_id


def drop_existing(cursor, method_id: int) -> None:
    cursor.execute(
        """
        DELETE FROM substance_group_memberships
        WHERE substance_group_id IN (SELECT id FROM substance_groups WHERE method_id = %s)
        """,
        (method_id,),
    )
    cursor.execute("DELETE FROM substance_groups WHERE method_id = %s", (method_id,))
    cursor.execute("DELETE FROM cutoff_processing_groups WHERE method_id = %s", (method_id,))
    cursor.execute("DELETE FROM method_substances WHERE method_id = %s", (method_id,))
    cursor.execute("DELETE FROM methods WHERE id = %s", (method_id,))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_file")
    parser.add_argument("--name", default="Method", help="name of the method to create")
    parser.add_argument("--replace", action="store_true", help="replace it if the name already exists")
    parser.add_argument(
        "--drop-orphans", action="store_true",
        help="exclude members of groups the export cut in half, instead of "
             "promoting them to individually analysed substances",
    )
    args = parser.parse_args()

    payload = json.loads(Path(args.json_file).read_text(encoding="utf-8"))
    substances = payload["substances"]
    source = payload.get("source", {})

    # A substance typed 'SUBSTANCE_GROUP' is reached ONLY through its group. If
    # the export window cut its group in half, leaving the type alone would mean
    # the substance is imported, counted, and never analysed - a silent hole in
    # the method. Two honest ways out, and no third:
    #   default        promote it to an individually analysed substance. The
    #                  count stays right and the compound is judged on its own,
    #                  instead of jointly with the fragment it lost.
    #   --drop-orphans leave it out entirely, so the count reflects what is
    #                  actually analysed.
    orphan_ids = {
        member
        for group in payload.get("partial_groups", [])
        for member in group["member_ids"]
    }
    if orphan_ids and args.drop_orphans:
        substances = [s for s in substances if s["id"] not in orphan_ids]

    conn = analysis_conn()
    try:
        cursor = conn.cursor(dictionary=True)

        cursor.execute("SELECT id FROM methods WHERE name = %s", (args.name,))
        existing = cursor.fetchone()
        if existing:
            if not args.replace:
                print(f"a method named '{args.name}' already exists (id {existing['id']}). "
                      f"Use --replace to rebuild it.")
                return 1
            drop_existing(conn.cursor(), int(existing["id"]))
            print(f"replaced the existing '{args.name}' (id {existing['id']})")

        # The analysis rules are matched by code - see the module docstring.
        cursor.execute("SELECT id, code FROM analysis_rules WHERE deleted_at IS NULL")
        rule_by_code = {r["code"]: int(r["id"]) for r in cursor.fetchall()}

        writer = conn.cursor()
        writer.execute(
            "INSERT INTO methods (name, version) VALUES (%s, %s)",
            (args.name, str(source.get("method_version") or "1")),
        )
        method_id = int(writer.lastrowid)

        filter_cache: dict[str, int] = {}
        id_map: dict[int, int] = {}
        unknown_rules: Counter = Counter()

        columns = ", ".join(COPIED)
        placeholders = ", ".join("%s" for _ in COPIED)

        for row in substances:
            code = row.get("analysis_rule_code")
            rule_id = rule_by_code.get(code) if code else None
            if code and rule_id is None:
                unknown_rules[code] += 1

            # See the orphan note above: promoted so that it is actually analysed.
            substance_type = str(row["type"]).lower()
            if row["id"] in orphan_ids and not args.drop_orphans:
                substance_type = "substance"

            writer.execute(
                f"""
                INSERT INTO method_substances
                    (method_id, type, analysis_rule_id, scan_filter_id, {columns})
                VALUES (%s, %s, %s, %s, {placeholders})
                """,
                (
                    method_id,
                    substance_type,
                    rule_id,
                    get_or_create_scan_filter(writer, row.get("scan_filter"), filter_cache),
                    *(row.get(column) for column in COPIED),
                ),
            )
            id_map[int(row["id"])] = int(writer.lastrowid)

        # Groups: fragments of the same molecule. A member missing from id_map
        # would mean the export and the group data disagree - better to stop than
        # to build a half group that quietly changes the aggregation.
        for group in payload.get("groups", []):
            members = [id_map[m] for m in group["member_ids"] if m in id_map]
            if len(members) != len(group["member_ids"]):
                raise RuntimeError(
                    f"group {group['source_group_id']} references a substance that "
                    f"was not exported - refusing to import a partial group"
                )
            writer.execute("INSERT INTO substance_groups (method_id) VALUES (%s)", (method_id,))
            group_id = int(writer.lastrowid)
            for member in members:
                writer.execute(
                    "INSERT INTO substance_group_memberships (substance_group_id, method_substance_id) "
                    "VALUES (%s, %s)",
                    (group_id, member),
                )

        gate_source = payload.get("is_gate_substance_id")
        gate_local = id_map.get(gate_source) if gate_source else None
        if gate_local:
            writer.execute(
                "INSERT INTO cutoff_processing_groups (method_id, internal_standard) VALUES (%s, %s)",
                (method_id, gate_local),
            )

        writer.close()
        cursor.close()
        conn.commit()
    finally:
        conn.close()

    # ---- report ----------------------------------------------------------
    by_type = Counter(
        "substance" if (s["id"] in orphan_ids and not args.drop_orphans) else str(s["type"]).lower()
        for s in substances
    )
    print(f"\nmethod '{args.name}' created (id {method_id}) from {source.get('method_name')}")
    print(f"  substances: {len(substances)}  {dict(by_type)}")
    print(f"  groups: {len(payload.get('groups', []))}")
    gate_name = next((s["name"] for s in substances if s["id"] == gate_source), None)
    print(f"  IS gate: {gate_name or 'NONE - nothing will block a failed run'}")
    analysed = by_type["substance"] + by_type["substance_group"]
    print(f"  analysed (substances + group members): {analysed}")
    print(f"  internal standards: {by_type['internal_standard']}")
    print(f"  scan filters: {len(filter_cache)}")
    for group in payload.get("partial_groups", []):
        action = "left out" if args.drop_orphans else "promoted to individual substances"
        names = ", ".join(
            next(s["name"] for s in payload["substances"] if s["id"] == m)
            for m in group["member_ids"]
        )
        print(f"\n  NOTE group {group['source_group_id']} was split by the export window.")
        print(f"    {names} -> {action}.")
        print(f"    In the source system it is judged jointly with the member(s) "
              f"left outside ({group['excluded_ids']}); here it is judged alone.")
    if unknown_rules:
        print(f"  WARNING unmapped analysis rules (they fall back to 'standard'): {dict(unknown_rules)}")

    # The classic cause of an empty result: a scan filter the ingested samples
    # never produced. Worth saying now rather than after a run returns nothing.
    lake = lake_conn()
    try:
        lake_cursor = lake.cursor()
        lake_cursor.execute("SELECT scan_filter FROM scan_filters")
        observed = {r[0] for r in lake_cursor.fetchall()}
        lake_cursor.close()
    finally:
        lake.close()

    if observed:
        unseen = sorted(set(filter_cache) - observed)
        if unseen:
            print(f"\n  {len(unseen)} scan filter(s) not observed in any ingested sample:")
            for name in unseen:
                print(f"    - {name}")
            print("  Substances on those filters will come out INVALID (and produce no row).")
        else:
            print("\n  every scan filter of the method was observed in the ingested samples")
    else:
        print("\n  no sample ingested yet - scan filters could not be cross-checked")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
