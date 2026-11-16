"""Integration: process_single_substance end to end - extraction, NQC
subtraction, the 3000-point grid, the article's Steps 1-4, the nearest internal
standard and the INSERT into sample_results - on a temporary SQLite database
with synthetic chromatograms. No ClickHouse, no Docker, no real data."""
import numpy as np
import pytest

from app.db import schema
from app.db.sqlite import connect
from app.services import virtual_analyst as va


class FakeRaw:
    """Answers GetChroData like MSFileReader / DbDataProvider: one Gaussian
    peak per m/z, scans every 0.004 min."""

    def __init__(self, peaks):
        self.peaks = peaks          # mz -> (rt_centre, amplitude, sigma)

    def GetChroData(self, startTime, endTime, massRange1, scanFilter, smoothingType, smoothingValue):
        lo, hi = (float(v) for v in massRange1.split("-"))
        mz = round((lo + hi) / 2, 4)
        rt = np.arange(float(startTime), float(endTime) + 1e-9, 0.004)
        y = np.zeros_like(rt)
        if mz in self.peaks:
            c, a, s = self.peaks[mz]
            y = a * np.exp(-0.5 * ((rt - c) / s) ** 2)
        return [[rt, y]]


class FakeCH:
    def insert(self, *args, **kwargs):
        pass


SUBSTANCES = {
    # name: (mz, rt window start, rule, cut_off, fortified)
    "match_shape":   (300.0, 3.0, "standard", None, None),
    "absent":        (310.0, 4.0, "standard", None, None),
    "above_pqc":     (320.0, 5.0, "standard", None, None),
    "mrl_substance": (330.0, 7.1, "standard", 5.0, 10.0),
    "any_points":    (340.0, 9.0, "any_points", None, None),
}
IS_NEAR = (400.0, 7.6)   # closest in tR to mrl_substance (window centre 7.5)
IS_FAR = (410.0, 1.6)


@pytest.fixture()
def db(tmp_path):
    conn = connect(tmp_path / "analysis.sqlite")
    cur = conn.cursor(dictionary=True)
    for ddl in schema.ANALYSIS_TABLES:
        cur.execute(ddl)
    for code, name, desc in schema.SEED_RULES:
        cur.execute("INSERT INTO analysis_rules (code, name, description) VALUES (%s, %s, %s)", (code, name, desc))
    cur.execute("INSERT INTO methods (name) VALUES ('test')")
    for name, (mz, start) in {"IS_near": IS_NEAR, "IS_far": IS_FAR}.items():
        cur.execute("""INSERT INTO method_substances (name, method_id, type, start_time, end_time, mz)
                       VALUES (%s, 1, 'internal_standard', %s, %s, %s)""", (name, start, start + 0.8, mz))
    for name, (mz, start, rule, cut, fort) in SUBSTANCES.items():
        cur.execute("""INSERT INTO method_substances
                       (name, method_id, analysis_rule_id, type, start_time, end_time, mz,
                        cut_off, fortified_concentration, r2_threshold, dtw_limit)
                       VALUES (%s, 1, (SELECT id FROM analysis_rules WHERE code = %s), 'substance',
                               %s, %s, %s, %s, %s, -50.0, 1.0)""",
                    (name, rule, start, start + 0.8, mz, cut, fort))
    cur.execute("""INSERT INTO sample_processing_tasks (samples, method_id, batch, method_snapshot, status)
                   VALUES ('[]', 1, 'synthetic', '{}', 'running')""")
    conn.commit()
    yield conn, cur
    conn.close()


def _pqc():
    peaks = {mz: (start + 0.4, 1e6, 0.03) for _, (mz, start, *_r) in SUBSTANCES.items()}
    peaks[IS_NEAR[0]] = (IS_NEAR[1] + 0.4, 5e5, 0.03)
    peaks[IS_FAR[0]] = (IS_FAR[1] + 0.4, 5e5, 0.03)
    return FakeRaw(peaks)


def _sample():
    s = SUBSTANCES
    return FakeRaw({
        s["match_shape"][0]: (s["match_shape"][1] + 0.4, 3e5, 0.03),
        s["above_pqc"][0]: (s["above_pqc"][1] + 0.4, 2e6, 0.03),
        # C_hat = (0.6e6 / 5e5) / (1e6 / 5e5) * 10 = 6.0 -> ratio 1.2 -> Low.
        # With the FAR internal standard (twice as intense) C_hat would be 3.0
        # -> Negative: the test fails if the wrong IS is picked.
        s["mrl_substance"][0]: (s["mrl_substance"][1] + 0.4, 6e5, 0.03),
        s["any_points"][0]: (s["any_points"][1] + 0.25, 2e4, 0.01),
        IS_NEAR[0]: (IS_NEAR[1] + 0.4, 5e5, 0.03),
        IS_FAR[0]: (IS_FAR[1] + 0.4, 1e6, 0.03),
    })


def _run(cur, name, tmp_path):
    cur.execute("""SELECT ms.*, ar.code AS analysis_rule_code, NULL AS scan_filter
                   FROM method_substances ms LEFT JOIN analysis_rules ar ON ar.id = ms.analysis_rule_id
                   WHERE ms.name = %s""", (name,))
    sub = va._attach_mass_range(cur.fetchone())
    blank = FakeRaw({})
    return va.process_single_substance(
        cur, 1, "db://synthetic/Data/S1.raw", _sample(), _pqc(), _pqc(), blank, sub,
        str(tmp_path / "test.log"), ch_client=FakeCH())


@pytest.mark.parametrize("name, verdict", [
    ("match_shape", "SUSPECT_VERY_HIGH"),   # Step 2: Pseudo R2 = 1 -> Mamdani -> Very High
    ("absent", "NEGATIVE"),
    ("above_pqc", "SUSPECT_HIGH"),          # Step 1: Protocol #1
    ("mrl_substance", "SUSPECT_LOW"),       # Step 4: C_hat/cut-off = 1.2
    ("any_points", "SUSPECT_LOW"),          # Step 1: Protocol #2 (PQC/sample = 50)
])
def test_substance_verdicts(db, tmp_path, name, verdict):
    conn, cur = db
    result_id, result, _, errors = _run(cur, name, tmp_path)
    assert errors == []
    assert result == verdict
    cur.execute("SELECT result, concentration_value FROM sample_results WHERE id = %s", (result_id,))
    row = cur.fetchone()
    assert row["result"] == verdict
    if name == "mrl_substance":
        assert row["concentration_value"] == pytest.approx(6.0, rel=1e-3)


def test_nearest_internal_standard_is_chosen_by_tr(db):
    _, cur = db
    cur.execute("SELECT * FROM method_substances WHERE name = 'mrl_substance'")
    assert va.nearest_internal_standard(cur, cur.fetchone())["name"] == "IS_near"


@pytest.mark.parametrize("members, expected", [
    (["SUSPECT_LOW", "SUSPECT_VERY_HIGH"], "SUSPECT"),   # every member suspect, any level
    (["SUSPECT_HIGH", "NEGATIVE"], "NEGATIVE"),
    (["SUSPECT_MED", "INVALID"], "NEGATIVE"),
    ([], "NEGATIVE"),
])
def test_group_result_needs_every_member_suspect(members, expected):
    assert va.group_result(members) == expected
