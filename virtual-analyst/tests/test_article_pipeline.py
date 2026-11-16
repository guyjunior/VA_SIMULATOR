"""Unit tests: each step of article_pipeline against the article's definitions."""
import numpy as np
import pytest
from sklearn.metrics import r2_score

from app.services import article_pipeline as ap


def gauss(x, mu, sigma, amp=1.0):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# --------------------------------------------------------------- pre-processing
def test_grid_has_3000_points():
    g = ap.make_grid()
    assert g.size == 3000
    assert g[0] == ap.GRID_RT_MIN and g[-1] == ap.GRID_RT_MAX


def test_window_of_0_8_min_on_grid():
    g = ap.make_grid()
    m = ap.window_mask(g, 5.0, 5.8)
    assert 195 <= m.sum() <= 202          # 0.8 min / (12 min / 2999)


def test_to_grid_is_zero_outside_measured_span():
    g = np.linspace(0, 10, 11)
    y = ap.to_grid([4.0, 5.0, 6.0], [1.0, 2.0, 1.0], g)
    assert y[5] == 2.0 and y[0] == 0.0 and y[10] == 0.0


# --------------------------------------------------------------- Pseudo R2
def test_pseudo_r2_is_r2_score_sample_true_pqc_pred_on_max_normalised():
    x = np.linspace(0, 1, 200)
    s = gauss(x, 0.5, 0.05, 3.0) + 0.1 * x
    p = gauss(x, 0.52, 0.06, 7.0)
    expected = r2_score(s / s.max(), p / p.max())
    assert ap.pseudo_r2(s, p) == pytest.approx(expected)
    # argument order matters: y_true is the SAMPLE
    assert ap.pseudo_r2(s, p) != pytest.approx(r2_score(p / p.max(), s / s.max()))


def test_pseudo_r2_is_scale_free():
    x = np.linspace(0, 1, 200)
    p = gauss(x, 0.5, 0.05, 1e6)
    assert ap.pseudo_r2(0.01 * p, p) == pytest.approx(1.0)


# --------------------------------------------------------------- Step 2 (Table 1)
def test_table1_has_five_membership_functions():
    assert list(ap.INPUT_MF) == ["VERY_LOW", "LOW", "MED", "HIGH", "VERY_HIGH"]
    assert list(ap.OUTPUT_MF) == ["VERY_LOW", "LOW", "MED", "HIGH", "VERY_HIGH"]


@pytest.mark.parametrize("r2, level", [
    (1.0, "VERY_HIGH"),
    (0.6, "VERY_HIGH"),
    (0.0, "HIGH"),
    (-1.0, "HIGH"),     # equal firing of MED and HIGH -> centroid 3.5 rounds to HIGH
    (-2.0, "MED"),
    (-4.0, "LOW"),
    (-6.0, "LOW"),
    (-9.0, "VERY_LOW"),
    (-12.0, "VERY_LOW"),
])
def test_mamdani_levels(r2, level):
    assert ap.mamdani_r2(r2).level == level


def test_mamdani_defuzzifies_by_centroid():
    # R2 = 0 fires only HIGH (degree 1): output = trimf [4,4,5], centroid 13/3
    res = ap.mamdani_r2(0.0)
    assert res.firing["HIGH"] == pytest.approx(1.0)
    assert sum(v for k, v in res.firing.items() if k != "HIGH") == pytest.approx(0.0)
    assert res.centroid == pytest.approx(13 / 3, abs=1e-3)


def test_mamdani_outside_domain_does_not_fire():
    assert ap.mamdani_r2(-12.5).level is None
    assert ap.mamdani_r2(None).level is None


# --------------------------------------------------------------- Step 3 (Eq. 4)
def test_dtw_is_normalised_by_n_plus_m():
    # cumulative cost 1 on the optimal path, n + m = 4
    assert ap.dtw_normalised(np.array([1.0, 0.0]), np.array([1.0, 1.0])) == pytest.approx(0.25)


def test_dtw_identical_is_zero():
    x = np.linspace(0, 1, 100)
    p = gauss(x, 0.5, 0.05)
    assert ap.dtw_normalised(p, 3 * p) == pytest.approx(0.0)


@pytest.mark.parametrize("d, level", [
    (0.0, "HIGH"), (0.0019, "HIGH"), (0.002, "MED"), (0.0079, "MED"),
    (0.008, "LOW"), (0.120, "LOW"), (0.1201, None),
])
def test_dtw_bands(d, level):
    assert ap.dtw_level(d) == level


# --------------------------------------------------------------- Step 4
@pytest.mark.parametrize("ratio, verdict", [
    (0.5, "NEGATIVE"), (0.97, "NEGATIVE"),
    (1.15, "SUSPECT_LOW"), (1.32, "SUSPECT_LOW"),
    (1.45, "SUSPECT_MED"), (1.52, "SUSPECT_MED"),
    (1.58, "SUSPECT_HIGH"),   # overlap 1.55-1.60: mu_high 0.20 > mu_med 0.13
    (1.70, "SUSPECT_HIGH"), (3.0, "SUSPECT_HIGH"), (25.0, "SUSPECT_HIGH"),
])
def test_cutoff_bands(ratio, verdict):
    assert ap.cutoff_level(ratio) == verdict


# --------------------------------------------------------------- Steps 1-4 together
X = np.linspace(0, 0.8, 200)
PQC = gauss(X, 0.4, 0.03, 1000.0)


def test_protocol1_sample_above_pqc_is_high():
    assert ap.evaluate_peak(1.2 * PQC, PQC).verdict == "SUSPECT_HIGH"


def test_protocol2_any_points_low_only_without_mrl():
    sample = gauss(X, 0.25, 0.01, 20.0)          # odd shape, 50x below the PQC
    d = ap.evaluate_peak(sample, PQC, any_points=True)
    assert (d.verdict, d.step) == ("SUSPECT_LOW", "protocol2")
    with_mrl = ap.evaluate_peak(sample, PQC, any_points=True,
                                cutoff_estimator=lambda: (0.1, 1.0))
    assert with_mrl.step != "protocol2"


def test_fuzzy_on_r2_when_shapes_match():
    d = ap.evaluate_peak(0.5 * PQC, PQC)
    assert d.step == "fuzzy_r2" and d.dtw is None
    assert d.verdict == "SUSPECT_VERY_HIGH"


def test_dtw_rescues_when_r2_fails():
    shifted = gauss(X, 0.43, 0.03, 500.0)        # same shape, shifted 0.03 min
    r2 = ap.pseudo_r2(shifted, PQC)
    d = ap.evaluate_peak(shifted, PQC, r2_threshold=r2 + 0.01)   # force the R2 to fail
    assert d.step == "dtw" and d.dtw is not None
    assert d.verdict == f"SUSPECT_{ap.dtw_level(d.dtw)}"


def test_negative_when_r2_and_dtw_fail():
    noise = np.abs(np.sin(np.linspace(0, 60, X.size))) * 5.0
    d = ap.evaluate_peak(noise, PQC, r2_threshold=0.9)
    assert d.verdict == "NEGATIVE" and d.step == "all_failed"


def test_mrl_level_comes_from_step4():
    d = ap.evaluate_peak(0.5 * PQC, PQC, cutoff_estimator=lambda: (6.0, 5.0))
    assert d.verdict == "SUSPECT_LOW" and d.ratio_cutoff == pytest.approx(1.2)
    d = ap.evaluate_peak(0.5 * PQC, PQC, cutoff_estimator=lambda: (4.0, 5.0))
    assert d.verdict == "NEGATIVE"


def test_no_pqc_signal_is_invalid_and_no_sample_signal_negative():
    assert ap.evaluate_peak(PQC, np.zeros_like(PQC)).verdict == "INVALID"
    assert ap.evaluate_peak(np.zeros_like(PQC), PQC).verdict == "NEGATIVE"


def test_two_peak_split_at_valley():
    pqc2 = gauss(X, 0.25, 0.02, 1000.0) + gauss(X, 0.55, 0.02, 800.0)
    (lo1, hi1), (lo2, hi2) = ap.split_two_peaks(X, pqc2)
    valley = hi1 - 1
    assert lo1 == 0 and hi2 == X.size and lo2 == valley
    assert 0.3 < X[valley] < 0.5
    assert ap.split_two_peaks(X, PQC) is None


def test_output_levels_are_the_articles():
    assert set(ap.SEVERITY) == {"NEGATIVE", "SUSPECT_VERY_LOW", "SUSPECT_LOW", "SUSPECT_MED",
                                "SUSPECT_HIGH", "SUSPECT_VERY_HIGH", "INVALID"}
