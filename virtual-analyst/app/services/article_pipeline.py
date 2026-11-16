"""
article_pipeline.py - the interpretation steps exactly as the article describes
them (Sections 2.3 and 2.5-2.8, Table 1, Eq. 4).

Everything here is pure: arrays in, verdict out. No database, no files. The
caller (virtual_analyst.py) extracts the chromatograms, subtracts the NQC and
hands the traces over; this module decides.

    Pre-processing (Section 2.3)
        * every trace is resampled onto a uniform grid of GRID_POINTS = 3000
          points over the chromatographic axis [GRID_RT_MIN, GRID_RT_MAX];
        * 7-point smoothing with the MSFileReader-equivalent kernel, set per
          substance - applied upstream by the data provider, before the grid;
        * mass tolerance of +/-6 ppm - per substance, applied upstream;
        * NQC subtraction - applied upstream (adjust_intensity);
        * everything below is evaluated inside the substance's retention-time
          window (0.8 min wide in the screening method).

    Step 1 - centroid protocols
        Protocol #1: sample apex above the PQC apex          -> High
        Protocol #2: substances with no MRL flagged "any points": PQC/sample
                     apex ratio <= ratio limit (100)          -> Low
        Two-peak rule: co-eluting isomers - each PQC peak is evaluated on its
                     own and the most severe verdict is kept.

    Step 2 - Mamdani inference (Scikit-Fuzzy) over the Pseudo R2 ONLY
        Pseudo R2 = r2_score(y_true=sample, y_pred=PQC) on max-normalised
        traces. Five membership functions (Very Low ... Very High, Table 1),
        one rule per level, min implication, max aggregation, CENTROID
        defuzzification.

    Step 3 - DTW as a RESCUE when the Pseudo R2 fails
        DTW distance normalised by (n + m) (Eq. 4). Bands:
            DTW <  0.002 -> High;  < 0.008 -> Medium;  <= 0.120 -> Low;
            above 0.120  -> Negative.

    Step 4 - substances with an MRL (cut-off)
        C_hat estimated against the internal standard closest in tR;
        ratio = C_hat / cut-off, bands 1-1.35 (Low) / 1.30-1.60 (Medium) /
        >= 1.55 (High); below 1 -> Negative.

    Output levels: NEGATIVE, SUSPECT_VERY_LOW, SUSPECT_LOW, SUSPECT_MED,
                   SUSPECT_HIGH, SUSPECT_VERY_HIGH (INVALID when the PQC itself
                   carries no signal, i.e. there is no reference to compare to).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import skfuzzy as fuzz
from sklearn.metrics import r2_score

# ---------------------------------------------------------------------------
# Section 2.3 - the grid
# ---------------------------------------------------------------------------
GRID_POINTS = 3000
GRID_RT_MIN = 0.0   # min
GRID_RT_MAX = 12.0  # min - the screening run (substances elute 0.5-11.1 min)

# ---------------------------------------------------------------------------
# Step 1 - centroid protocols
# ---------------------------------------------------------------------------
PROTOCOL2_RATIO_LIMIT = 100.0   # PQC / sample apex ratio
TWO_PEAK_MIN_SEPARATION = 0.08  # min - two PQC apexes closer than this are one peak

# ---------------------------------------------------------------------------
# Step 2 - Table 1: membership functions of the Mamdani system
# ---------------------------------------------------------------------------
R2_DOMAIN = (-12.0, 1.0)
R2_UNIVERSE = np.round(np.arange(R2_DOMAIN[0], R2_DOMAIN[1] + 1e-9, 0.001), 6)
OUT_UNIVERSE = np.round(np.arange(0.0, 6.5 + 1e-9, 0.001), 6)

LEVELS = ("VERY_LOW", "LOW", "MED", "HIGH", "VERY_HIGH")

# input: Pseudo R2
INPUT_MF = {
    "VERY_LOW":  ("trap", [-12.0, -12.0, -8.0, -4.0]),
    "LOW":       ("tri",  [-8.0, -4.0, -2.0]),
    "MED":       ("tri",  [-4.0, -2.0, 0.0]),
    "HIGH":      ("tri",  [-2.0, 0.0, 1.0]),
    "VERY_HIGH": ("tri",  [0.0, 1.0, 2.0]),
}
# output: suspicion level, 1 (Very Low) ... 5 (Very High)
OUTPUT_MF = {
    "VERY_LOW":  [1.0, 1.0, 2.0],
    "LOW":       [2.0, 2.0, 3.0],
    "MED":       [3.0, 3.0, 4.0],
    "HIGH":      [4.0, 4.0, 5.0],
    "VERY_HIGH": [5.0, 5.0, 6.0],
}
# rule base: IF Pseudo R2 is <level> THEN suspicion is <level>
RULES = [(lvl, lvl) for lvl in LEVELS]

# ---------------------------------------------------------------------------
# Step 3 - DTW bands (Eq. 4 distance)
# ---------------------------------------------------------------------------
DTW_HIGH = 0.002
DTW_MED = 0.008
DTW_LOW = 0.120

# ---------------------------------------------------------------------------
# Step 4 - C_hat / cut-off membership functions
# ---------------------------------------------------------------------------
CUTOFF_MF = {
    "NEGATIVE": ("trap", [0.00, 0.00, 0.95, 1.00]),
    "LOW":      ("tri",  [1.00, 1.15, 1.35]),
    "MED":      ("tri",  [1.30, 1.45, 1.60]),
    "HIGH":     ("trap", [1.55, 1.70, 10.10, 10.10]),
}

SEVERITY = {
    "SUSPECT_VERY_HIGH": 5,
    "SUSPECT_HIGH": 4,
    "SUSPECT_MED": 3,
    "SUSPECT_LOW": 2,
    "SUSPECT_VERY_LOW": 1,
    "NEGATIVE": 0,
    "INVALID": -1,
}


def most_severe(verdicts):
    return max(verdicts, key=lambda v: SEVERITY.get(v, -2))


# ===========================================================================
# Pre-processing
# ===========================================================================
def make_grid(points: int = GRID_POINTS, rt_min: float = GRID_RT_MIN,
              rt_max: float = GRID_RT_MAX) -> np.ndarray:
    return np.linspace(rt_min, rt_max, points)


def window_mask(grid: np.ndarray, start_time: float, end_time: float) -> np.ndarray:
    return (grid >= float(start_time)) & (grid <= float(end_time))


def to_grid(rt, intensity, grid_rt: np.ndarray) -> np.ndarray:
    """Resample one (rt, intensity) trace onto `grid_rt`. Outside the span of
    the measured trace the signal is 0 (nothing was measured there)."""
    rt = np.asarray(rt, dtype=float)
    y = np.asarray(intensity, dtype=float)
    if rt.size == 0:
        return np.zeros_like(grid_rt, dtype=float)
    order = np.argsort(rt)
    rt, y = rt[order], np.nan_to_num(y[order])
    return np.interp(grid_rt, rt, y, left=0.0, right=0.0)


def normalise_to_max(y: np.ndarray) -> Optional[np.ndarray]:
    m = float(np.max(y)) if y.size else 0.0
    if not np.isfinite(m) or m <= 0:
        return None
    return y / m


# ===========================================================================
# Similarity measures
# ===========================================================================
def pseudo_r2(sample: np.ndarray, pqc: np.ndarray) -> Optional[float]:
    """r2_score(y_true=sample, y_pred=PQC) on traces normalised to their max."""
    s = normalise_to_max(sample)
    p = normalise_to_max(pqc)
    if s is None or p is None or s.size < 2:
        return None
    if np.allclose(s, s[0]):          # no variance in y_true -> undefined
        return None
    value = float(r2_score(s, p))
    return value if np.isfinite(value) else None


def dtw_normalised(a: np.ndarray, b: np.ndarray) -> Optional[float]:
    """Eq. 4: DTW(A, B) / (n + m), absolute-difference local cost, on the same
    max-normalised traces used by the Pseudo R2."""
    A = normalise_to_max(np.asarray(a, dtype=float))
    B = normalise_to_max(np.asarray(b, dtype=float))
    if A is None or B is None:
        return None
    n, m = len(A), len(B)
    prev = np.full(m + 1, np.inf)
    prev[0] = 0.0
    for i in range(1, n + 1):
        cost = np.abs(A[i - 1] - B)
        cur = np.empty(m + 1)
        cur[0] = np.inf
        diag_up = np.minimum(prev[1:], prev[:-1])   # min(D[i-1][j], D[i-1][j-1])
        for j in range(1, m + 1):                     # D[i][j-1] is sequential
            left = cur[j - 1]
            best = diag_up[j - 1] if diag_up[j - 1] < left else left
            cur[j] = cost[j - 1] + best
        prev = cur
    return float(prev[m] / (n + m))


# ===========================================================================
# Step 2 - Mamdani over the Pseudo R2
# ===========================================================================
def _mf_array(universe, kind, params):
    return fuzz.trapmf(universe, params) if kind == "trap" else fuzz.trimf(universe, params)


_IN_ARRAYS = {lvl: _mf_array(R2_UNIVERSE, *INPUT_MF[lvl]) for lvl in LEVELS}
_OUT_ARRAYS = {lvl: fuzz.trimf(OUT_UNIVERSE, OUTPUT_MF[lvl]) for lvl in LEVELS}


@dataclass
class FuzzyResult:
    level: Optional[str]          # e.g. "HIGH"; None when no rule fires
    centroid: Optional[float]
    firing: dict = field(default_factory=dict)


def mamdani_r2(r2: Optional[float]) -> FuzzyResult:
    """Fuzzify the Pseudo R2, fire the 5 rules (min implication), aggregate
    (max) and defuzzify by centroid. The crisp centroid is mapped back to the
    nearest output level (1 = Very Low ... 5 = Very High)."""
    if r2 is None or not np.isfinite(r2) or r2 < R2_DOMAIN[0]:
        return FuzzyResult(None, None, {})
    x = min(float(r2), R2_DOMAIN[1])
    firing = {lvl: float(fuzz.interp_membership(R2_UNIVERSE, _IN_ARRAYS[lvl], x)) for lvl in LEVELS}
    aggregated = np.zeros_like(OUT_UNIVERSE)
    for antecedent, consequent in RULES:
        aggregated = np.fmax(aggregated, np.fmin(firing[antecedent], _OUT_ARRAYS[consequent]))
    if aggregated.sum() <= 0:
        return FuzzyResult(None, None, firing)
    centroid = float(fuzz.defuzz(OUT_UNIVERSE, aggregated, "centroid"))
    idx = int(np.floor(centroid + 0.5)) - 1
    idx = max(0, min(len(LEVELS) - 1, idx))
    return FuzzyResult(LEVELS[idx], centroid, firing)


# ===========================================================================
# Step 3 - DTW rescue bands
# ===========================================================================
def dtw_level(d: Optional[float]) -> Optional[str]:
    if d is None or not np.isfinite(d):
        return None
    if d < DTW_HIGH:
        return "HIGH"
    if d < DTW_MED:
        return "MED"
    if d <= DTW_LOW:
        return "LOW"
    return None


# ===========================================================================
# Step 4 - C_hat / cut-off
# ===========================================================================
def _mf_scalar(kind, params, x):
    a = np.array([x], dtype=float)
    return float(_mf_array(a, kind, params)[0])


def cutoff_level(ratio: Optional[float]) -> str:
    """Bands 1-1.35 / 1.30-1.60 / >= 1.55 (overlapping memberships; the level
    with the highest membership wins, ties go to the more severe level)."""
    if ratio is None or not np.isfinite(ratio):
        return "NEGATIVE"
    if ratio <= 0:
        return "NEGATIVE"
    if ratio >= 10.0:
        return "SUSPECT_HIGH"
    mu = {k: _mf_scalar(kind, p, ratio) for k, (kind, p) in CUTOFF_MF.items()}
    if mu["HIGH"] >= mu["MED"] and mu["HIGH"] >= mu["LOW"] and mu["HIGH"] >= mu["NEGATIVE"]:
        return "SUSPECT_HIGH"
    if mu["MED"] >= mu["LOW"] and mu["MED"] >= mu["NEGATIVE"]:
        return "SUSPECT_MED"
    if mu["LOW"] >= mu["NEGATIVE"]:
        return "SUSPECT_LOW"
    return "NEGATIVE"


# ===========================================================================
# One peak through Steps 1-4
# ===========================================================================
@dataclass
class PeakDecision:
    verdict: str
    step: str                     # which step decided
    sample_max: float = 0.0
    pqc_max: float = 0.0
    r2: Optional[float] = None
    fuzzy: Optional[FuzzyResult] = None
    dtw: Optional[float] = None
    ratio_cutoff: Optional[float] = None
    concentration: Optional[float] = None

    def describe(self) -> str:
        parts = [f"step={self.step}", f"smp={self.sample_max:.3e}", f"pqc={self.pqc_max:.3e}"]
        if self.r2 is not None:
            parts.append(f"pseudoR2={self.r2:.4f}")
        if self.fuzzy is not None and self.fuzzy.centroid is not None:
            parts.append(f"centroid={self.fuzzy.centroid:.3f}")
        if self.dtw is not None:
            parts.append(f"DTW={self.dtw:.5f}")
        if self.ratio_cutoff is not None:
            parts.append(f"C/cutoff={self.ratio_cutoff:.3f}")
        return " ".join(parts) + f" -> {self.verdict}"


def evaluate_peak(
    sample: np.ndarray,
    pqc: np.ndarray,
    *,
    r2_threshold: Optional[float] = None,
    any_points: bool = False,
    ratio_limit: float = PROTOCOL2_RATIO_LIMIT,
    cutoff_estimator: Optional[Callable[[], tuple]] = None,
) -> PeakDecision:
    """Steps 1-4 for one peak. `sample` and `pqc` are the NQC-corrected traces
    on the grid, restricted to the evaluation window.

    `cutoff_estimator` is given only for substances with an MRL; it returns
    (C_hat, cut_off) - or (None, cut_off) when C_hat cannot be estimated.
    """
    has_mrl = cutoff_estimator is not None
    smp_max = float(np.max(sample)) if sample.size else 0.0
    pqc_max = float(np.max(pqc)) if pqc.size else 0.0
    d = PeakDecision("NEGATIVE", "step1", smp_max, pqc_max)

    if not np.isfinite(pqc_max) or pqc_max <= 0:
        d.verdict, d.step = "INVALID", "no_pqc_signal"
        return d
    if not np.isfinite(smp_max) or smp_max <= 0:
        d.verdict, d.step = "NEGATIVE", "no_sample_signal"
        return d

    # ---- Step 1 ----
    if smp_max > pqc_max:
        d.verdict, d.step = "SUSPECT_HIGH", "protocol1"
        return d
    if any_points and not has_mrl and (pqc_max / smp_max) <= float(ratio_limit):
        d.verdict, d.step = "SUSPECT_LOW", "protocol2"
        return d

    # ---- Step 2 ----
    d.r2 = pseudo_r2(sample, pqc)
    passes_threshold = d.r2 is not None and (r2_threshold is None or d.r2 >= float(r2_threshold))
    level = None
    if passes_threshold:
        d.fuzzy = mamdani_r2(d.r2)
        level = d.fuzzy.level
        if level is not None:
            d.step = "fuzzy_r2"

    # ---- Step 3 (rescue) ----
    if level is None:
        d.dtw = dtw_normalised(sample, pqc)
        level = dtw_level(d.dtw)
        d.step = "dtw" if level is not None else "all_failed"

    if level is None:
        d.verdict = "NEGATIVE"
        return d

    # ---- Step 4 ----
    if has_mrl:
        c_hat, cut_off = cutoff_estimator()
        d.concentration = c_hat
        if c_hat is None or cut_off is None or float(cut_off) <= 0:
            d.verdict, d.step = "NEGATIVE", d.step + "+cutoff_no_estimate"
            return d
        d.ratio_cutoff = float(c_hat) / float(cut_off)
        d.verdict = cutoff_level(d.ratio_cutoff)
        d.step += "+cutoff"
        return d

    d.verdict = f"SUSPECT_{level}"
    return d


# ===========================================================================
# Two-peak rule
# ===========================================================================
def split_two_peaks(grid_rt: np.ndarray, pqc: np.ndarray,
                    min_sep: float = TWO_PEAK_MIN_SEPARATION):
    """Return [(lo, hi), (lo, hi)] index ranges (relative to the window) for
    the two PQC peaks, split at the valley between them; None if the PQC does
    not show two separate peaks."""
    y = np.asarray(pqc, dtype=float)
    if y.size < 5:
        return None
    idx = np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]))[0] + 1
    idx = idx[y[idx] > 0]
    if idx.size < 2:
        return None
    apexes = []
    for i in idx[np.argsort(y[idx])[::-1]]:
        if all(abs(grid_rt[i] - grid_rt[j]) >= min_sep for j in apexes):
            apexes.append(int(i))
        if len(apexes) == 2:
            break
    if len(apexes) < 2:
        return None
    a, b = sorted(apexes)
    valley = a + int(np.argmin(y[a:b + 1]))
    return [(0, valley + 1), (valley, y.size)]
