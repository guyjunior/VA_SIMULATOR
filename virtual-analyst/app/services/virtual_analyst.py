"""
virtual_analyst.py - the chemistry core of the processing.

FAITHFUL PORT of the pipeline running in production. The chemistry was not
touched: the same functions, in the same order, with the same thresholds. Only
the plumbing changed - MariaDB became SQLite (through the compatibility layer in
app/db/sqlite.py, which accepts the same `%s` dialect) and ClickHouse is the one
from this project's docker-compose. Changing the chemistry here changes the
result; if you must, change it knowing that.

THIS IS THE PIPELINE THAT DECIDES WHETHER A SAMPLE IS SUSPECT OR NEGATIVE FOR
EACH SUBSTANCE IN THE METHOD. Identical to the V1 project's, with one critical
difference:

    V1   - takes .raw file paths and opens them via pymsfilereader.MSFileReader.
    HERE - takes sample ids from the lake and uses DbDataProvider, which passes
           for MSFileReader but fetches chromatograms from the database.

So the entire chemistry was reused WITH NO CHANGES - we only calibrated the
data_provider so that `GetChroData` returns XICs identical to MSFileReader's.

============================================================================
ONE SAMPLE, END TO END - process_sample_from_db(sample, controls, method_id)
============================================================================

    1. Create 4 DbDataProviders:
           raw_sample, raw_cqcl, raw_cqcl_reinj, raw_cqn
       (they look exactly like MSFileReader - any function taking `raw_*`
       cannot tell the difference)

    2. EVALUATE THE INTERNAL STANDARDS (gate)
       For every IS configured in method_substances (type='internal_standard'):
         - extract the XIC of the sample and of the CQCL
         - compute DTW + R2 over the CQCL peak window
         - if it passes: carry on | if it fails:
             * if it is the IS gate -> BLOCK the whole sample
             * if it is another IS  -> record it in internal_standards_audit
               (an alert in the UI)

    3. PER-SAMPLE FLAGS (they come from the lake's samples table)
         competition  -> False filters out substances/groups with competition=1
         beta_block   -> False filters out substances/groups with beta_blocker=1
       Both default-safe: NULL -> True. The 2 flags are orthogonal: a sample can
       pass or fail one without affecting the other. Intra-group aggregation
       (ANY vs ALL) is still governed by competition alone.

    4. PROCESS EACH INDIVIDUAL SUBSTANCE
       Loop over method_substances.type='substance' (filtered by competition and
       beta_blocker according to the sample's flags). For each one it calls
       process_single_substance, which:
         a. Extracts the XIC of sample, cqcl, cqcl_reinj, cqcn (4 chromatograms)
         b. Applies adjust_intensity (subtracts the CQN baseline)
         c. Locates the peak in the corrected CQCL (that gives the reference RT)
         d. Applies the analysis rule (analysis_rule_code):
              "standard"  -> R2 + DTW over the peak window (with a fallback to
                             the reinjection)
              "two_peaks" -> R2 + DTW over TWO peaks (bimodal substances)
              "any_points"-> sample/control height ratio over a fixed window
         e. If a cut_off is configured -> computes the concentration via the
            internal standard and classifies with fuzzy logic (SUSPECT_LOW/MED/
            HIGH or NEGATIVE)
         f. INSERTs into sample_results + stores the XICs in the analysis
            ClickHouse

    5. PROCESS THE SUBSTANCE GROUPS (substance_groups)
       Grouped substances are fragments of the same molecule. Each member
       becomes a row in sample_results, and then an update applies the
       consolidated `result_group`:
         - in_competition=True:  the group is SUSPECT if ANY member is SUSPECT.
         - in_competition=False: the group is SUSPECT only if ALL members are.

============================================================================
QUICK GLOSSARY
============================================================================
    XIC = Extracted Ion Chromatogram. A chromatogram filtered to an m/z window.
          X axis = retention time, Y axis = the summed intensity of the
          centroids inside that m/z window in each scan.

    CQCL = the "low" positive control (substance present at a low, known
           concentration). It is the reference - if the shape of the sample's
           peak resembles the CQCL's at the SAME RT, it is a SUSPECT candidate.
    CQCL_reinj = the CQCL re-injected (a second chance, method robustness).
    CQN  = the negative control (clean matrix). Used to subtract the baseline.
    BR   = lab blank - discarded.

    R2  = the coefficient of determination between the sample XIC and the CQCL
          XIC over the peak window. Close to 1 = perfectly correlated.
    DTW = Dynamic Time Warping. The distance between the two normalised curves;
          low = similar curves. Robust to small RT shifts.

    The standard (joint) rule: SUSPECT if R2 > r2_threshold AND DTW < dtw_limit.
    Every substance can carry its own thresholds - read from method_substances.
"""
import os
import sqlite3
import uuid
import traceback
from datetime import datetime
import time
import json
from typing import Optional
import numpy as np
import pandas as pd
from dotenv import load_dotenv
from sklearn.metrics import r2_score
from dtaidistance import dtw
from pathlib import Path

from app import config
from app.db.clickhouse import analysis_ch
from app.db.sqlite import analysis_conn

load_dotenv()

# Project root - the log/ folder is anchored here, never to the process cwd.
_REPO_ROOT = str(config.ROOT)


def maria_conn():
    """The ANALYSIS database (method, substances, tasks, results).

    The name stayed for fidelity with the original - in production this is
    MariaDB; here it is the SQLite file `data/analysis.sqlite`, with the same
    cursor interface.
    """
    return analysis_conn()


def click_conn():
    """The ANALYSIS ClickHouse database - where each run's XICs are stored."""
    return analysis_ch()


def compose_mass_range(mz, ppm) -> str:
    """From `mz` (the exact mass) and `ppm` (mass_error_ppm), compute the window
    [mz - mz*ppm/1e6, mz + mz*ppm/1e6] and return it in the "min-max" string
    format that `DbDataProvider.GetChroData` expects.

    Precision is preserved with `Decimal` (12 decimal places in the database).
    The final string does NOT truncate digits - it returns the exact value
    derived from the Decimal arithmetic (the same format mass_range used to be
    stored in).
    """
    from decimal import Decimal, getcontext
    getcontext().prec = 40

    mz_dec = Decimal(str(mz))
    ppm_dec = Decimal(str(ppm))
    delta = mz_dec * ppm_dec / Decimal(1_000_000)
    lo = mz_dec - delta
    hi = mz_dec + delta
    # format(..., 'f') avoids scientific notation; Decimal keeps every digit.
    return f"{format(lo, 'f')}-{format(hi, 'f')}"


def _attach_mass_range(row: dict) -> dict:
    """Take a `method_substances` row (a dict from the database) carrying `mz`
    and `mass_error_ppm`, inject the computed `mass_range` key and return it.

    The rest of the pipeline (process_chrom_data, get_data_and_process and so
    on) keeps reading `substance["mass_range"]` exactly as before - the only
    difference is that the value is now computed at runtime.
    """
    if row is None:
        return row
    if "mz" in row and "mass_error_ppm" in row and row.get("mz") is not None:
        row["mass_range"] = compose_mass_range(row["mz"], row["mass_error_ppm"])
    return row


def get_data_and_process(raw_file_path, start_time, end_time, mass_range, scan_filter, smoothing_type, smoothing_value):
    """A thin wrapper over `raw.GetChroData(...)`.

    It takes a raw object (which can be a DbDataProvider OR a real MSFileReader
    - the interface is the same) plus a substance's parameters, and returns a
    DataFrame with the columns (retention_time, intensity).

    This is where the "MSFileReader-like" abstraction pays off: everything
    downstream works with pandas DataFrames, without caring where they came
    from.
    """
    result = raw_file_path.GetChroData(
        startTime=start_time,
        endTime=end_time,
        massRange1=mass_range,
        scanFilter=scan_filter,
        smoothingType=smoothing_type,
        smoothingValue=smoothing_value
    )
    return pd.DataFrame({'retention_time': result[0][0], 'intensity': result[0][1]})

_SEVERITY_RANK = {
    "SUSPECT_VERY_HIGH": 5,
    "SUSPECT_HIGH":      4,
    "SUSPECT_MED":       3,
    "SUSPECT_LOW":       2,
    "SUSPECT":           1,
    "NEGATIVE":          0,
    "INVALID":          -1,
}


def verify_two_peaks(
    corrected_sample: pd.DataFrame,
    corrected_cqcl: pd.DataFrame,
    s_name: str,
    rt_window: float = 0.10,
    min_rt_sep: float = 0.08,
    r2_threshold: float = -50.0,
    dtw_limit: float = 1.0,
    log_file: Optional[str] = None,
) -> str:
    """
    The `two_peaks` rule - for isomeric substances (same molecular formula with
    distinct fragments, e.g. '5-APDB / Mephedrone'). It evaluates BOTH CQCL
    peaks and returns the MOST SEVERE level of the two.

    For each peak:
      - int1 = max sample intensity in the +/-rt_window window (resistant to RT
        shift)
      - int2 = max cqcl intensity in the +/-rt_window window
      - Gate: passes_joint (R2 > threshold AND DTW < limit) OR int1 > int2
      - Classification:
          int1 > int2          -> SUSPECT_VERY_HIGH
          passes_joint         -> classify_joint_fuzzy(R2, DTW) (LOW/MED/HIGH)
          otherwise            -> NEGATIVE

    Aggregation: max(severity) across the 2 peaks. It used to be "the first
    suspect found" - this reflects the best evidence from both peaks instead.
    """

    # ---------- find the 2 peaks ----------
    if corrected_cqcl is None or corrected_cqcl.empty:
        return "NEGATIVE"

    rt = corrected_cqcl["retention_time"].to_numpy()
    y  = corrected_cqcl["intensity"].to_numpy()

    if len(y) < 5:
        return "NEGATIVE"

    idx = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    if len(idx) == 0:
        idx = np.argsort(y)[-2:]

    idx_sorted = idx[np.argsort(y[idx])[::-1]]

    peaks_rt = []
    for i in idx_sorted:
        if y[i] <= 0:
            continue
        if all(abs(rt[i] - p) >= min_rt_sep for p in peaks_rt):
            peaks_rt.append(float(rt[i]))
        if len(peaks_rt) == 2:
            break

    if not peaks_rt:
        msg = f"ℹ️ [{s_name}] two_peaks: no valid peak in the CQCL -> NEGATIVE"
        if log_file:
            log_message(msg, log_file)
        else:
            print(msg)
        return "NEGATIVE"

    # ---------- inner helper: evaluate 1 peak, classify its level ----------
    def eval_one_peak(rt_peak: float, idx_peak: int) -> str:
        sample_w, cqcl_w = windows_max_int_cqcl(corrected_sample, corrected_cqcl, rt_peak, rt_window=rt_window)

        if sample_w is None or cqcl_w is None or len(sample_w) < 2 or len(cqcl_w) < 2:
            msg = f"⚠️ [{s_name}] two_peaks#{idx_peak} peak@{rt_peak:.4f} window too small -> NEGATIVE"
            if log_file:
                log_message(msg, log_file)
            else:
                print(msg)
            return "NEGATIVE"

        r_squared_float, r_squared_str = calculate_r_squared(sample_w, cqcl_w)
        if r_squared_float is None or r_squared_float in (0.0, 1.0):
            msg = f"⚠️ [{s_name}] two_peaks#{idx_peak} peak@{rt_peak:.4f} invalid R2 ({r_squared_str}) -> NEGATIVE"
            if log_file:
                log_message(msg, log_file)
            else:
                print(msg)
            return "NEGATIVE"

        # int1/int2 as the window max (consistent with evaluate_against_reference_fast)
        int1 = float(sample_w["intensity"].max()) if not sample_w.empty else 0.0
        int2 = float(cqcl_w["intensity"].max())   if not cqcl_w.empty   else 0.0

        # DTW
        A_y = sample_w["intensity"].values
        C_y = cqcl_w["intensity"].values
        A_y_norm = (A_y - np.min(A_y)) / (np.max(A_y) - np.min(A_y) + 1e-9)
        C_y_norm = (C_y - np.min(C_y)) / (np.max(C_y) - np.min(C_y) + 1e-9)
        dtw_dist = float(dtw.distance(A_y_norm, C_y_norm))

        passes_joint = (r_squared_float > r2_threshold and dtw_dist < dtw_limit)

        # quality (for the log + the fuzzy step)
        q_r2  = (r_squared_float - r2_threshold) / (1.0 - r2_threshold) if (1.0 - r2_threshold) != 0 else 0.0
        q_dtw = 1.0 - dtw_dist / dtw_limit if dtw_limit > 0 else 0.0
        int_cmp = ">" if int1 > int2 else "≤"

        log_prefix = (
            f"[{s_name}] two_peaks#{idx_peak} peak@{rt_peak:.4f} "
            f"sample={int1:.2e} cqcl={int2:.2e} (int1{int_cmp}int2) | "
            f"R²={r_squared_float:.3f} (q={q_r2:.3f}) | "
            f"DTW={dtw_dist:.4f} (q={q_dtw:.3f})"
        )

        if not (int1 > int2 or passes_joint):
            msg = f"ℹ️ {log_prefix} | gate=FAIL → NEGATIVE"
            if log_file:
                log_message(msg, log_file)
            else:
                print(msg)
            return "NEGATIVE"

        if int1 > int2:
            result = "SUSPECT_VERY_HIGH"
        else:
            result = classify_joint_fuzzy(r_squared_float, r2_threshold, dtw_dist, dtw_limit)

        msg = f"ℹ️ {log_prefix} | gate=PASS → {result}"
        if log_file:
            log_message(msg, log_file)
        else:
            print(msg)
        return result

    # ---------- evaluate BOTH peaks, return the most severe ----------
    # Change of behaviour: it used to be "the first SUSPECT found"; now it takes
    # the highest level of the 2. That stops a HIGH peak from being swallowed by
    # a MED one, and gives full visibility in the log.
    results = [eval_one_peak(rt_p, i + 1) for i, rt_p in enumerate(peaks_rt)]
    best = max(results, key=lambda r: _SEVERITY_RANK.get(r, -2))
    return best

def any_points_suspect_fixed_window_on_control_peak(
    corrected_sample: pd.DataFrame,
    corrected_control: pd.DataFrame,
    *,
    rt_half_window: float = 0.1,  # fixed window
    min_points: int = 5,
    ratio_limit: float = 100.0,     # control/sample <= 100 => SUSPECT
    min_intensity_sample: float = 1.0,   # keeps noise out
    min_intensity_control: float = 1.0,  # keeps a "flat" control out
    pad_rt: float = 0.0                 # optional: extra slack
):
    """
    A fixed window based on the RT of the CONTROL's largest peak:
      rt_start = peak_rt - rt_half_window
      rt_end   = peak_rt + rt_half_window

    Inside that window:
      sample_max  = max(sample intensity)
      control_max = max(control intensity)

    The rule:
      ratio = control_max / sample_max
      if ratio <= ratio_limit => SUSPECT_LOW
      otherwise => NEGATIVE
    """

    debug = {
        "peak_rt_control": None,
        "peak_int_control": None,
        "rt_start": None,
        "rt_end": None,
        "sample_max_in_window": None,
        "control_max_in_window": None,
        "ratio_control_over_sample": None,
        "ratio_limit": ratio_limit,
        "sample_apex_in_window": None,
        "reason": None
    }

    if corrected_control is None or corrected_control.empty:
        debug["reason"] = "empty_control"
        return "INVALID", (None, None), debug

    if corrected_sample is None or corrected_sample.empty:
        debug["reason"] = "empty_sample"
        return "INVALID", (None, None), debug

    control = corrected_control.sort_values("retention_time").reset_index(drop=True)
    sample  = corrected_sample.sort_values("retention_time").reset_index(drop=True)

    if len(control) < min_points:
        debug["reason"] = "control_too_few_points"
        return "INVALID", (None, None), debug

    # the control's tallest peak (it defines the centre RT)
    y_c = pd.to_numeric(control["intensity"], errors="coerce").fillna(0).to_numpy(dtype=float)
    rt_c = pd.to_numeric(control["retention_time"], errors="coerce").to_numpy(dtype=float)

    peak_idx = int(np.argmax(y_c))
    peak_int_control = float(y_c[peak_idx])
    peak_rt_control = float(rt_c[peak_idx])

    debug["peak_rt_control"] = peak_rt_control
    debug["peak_int_control"] = peak_int_control

    if not np.isfinite(peak_int_control) or peak_int_control < float(min_intensity_control):
        debug["reason"] = "control_peak_below_min_intensity"
        return "INVALID", (None, None), debug

    # Fixed window around the RT of the control's peak
    rt_start = peak_rt_control - float(rt_half_window) - float(pad_rt)
    rt_end   = peak_rt_control + float(rt_half_window) + float(pad_rt)

    debug["rt_start"] = rt_start
    debug["rt_end"] = rt_end

    cw = control[(control["retention_time"] >= rt_start) & (control["retention_time"] <= rt_end)]
    sw = sample[(sample["retention_time"] >= rt_start) & (sample["retention_time"] <= rt_end)]

    if cw.empty:
        debug["reason"] = "no_control_points_in_window"
        return "INVALID", (rt_start, rt_end), debug

    control_max = float(pd.to_numeric(cw["intensity"], errors="coerce").fillna(0).max())
    debug["control_max_in_window"] = control_max

    if not np.isfinite(control_max) or control_max < float(min_intensity_control):
        debug["reason"] = "control_no_signal_in_window"
        return "INVALID", (rt_start, rt_end), debug

    if sw.empty:
        debug["reason"] = "no_sample_points_in_window"
        return "NEGATIVE", (rt_start, rt_end), debug

    # Pre-check: the sample's peak (apex/local max) has to fall inside the
    # control's window. Otherwise what shows up in the window is just the
    # shoulder of another peak, or an interferent.
    y_s = pd.to_numeric(sample["intensity"], errors="coerce").fillna(0).to_numpy(dtype=float)
    rt_s = pd.to_numeric(sample["retention_time"], errors="coerce").to_numpy(dtype=float)

    apex_in_window = False
    if len(y_s) >= 3:
        is_local_max = (y_s[1:-1] > y_s[:-2]) & (y_s[1:-1] > y_s[2:])
        apex_idx = np.where(is_local_max)[0] + 1
        apex_rts = rt_s[apex_idx]
        apex_intensities = y_s[apex_idx]

        # keep only non-trivial apexes (intensity above the minimum)
        valid = apex_intensities >= float(min_intensity_sample)
        apex_rts = apex_rts[valid]

        apex_in_window = bool(np.any((apex_rts >= rt_start) & (apex_rts <= rt_end)))

    debug["sample_apex_in_window"] = apex_in_window

    if not apex_in_window:
        debug["reason"] = "no_sample_apex_in_window"
        return "NEGATIVE", (rt_start, rt_end), debug

    sample_max = float(pd.to_numeric(sw["intensity"], errors="coerce").fillna(0).max())
    debug["sample_max_in_window"] = sample_max

    if not np.isfinite(sample_max) or sample_max < float(min_intensity_sample):
        debug["reason"] = "sample_no_signal_in_window"
        return "NEGATIVE", (rt_start, rt_end), debug

    # The rule: a sample up to 100x smaller => control/sample <= 100
    ratio = control_max / sample_max
    debug["ratio_control_over_sample"] = float(ratio)

    if ratio <= float(ratio_limit):
        debug["reason"] = "ratio<=limit => SUSPECT"
        return "SUSPECT_LOW", (rt_start, rt_end), debug

    debug["reason"] = "ratio>limit => NEGATIVE"
    return "NEGATIVE", (rt_start, rt_end), debug

def is_suspect_label(label: str) -> bool:
    return isinstance(label, str) and "SUSPECT" in label.upper()

def calculate_r_squared(sample, cqcl):
    try:
        min_length = min(len(sample), len(cqcl))
        sample = sample.iloc[:min_length]
        cqcl = cqcl.iloc[:min_length]

        r_squared = r2_score(sample.iloc[:, 1], cqcl.iloc[:, 1])
        r_squared_str = "{:.6f}".format(r_squared)[:6]
        r_squared_float = float(r_squared_str)

        return r_squared_float, r_squared_str
    except Exception as e:
        print(f"Error calculating r_squared: {e}")
        return None, None

def classify_joint_fuzzy(r2: float, r2_threshold: float, dtw_val: float, dtw_limit: float) -> str:
    """Classify into SUSPECT_LOW/MED/HIGH using Mamdani fuzzy logic over R2 and DTW.

    Precondition: passes_joint == True (R2 > r2_threshold AND DTW < dtw_limit).
    This function ONLY classifies the LEVEL of suspicion - it never returns
    NEGATIVE. The suspect-vs-negative decision belongs to the GATE, over in
    `evaluate_against_reference_fast`.

    Strategy:
      1. Normalise R2 and DTW into a "quality in [0, 1]" relative to the
         substance's own threshold:
            quality_r2  = (R2 - r2_threshold) / (1 - r2_threshold)
            quality_dtw = 1 - DTW / dtw_limit
         That way the scale adapts automatically to whatever calibration the
         user picks (strict or loose) - quality=0 means "right at the gate" and
         quality=1 means "perfect match".
      2. Apply 3 fuzzy memberships (LOW/MED/HIGH) to each quality.
      3. Run 9 Mamdani rules (3x3) with OPTIMISTIC compensation - a strong
         indicator rescues a weak one (R2 HIGH x DTW LOW -> MED, and vice versa).
      4. Winner takes all: the category with the highest aggregated strength wins.
    """
    # Normalisation: the gate already guaranteed both qualities are > 0
    q_r2  = (r2 - r2_threshold) / (1.0 - r2_threshold)
    q_dtw = 1.0 - dtw_val / dtw_limit
    # defensive clamp (the gate already filtered, but numbers drift)
    q_r2  = max(0.0, min(1.0, q_r2))
    q_dtw = max(0.0, min(1.0, q_dtw))

    def _trap(x, a, b, c, d):
        if x <= a or x >= d: return 0.0
        if b <= x <= c: return 1.0
        if x < b: return (x - a) / (b - a)
        return (d - x) / (d - c)

    def _tri(x, a, b, c):
        if x <= a or x >= c: return 0.0
        if x == b: return 1.0
        if x < b: return (x - a) / (b - a)
        return (c - x) / (c - b)

    # Memberships over quality in [0, 1]
    def _mfs(q):
        return {
            "LOW":  _trap(q, 0.0, 0.0, 0.25, 0.40),
            "MED":  _tri (q, 0.30, 0.50, 0.70),
            "HIGH": _trap(q, 0.60, 0.75, 1.00, 1.00),
        }

    mf_r2 = _mfs(q_r2)
    mf_dtw = _mfs(q_dtw)

    # The 3x3 rule table (optimistic compensation - a strong indicator rescues
    # a weak one)
    #              DTW: HIGH    MED     LOW
    # R2 HIGH:          HIGH    HIGH    MED
    # R2 MED:           HIGH    MED     LOW
    # R2 LOW:           MED     LOW     LOW
    rules = [
        ("HIGH", "HIGH", "HIGH"),
        ("HIGH", "MED",  "HIGH"),
        ("HIGH", "LOW",  "MED"),
        ("MED",  "HIGH", "HIGH"),
        ("MED",  "MED",  "MED"),
        ("MED",  "LOW",  "LOW"),
        ("LOW",  "HIGH", "MED"),
        ("LOW",  "MED",  "LOW"),
        ("LOW",  "LOW",  "LOW"),
    ]

    scores = {"LOW": 0.0, "MED": 0.0, "HIGH": 0.0}
    for level_r2, level_dtw, output in rules:
        strength = min(mf_r2[level_r2], mf_dtw[level_dtw])
        if strength > scores[output]:
            scores[output] = strength

    # Winner-takes-all
    winner = max(scores, key=scores.get)
    return f"SUSPECT_{winner}"


def classify_by_cutoff_zones(est_conc: float, cut_off: float, fortified: float) -> str:
    """Classify into SUSPECT_LOW/MED/HIGH/VERY_HIGH by zones inside the
    [cut_off, fortified] interval.

    Precondition: est_conc >= cut_off (the suspect-vs-negative decision belongs
    to `evaluate_cutoff`, not to this function).

    Zones (t = the relative position inside [cut_off, fortified]):
        est_conc > fortified         -> SUSPECT_VERY_HIGH (above the fortified level)
        t in [0.50, 1.00]            -> SUSPECT_HIGH
        t in [0.25, 0.50)            -> SUSPECT_MED
        t in [0.00, 0.25)            -> SUSPECT_LOW

    The scale adjusts per substance automatically - the range is defined by the
    calibration itself (cut_off and fortified come from the database).
    """
    if est_conc > fortified:
        return "SUSPECT_VERY_HIGH"

    span = fortified - cut_off
    if span <= 0:
        # Degenerate configuration (fortified <= cut_off). It should not happen,
        # since the database pairs them. Defensive: classify as HIGH if it got
        # past the cut_off.
        return "SUSPECT_HIGH"

    t = (est_conc - cut_off) / span
    if t < 0.25:
        return "SUSPECT_LOW"
    if t < 0.50:
        return "SUSPECT_MED"
    return "SUSPECT_HIGH"


def evaluate_cutoff(cursor, raw_sample, raw_cqcl, substance, data_sample, data_cqcl, rt2):
    """
    Run the complete cut-off logic and return:
      (result: str, concentration: float|None, errors: list[str], debug: dict)

    Notes:
    - It keeps the original rule of using windows WITHOUT the CQN correction
      (data_sample/data_cqcl).
    - It looks the internal standard up through the cutoff_processing_groups table.
    - `debug` carries the 4 areas + estConc + cutoff + fortified for the enriched log.
    """
    errors = []
    concentration = None
    debug = {"areas": {}, "cut_off": None, "fortified": None, "est_conc": None, "zero_area": None}

    try:
        cutoff_value = float(substance.get("cut_off") or 0)
        fortified = float(substance.get("fortified_concentration") or 0)
        debug["cut_off"] = cutoff_value
        debug["fortified"] = fortified

        # The substance window (uncorrected, as in the original)
        sample_window, cqcl_window = windows_max_int_cqcl(data_sample, data_cqcl, rt2)

        # Look up the cut-off's internal standard
        cursor.execute("""
            SELECT 
                cutoff_processing_groups.id AS cutoff_group_id,
                cutoff_processing_groups.method_id AS cutoff_group_method_id,
                cutoff_processing_groups.internal_standard AS cutoff_group_internal_standard_id,
                method_substances.id AS internal_standard,
                method_substances.name AS internal_standard_name,
                method_substances.start_time AS internal_standard_start_time,
                method_substances.end_time AS internal_standard_end_time,
                method_substances.mz AS internal_standard_mz,
                method_substances.mass_error_ppm AS internal_standard_mass_error_ppm,
                sf.name AS internal_standard_scan_filter,
                method_substances.smoothing_type AS internal_standard_smoothing_type,
                method_substances.smoothing_value AS internal_standard_smoothing_value
            FROM cutoff_processing_groups
            INNER JOIN method_substances
                ON cutoff_processing_groups.internal_standard = method_substances.id
            LEFT JOIN scan_filters sf
                ON sf.id = method_substances.scan_filter_id
            WHERE cutoff_processing_groups.method_id = %s
              AND method_substances.deleted_at IS NULL
            LIMIT 1
        """, (substance["method_id"],))

        cutoff_group = cursor.fetchone()
        if not cutoff_group:
            errors.append("No internal standard defined for this method's cut-off.")
            return "NEGATIVE", None, errors, debug

        internal_std = {
            "start_time": cutoff_group["internal_standard_start_time"],
            "end_time": cutoff_group["internal_standard_end_time"],
            "mass_range": compose_mass_range(
                cutoff_group["internal_standard_mz"],
                cutoff_group["internal_standard_mass_error_ppm"],
            ),
            "scan_filter": cutoff_group["internal_standard_scan_filter"],
            "smoothing_type": cutoff_group["internal_standard_smoothing_type"],
            "smoothing_value": cutoff_group["internal_standard_smoothing_value"],
        }

        # Local helper (the same logic as process_chrom_data)
        def process_chrom_data(raw, sub):
            return get_data_and_process(
                raw,
                sub["start_time"],
                sub["end_time"],
                sub["mass_range"],
                sub["scan_filter"],
                sub["smoothing_type"],
                sub["smoothing_value"]
            )

        data_internal_cqcl = process_chrom_data(raw_cqcl, internal_std)
        data_internal_sample = process_chrom_data(raw_sample, internal_std)

        concentration, zero_area, areas = cut_off_verification(
            sample_window, cqcl_window,
            data_internal_sample, data_internal_cqcl,
            fortified
        )
        debug["areas"] = areas
        debug["zero_area"] = zero_area
        debug["est_conc"] = concentration

        if concentration is None:
            errors.append(f"cut_off_verification: zero area in '{zero_area}'")
            return "NEGATIVE", None, errors, debug

        if concentration < cutoff_value:
            return "NEGATIVE", concentration, errors, debug

        result = classify_by_cutoff_zones(concentration, cutoff_value, fortified)
        return result, concentration, errors, debug

    except Exception as e:
        errors.append(f"Cut-off computation failed: {e}")
        return "NEGATIVE", None, errors, debug

def evaluate_against_reference_fast(
    cursor,
    raw_sample,
    raw_reference,
    substance,
    corrected_sample,
    corrected_reference,
    data_sample,
    data_reference,
    log_file,
    dtw_limit=1.0,
):
    errors = []
    concentration = None
    s_name = substance["name"]

    # R2 threshold (NULL in the database -> fallback of -50).
    r2_threshold_raw = substance.get("r2_threshold")
    r2_threshold = float(r2_threshold_raw) if r2_threshold_raw is not None else -50.0

    # The tallest peak in the corrected reference
    i2 = corrected_reference["intensity"].idxmax()
    rt2 = float(corrected_reference.loc[i2, "retention_time"])
    int2 = float(corrected_reference.loc[i2, "intensity"])

    # Windows around the peak
    sample_window_df, ref_window_df = windows_max_int_cqcl(corrected_sample, corrected_reference, rt2)

    # int1 = the TALLEST sample point INSIDE the window (rt2 +/- rt_window).
    # It used to be the sample intensity at the CQCL's exact RT - which punished
    # cases where the sample had a small RT shift relative to the control (even
    # with a visibly taller peak, int1 landed on the peak's downslope and came
    # out below int2). The +/-0.1 min window absorbs the natural shift without
    # loosening the criterion (it still only looks at the substance's
    # chromatographic region).
    if sample_window_df.empty:
        int1 = 0.0
    else:
        int1 = float(sample_window_df["intensity"].max())

    if len(sample_window_df) < 2 or len(ref_window_df) < 2:
        print(f"⚠️ [{s_name}] Invalid window ({sample_window_df}) - ignoring.")
        return "INVALID", None, [], {
            "ref": s_name, "reason": "insufficient_window",
            "rt2": rt2, "int1": int1, "int2": int2, "r2": None, "dtw": None
        }

    # R²
    r_squared_float, r_squared_str = calculate_r_squared(sample_window_df, ref_window_df)
    if r_squared_float is None or r_squared_float in (0.0, 1.0):
        log_message(f"⚠️ [{s_name}] Invalid R2 ({r_squared_str}) - ignoring.", log_file)
        return "INVALID", None, [], {
            "ref": s_name, "reason": f"invalid_r2({r_squared_str})",
            "rt2": rt2, "int1": int1, "int2": int2, "r2": r_squared_str, "dtw": None
        }

    # DTW
    A_y = sample_window_df["intensity"].values
    C_y = ref_window_df["intensity"].values
    A_y_norm = (A_y - np.min(A_y)) / (np.max(A_y) - np.min(A_y) + 1e-9)
    C_y_norm = (C_y - np.min(C_y)) / (np.max(C_y) - np.min(C_y) + 1e-9)
    dtw_dist = float(dtw.distance(A_y_norm, C_y_norm))

    # Gate: SUSPECT if the sample peak > the control's, OR if R2 and DTW are both OK
    passes_joint = (r_squared_float > r2_threshold and dtw_dist < dtw_limit)

    # quality_r2/quality_dtw are values normalised to [0, 1] relative to the
    # substance's own threshold. Used by the fuzzy step below and logged for
    # debugging.
    q_r2  = (r_squared_float - r2_threshold) / (1.0 - r2_threshold) if (1.0 - r2_threshold) != 0 else 0.0
    q_dtw = 1.0 - dtw_dist / dtw_limit if dtw_limit > 0 else 0.0

    # Pre-build the shared part of the enriched log line
    has_cutoff = substance.get("cut_off") is not None
    int_cmp = ">" if int1 > int2 else "≤"
    log_prefix = (
        f"[{s_name}] sample={int1:.2e} cqcl={int2:.2e} (int1{int_cmp}int2) | "
        f"R²={r_squared_float:.3f} (q={q_r2:.3f}) | "
        f"DTW={dtw_dist:.4f} (q={q_dtw:.3f})"
    )

    if not (int1 > int2 or passes_joint):
        log_message(f"ℹ️ {log_prefix} | gate=FAIL → NEGATIVE", log_file)
        return "NEGATIVE", None, [], {
            "ref": s_name, "reason": "gate_failed",
            "rt2": rt2, "int1": int1, "int2": int2, "r2": r_squared_float, "dtw": dtw_dist
        }

    # Cut-off?
    if has_cutoff:
        # the cut-off uses data_sample/data_reference UNCORRECTED, but with the same rt2
        result, concentration, cutoff_errors, cutoff_debug = evaluate_cutoff(
            cursor=cursor,
            raw_sample=raw_sample,
            raw_cqcl=raw_reference,
            substance=substance,
            data_sample=data_sample,
            data_cqcl=data_reference,
            rt2=rt2
        )
        errors += cutoff_errors

        # Enriched log line for the cut-off path
        co_val = cutoff_debug.get("cut_off")
        fort = cutoff_debug.get("fortified")
        est_conc = cutoff_debug.get("est_conc")
        zero_area = cutoff_debug.get("zero_area")
        areas = cutoff_debug.get("areas") or {}
        areas_str = " ".join(
            f"{k}={(v if v is None else f'{v:.2e}')}" + ("❌" if k == zero_area else "")
            for k, v in areas.items()
        ) if areas else "n/a"
        conc_str = f"estConc={est_conc:.3f}" if est_conc is not None else "estConc=None"

        if zero_area is not None:
            log_message(
                f"⚠️ {log_prefix} | gate=PASS | cut_off={co_val} fortified={fort} "
                f"areas[{areas_str}] → cut_off_verification returned None "
                f"(zero area in: {zero_area}) → {result}",
                log_file
            )
        else:
            log_message(
                f"ℹ️ {log_prefix} | gate=PASS | cut_off={co_val} fortified={fort} "
                f"areas[{areas_str}] → {conc_str} → {result}",
                log_file
            )

        return result, concentration, errors, {
            "ref": s_name, "reason": "cutoff",
            "rt2": rt2, "int1": int1, "int2": int2, "r2": r_squared_float, "dtw": dtw_dist
        }

    # ==========================================
    # No cut-off: binary classification + fuzzy
    # ==========================================
    # int1 > int2 = the sample is physically taller than the control -> the
    # strongest possible signal; it becomes VERY_HIGH (only for substances
    # without a cut_off; the ones that have it took the branch above and are
    # classified by concentration).
    # passes_joint = R2 and DTW both cleared their thresholds -> it goes into the
    # fuzzy step, which classifies the LEVEL of suspicion (LOW/MED/HIGH).
    if int1 > int2:
        result = "SUSPECT_VERY_HIGH"
    elif passes_joint:
        result = classify_joint_fuzzy(r_squared_float, r2_threshold, dtw_dist, dtw_limit)
    else:
        result = "NEGATIVE"

    log_message(f"ℹ️ {log_prefix} | gate=PASS → {result}", log_file)

    return result, None, errors, {
            "ref": s_name, "reason": "joint",
            "rt2": rt2, "int1": int1, "int2": int2, "r2": r_squared_float, "dtw": dtw_dist
    }

def adjust_intensity(sample_df, cqn_df):
    """Subtract the baseline (CQN) from the target chromatogram.

    The CQN is a negative control (clean matrix, no analyte). Whatever shows up
    in it is the method's background: solvent contaminants, column background,
    satellite peaks of abundant species, and so on.

    For every point of the sample/cqcl/cqcl_reinj we interpolate the CQN at the
    same RT and subtract it. Clamping at 0 avoids negative intensity (it has no
    physical meaning - going negative just means the noise wobbled).

    This correction happens BEFORE the R2/DTW step. The XICs stored in
    processing_sampleparameter are kept UNCORRECTED (raw) - the correction is
    applied again in the front end when it needs to be displayed. That is what
    guarantees the UI and the pipeline compute over the same data.

    Edge cases:
      - empty sample -> returned as-is (nothing to correct).
      - empty cqn -> returns a copy of the sample (with no CQN there is no
        baseline to subtract). It happens when the substance has no point in the
        CQN - common when the lake only stores centroids that carry signal.
    """
    if sample_df is None or sample_df.empty:
        return sample_df
    if cqn_df is None or cqn_df.empty:
        return sample_df.copy()

    cqn_interp = np.interp(sample_df['retention_time'], cqn_df['retention_time'], cqn_df['intensity'])
    corrected_intensity = sample_df['intensity'] - cqn_interp
    corrected_intensity[corrected_intensity < 0] = 0
    corrected_sample_df = sample_df.copy()
    corrected_sample_df['intensity'] = corrected_intensity
    return corrected_sample_df

def cut_off_verification(sample_substance_window,cqcl_substance_window,sample_internal_window,cqcl_internal_window,known_concentration_cqcl):
    def get_peak_area(window_df):
        def detect_peak_limits(rt_array, intensity_array, threshold_ratio=0.05):
            peak_max = np.max(intensity_array)
            if peak_max == 0 or np.all(intensity_array < peak_max * threshold_ratio):
                return 0, len(intensity_array)

            peak_center = np.argmax(intensity_array)
            threshold = peak_max * threshold_ratio

            start = peak_center
            while start > 0 and intensity_array[start] > threshold:
                start -= 1

            end = peak_center
            while end < len(intensity_array) - 1 and intensity_array[end] > threshold:
                end += 1

            return start, end

        rt = pd.to_numeric(window_df['retention_time'], errors='coerce').dropna().to_numpy()
        intensity = pd.to_numeric(window_df['intensity'], errors='coerce').dropna().to_numpy()

        if len(rt) == 0 or len(intensity) == 0:
            return 0.0

        start, end = detect_peak_limits(rt, intensity)
        # np.trapezoid: numpy 2.x removed np.trapz.
        return np.trapezoid(intensity[start:end], rt[start:end])

    try:
        area_sample_substance = get_peak_area(sample_substance_window)
        area_cqcl_substance = get_peak_area(cqcl_substance_window)
        area_sample_internal = get_peak_area(sample_internal_window)
        area_cqcl_internal = get_peak_area(cqcl_internal_window)

        # The return says WHICH area failed (for the log). It used to return
        # only None, and the operator had no idea where it broke.
        if area_sample_internal == 0:
            return None, "sample_internal", {
                "sample_substance": area_sample_substance,
                "cqcl_substance": area_cqcl_substance,
                "sample_internal": area_sample_internal,
                "cqcl_internal": area_cqcl_internal,
            }
        if area_cqcl_internal == 0:
            return None, "cqcl_internal", {
                "sample_substance": area_sample_substance,
                "cqcl_substance": area_cqcl_substance,
                "sample_internal": area_sample_internal,
                "cqcl_internal": area_cqcl_internal,
            }
        if area_cqcl_substance == 0:
            return None, "cqcl_substance", {
                "sample_substance": area_sample_substance,
                "cqcl_substance": area_cqcl_substance,
                "sample_internal": area_sample_internal,
                "cqcl_internal": area_cqcl_internal,
            }

        sample_ratio = area_sample_substance / area_sample_internal
        cqcl_ratio = area_cqcl_substance / area_cqcl_internal

        # Cast to a native float: the areas come out of np.trapezoid() as
        # numpy.float64. Without the cast the INSERT fails with
        #   "Python 'float64' cannot be converted to a MySQL type"
        # (the driver does not convert numpy types automatically).
        estimated_concentration = float(
            (sample_ratio * known_concentration_cqcl) / cqcl_ratio
        )
        return estimated_concentration, None, {
            "sample_substance": float(area_sample_substance),
            "cqcl_substance": float(area_cqcl_substance),
            "sample_internal": float(area_sample_internal),
            "cqcl_internal": float(area_cqcl_internal),
        }

    except Exception as error:
        return None, f"exception:{error}", {}

def windows_max_int_cqcl(sample, cqcl, rt2, rt_window = 0.1):
    """Slice the [rt2 - rt_window, rt2 + rt_window] window out of both DataFrames.

    `rt2` is normally the RT of the CQCL's tallest peak (the reference). Slicing
    the same window out of sample and cqcl produces two synchronised
    chromatographic fragments that feed R2 and DTW.

    rt_window=0.1 min (= +/-6 seconds) is the default calibrated against the
    chromatographic peaks typical of the TVII method. For a wider peak, cutting
    more aggressively would throw signal away.
    """
    sample_window_df = sample[(sample['retention_time'] >= rt2 - rt_window) & (sample['retention_time'] <= rt2 + rt_window)]
    cqcl_window_df = cqcl[(cqcl['retention_time'] >= rt2 - rt_window) & (cqcl['retention_time'] <= rt2 + rt_window)]

    return sample_window_df, cqcl_window_df

def insert_sample_result(cursor, task_id, sample_name, substance_name, result, concentration_value):
    """
    Insert one row into the sample_results table.
    """
    errors = []
    now = datetime.utcnow()

    # Defence in depth: make sure concentration_value is a native Python float
    # (and not a numpy.float64). The driver does not convert numpy types
    # automatically - without this cast the INSERT would fail with
    #   "Python 'float64' cannot be converted to a MySQL type"
    # cut_off_verification already casts at the source, but this covers any
    # other path that may populate the field in the future.
    if concentration_value is not None:
        try:
            concentration_value = float(concentration_value)
        except (TypeError, ValueError):
            concentration_value = None

    try:
        cursor.execute("""
            INSERT INTO sample_results (
                sample_processing_task_id,
                sample_name,
                substance,
                result,
                concentration_value,
                created_at,
                updated_at
            ) VALUES (%s, %s, %s, %s, %s, %s, %s)
        """, (
            task_id,
            sample_name,
            substance_name,
            result,
            concentration_value,
            now,
            now
        ))

        sample_result_id = cursor.lastrowid

    except sqlite3.Error as e:
        tb = traceback.format_exc()
        errors.append(
            f"[SampleResult] Database error while inserting '{substance_name}': {str(e)}\n{tb}"
        )
        return None, errors

    return sample_result_id, errors

def insert_sample_parameters(
    sample_result_id, sample_type, data_sample, substance_name,
    max_retries=5, retry_delay=10,
    ch_client=None,
):
    """
    Insert chromatographic parameters into ClickHouse.
    On a memory error (MEMORY_LIMIT_EXCEEDED, code 241) it retries several times
    before giving up.

    `ch_client` (CONNECTION REUSE):
        - When passed: used directly - it does NOT open a new socket and does
          NOT close it at the end. Indispensable when this function is called
          hundreds or thousands of times in a loop (the real pipeline: N
          substances * N samples * 3-4 types per substance). Without it, Windows
          runs out of ephemeral ports (WinError 10048/10055).
        - When None: legacy behaviour - a new connection per call, closed at the
          end (only viable at low volume: ad-hoc scripts, tests, and so on).

    The memory retry works the same in both modes.
    """
    errors = []
    timestamp_now = datetime.utcnow()
    database = config.CH_ANALYSIS_DB

    try:
        # Make sure the id is usable
        if isinstance(sample_result_id, (tuple, list)):
            sample_result_id = sample_result_id[0]
        if sample_result_id is None:
            raise ValueError("sample_result_id cannot be None")

        sample_result_id = int(sample_result_id)

        # Clean and convert the data
        data_sample = data_sample.dropna(subset=["retention_time", "intensity"])
        data_sample = data_sample.astype({
            "retention_time": "float32",
            "intensity": "float32"
        })

        if data_sample.empty:
            return errors

        # Build every row at once (as in the original)
        rows = [
            (
                str(uuid.uuid4()),
                sample_result_id,
                sample_type,
                timestamp_now,
                float(row.retention_time),
                float(row.intensity),
            )
            for row in data_sample.itertuples(index=False)
        ]

        # Insert with automatic retries. The connection is reused when
        # `ch_client` came from the caller (the pipeline's hot path); a local one
        # is only created in the fallback (ad-hoc use).
        attempt = 1
        # Client ownership: if we created it here, we close it at the end. If it
        # came from the caller, we do NOT close it (the caller owns the
        # lifecycle).
        own_ch = ch_client is None
        ch = ch_client
        while attempt <= max_retries:
            try:
                if ch is None:
                    ch = click_conn()
                ch.insert(
                    f"{database}.processing_sampleparameter",
                    rows,
                    column_names=[
                        "id",
                        "sample_result_id",
                        "sample_type",
                        "timestamp",
                        "retention_time",
                        "intensity",
                    ],
                )
                break  # success -> leave the loop

            except Exception as e:
                err_text = str(e)
                tb = traceback.format_exc()

                # Detect a ClickHouse memory error (code 241)
                if "MEMORY_LIMIT_EXCEEDED" in err_text or "Code: 241" in err_text:
                    print(f"⚠️ Attempt {attempt} hit the memory limit. Waiting {retry_delay}s and trying again...")
                    time.sleep(retry_delay)
                    retry_delay = min(retry_delay * 2, 120)  # exponential backoff, up to 2 minutes
                    attempt += 1
                    continue  # try again

                # Any other error -> record it and stop
                errors.append(
                    f"[SampleParameter] Unrecoverable ClickHouse error "
                    f"for '{substance_name}': {err_text}\n{tb}"
                )
                break

        # If every attempt failed
        if attempt > max_retries:
            msg = f"[SampleParameter] Failed after {max_retries} attempts (persistent memory error) for '{substance_name}'."
            errors.append(msg)
            print(f"❌ {msg}")

    except Exception as e:
        tb = traceback.format_exc()
        errors.append(
            f"[SampleParameter] Unexpected failure for '{substance_name}': {str(e)}\n{tb}"
        )
    finally:
        # We only close the client when WE created it locally (legacy mode). On
        # the hot path the caller keeps it alive between calls.
        if 'own_ch' in locals() and own_ch and ch is not None:
            try:
                ch.close()
            except Exception:
                pass

    return errors

def _run_any_points_rule(ctx):
    # Per-substance ratio limit (NULL in the database -> fallback of 100.0)
    _ratio_raw = ctx["substance"].get("any_points_ratio_limit")
    ratio_limit = float(_ratio_raw) if _ratio_raw is not None else 100.0

    result, (_rt_start, _rt_end), dbg = any_points_suspect_fixed_window_on_control_peak(
        ctx["corrected_sample"],
        ctx["corrected_cqcl"],
        ratio_limit=ratio_limit,
    )
    if isinstance(result, str) and result.upper() == "INVALID":
        print(f"ℹ️ [{ctx['s_name']}] any_points: INVALID ({dbg.get('reason')}) -> skipping")
        return {"invalid": True, "errors": []}
    return {"result": result, "concentration": None, "errors": []}

def _run_two_peaks_rule(ctx):
    sub = ctx["substance"]
    r2_raw = sub.get("r2_threshold")
    dtw_raw = sub.get("dtw_limit")
    result = verify_two_peaks(
        ctx["corrected_sample"],
        ctx["corrected_cqcl"],
        ctx["s_name"],
        rt_window=0.10,
        min_rt_sep=0.08,
        r2_threshold=float(r2_raw) if r2_raw is not None else -50.0,
        dtw_limit=float(dtw_raw) if dtw_raw is not None else 1.0,
        log_file=ctx["log_file"],
    )
    return {"result": result, "concentration": None, "errors": []}

def _run_standard_rule(ctx):
    # DTW limit per substance (NULL in the database -> fallback of 1.0)
    _dtw_raw = ctx["substance"].get("dtw_limit")
    substance_dtw_limit = float(_dtw_raw) if _dtw_raw is not None else 1.0

    result, concentration, errs, _dbg = evaluate_against_reference_fast(
        ctx["cursor"], ctx["raw_sample"], ctx["raw_cqcl"], ctx["substance"],
        ctx["corrected_sample"], ctx["corrected_cqcl"],
        ctx["data_sample"], ctx["data_cqcl"], ctx["log_file"],
        dtw_limit=substance_dtw_limit,
    )

    if isinstance(result, str) and result.upper() == "INVALID":
        print(f"ℹ️ [{ctx['s_name']}] INVALID on the CQCL - skipping.")
        return {"invalid": True, "errors": errs}

    if result == "NEGATIVE":
        result2, concentration2, errs2, _dbg2 = evaluate_against_reference_fast(
            ctx["cursor"], ctx["raw_sample"], ctx["raw_cqcl_reinj"], ctx["substance"],
            ctx["corrected_sample"], ctx["corrected_reinj"],
            ctx["data_sample"], ctx["data_cqcl_reinj"], ctx["log_file"],
            dtw_limit=substance_dtw_limit,
        )

        if isinstance(result2, str) and result2.upper() == "INVALID":
            print(f"ℹ️ [{ctx['s_name']}] INVALID on the REINJ - skipping.")
            return {"invalid": True, "errors": errs2}

        if result2 != "NEGATIVE":
            result = result2
            concentration = concentration2

    return {"result": result, "concentration": concentration, "errors": []}

# ============================================================================
# Analysis rule dispatch.
#
# Every substance carries an `analysis_rule_code` in method_substances (through
# a JOIN with analysis_rules). This dict maps that code to the handler running
# the matching logic. Default = 'standard' (R2 + DTW + the reinjection fallback).
#
# To add a new rule:
#   1. Insert a row into `analysis_rules` with the new code.
#   2. Implement `_run_<rule>_rule(ctx)`, taking the ctx dict (assembled in
#      process_single_substance) and returning {result, concentration, errors}.
#   3. Register it here in the dict.
# ============================================================================
ANALYSIS_RULE_HANDLERS = {
    "standard": _run_standard_rule,
    "two_peaks": _run_two_peaks_rule,
    "any_points": _run_any_points_rule,
}

# The "do not interpret" rule (VISUALISATION-only substances).
#
# It is deliberately NOT in ANALYSIS_RULE_HANDLERS: handler dispatch happens
# AFTER the INVALID early returns (no peak in the CQCL, a window with fewer than
# 2 points), and in those cases the handler would never be called - the
# substance would come back INVALID having stored nothing. Since here we want
# the XIC stored regardless, the branch happens before them, in
# process_single_substance.
NO_ANALYZE_RULE_CODE = "no_analyze"

def process_single_substance(cursor, task_id, sample_path,raw_sample, raw_cqcl, raw_cqcl_reinj, raw_cqcn, substance, log_file, ch_client=None):
    """Process ONE substance for ONE sample. Called in a loop for every
    substance in the method.

    Steps:
      1. Extract the 4 raw chromatograms (sample, cqcl, cqcl_reinj, cqn) over
         the substance's window (RT range, m/z range, smoothing).
      2. Apply adjust_intensity (X - CQN) to the first 3 - that gives the
         "_corrected" ones.
      3. Check that the corrected CQCL has a peak (otherwise INVALID, early exit).
      4. Locate the RT of the CQCL peak (rt2) - the alignment reference.
      5. Dispatch to the appropriate analysis rule (standard / two_peaks /
         any_points) through ANALYSIS_RULE_HANDLERS.
      6. Persist the result:
            INSERT into sample_results (1 row in the analysis database)
            INSERT N points into processing_sampleparameter (analysis
                ClickHouse), in 4 categories (sample/cqcl/cqcl_reinj/cqcn) -
                RAW XICs, for reprocessing or displaying later.

    Exception: substances with analysis_rule = 'no_analyze' (VISUALISATION) run
    only steps 1-2 and 6 - they skip 3-5 entirely. They store
    result='NOT_INTERPRETED' and concentration=None. The CQN subtraction is
    applied at display time from the raw XICs stored here, exactly as in the
    normal flow.

    Returns: (sample_result_id, result, _, all_errors)
        sample_result_id: the id of the row inserted into sample_results, or None.
        result: 'NEGATIVE' | 'SUSPECT' | 'SUSPECT_LOW/MED/HIGH/VERY_HIGH'
                | 'INVALID' | 'NOT_INTERPRETED'.
    """
    try:
        concentration = None
        result_errors = []

        # === Helper for extracting chromatographic data ===
        # It returns the XIC with smoothing already applied, ready to analyse.
        def process_chrom_data(raw, sub):
            return get_data_and_process(
                raw,
                sub["start_time"],
                sub["end_time"],
                sub["mass_range"],
                sub["scan_filter"],
                sub["smoothing_type"],
                sub["smoothing_value"]
            )

        s_name = substance["name"]

        # === Extract the chromatographic data ===
        # 4 RAW XICs (no baseline subtraction), all over the same RT/mz window.
        # These are the ones stored in ClickHouse at the end, so the UI can
        # re-apply the correction when the user looks at the chromatograms.
        data_sample = process_chrom_data(raw_sample, substance)
        data_cqcl = process_chrom_data(raw_cqcl, substance)
        data_cqcl_reinj = process_chrom_data(raw_cqcl_reinj, substance)
        data_cqcn = process_chrom_data(raw_cqcn, substance)

        # === Correct the intensities with the CQN ===
        # Subtract the baseline. The "_corrected" ones are what feed the R2/DTW
        # computation (not the RAW ones). Only sample/cqcl/cqcl_reinj - the CQN
        # is the baseline itself, subtracting it from itself makes no sense.
        corrected_sample = adjust_intensity(data_sample, data_cqcn)
        corrected_cqcl = adjust_intensity(data_cqcl, data_cqcn)
        corrected_reinj  = adjust_intensity(data_cqcl_reinj, data_cqcn)

        rule_code = (substance.get("analysis_rule_code") or "standard").lower()

        if rule_code == NO_ANALYZE_RULE_CODE:
            # ================================================================
            # VISUALISATION MODE (analysis_rule = 'no_analyze')
            # ================================================================
            # A substance configured purely so the analyst can LOOK at the
            # chromatogram - there is no SUSPECT/NEGATIVE judgement and no
            # concentration computation.
            #
            # We deliberately skip EVERY interpretation step:
            #   - the CQCL peak check
            #   - the windows around the peak
            #   - the rule handler
            #
            # And, above all, we do NOT go through the INVALID early returns:
            # even with no peak in the CQCL we want to store whatever exists,
            # because the point is precisely to see the overlay (sample plus the
            # CQN-corrected controls). An "empty" XIC is information for the
            # analyst too. That is why the branch is here rather than a new
            # handler in ANALYSIS_RULE_HANDLERS (dispatch runs after the early
            # returns).
            #
            # Persistence (step 6) is the SAME - it falls into the common block
            # below.
            result = "NOT_INTERPRETED"
            concentration = None
        else:
            # No CQCL data -> there is no reference point for the window.
            # Common for substances with no signal in the control when the source
            # is the database (the lake ingests centroids, not profile data).
            if corrected_cqcl is None or corrected_cqcl.empty:
                print(f"⚠️ [{s_name}] CQCL has no points - ignoring.")
                return None, "INVALID", "0.0", []

            # === Find the tallest peak ===
            i2 = corrected_cqcl['intensity'].idxmax()
            rt2 = corrected_cqcl.loc[i2, 'retention_time']

            # === Define the windows around the peak ===
            sample_window_df, cqcl_window_df = windows_max_int_cqcl(corrected_sample, corrected_cqcl, rt2)

            if len(sample_window_df) < 2 or len(cqcl_window_df) < 2:
                print(f"⚠️ [{s_name}] Not enough data - ignoring.")
                return None, "INVALID", "0.0", []

            result = "NEGATIVE"
            concentration = None

            handler = ANALYSIS_RULE_HANDLERS.get(rule_code)

            if handler is None:
                print(f"⚠️ [{s_name}] analysis_rule '{rule_code}' has no registered handler - falling back to 'standard'")
                handler = _run_standard_rule

            ctx = {
                "cursor": cursor,
                "substance": substance,
                "s_name": s_name,
                "log_file": log_file,
                "raw_sample": raw_sample,
                "raw_cqcl": raw_cqcl,
                "raw_cqcl_reinj": raw_cqcl_reinj,
                "corrected_sample": corrected_sample,
                "corrected_cqcl": corrected_cqcl,
                "corrected_reinj": corrected_reinj,
                "data_sample": data_sample,
                "data_cqcl": data_cqcl,
                "data_cqcl_reinj": data_cqcl_reinj,
            }

            handler_result = handler(ctx)

            if handler_result.get("invalid"):
                return None, "INVALID", None, handler_result.get("errors", [])

            result = handler_result.get("result", "NEGATIVE")
            concentration = handler_result.get("concentration")
            result_errors += handler_result.get("errors", [])


        # === Persist the results ===
        sample_result_id, insert_errors = insert_sample_result(
            cursor, task_id, os.path.basename(sample_path),
            substance["name"], result, concentration
        )
        result_errors += insert_errors

        parameter_errors = []
        if sample_result_id:
            # `ch_client` is passed down to avoid opening a new socket on each
            # of the 4 inserts. Multiplied by N substances * N samples that
            # saves thousands of TCP connections per run.
            parameter_errors += insert_sample_parameters(sample_result_id, "sample", data_sample, substance["name"], ch_client=ch_client)
            parameter_errors += insert_sample_parameters(sample_result_id, "cqcl", data_cqcl, substance["name"], ch_client=ch_client)
            parameter_errors += insert_sample_parameters(sample_result_id, "cqcl_reinj", data_cqcl_reinj, substance["name"], ch_client=ch_client)
            parameter_errors += insert_sample_parameters(sample_result_id, "cqcn", data_cqcn, substance["name"], ch_client=ch_client)

        all_errors = result_errors + parameter_errors
        
        return sample_result_id, result, None, all_errors

    except Exception as e:
        tb = traceback.format_exc()
        print(f"❌ Error on {substance.get('name', 'unknown substance')}: {e}")
        return None, "INVALID", "0.0", [f"{e}\n{tb}"]

def log_message(message, log_file):
        """Write to the batch log.

        DELIBERATE DEVIATION FROM THE ORIGINAL: the write is guarded here. A log
        line must never take down an analysis, and one already did: the pipeline
        logs with emoji and Windows stdout assumes cp1252, which cannot encode
        them. `config` already reconfigures the streams to UTF-8; this try is the
        belt to that pair of braces.
        """
        try:
            print(message)
        except UnicodeEncodeError:
            print(message.encode("ascii", "replace").decode("ascii"))
        try:
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")
        except OSError:
            pass

# ============================================================================
# The adaptation that reads from the database instead of .raw files
# ============================================================================

from app.services.data_provider import DbDataProvider


def process_sample_from_db(
    *,
    src_mdb_conn,
    src_ch_client,
    sample_id: int,
    sample_name: str,
    cqcl_id: int,
    cqcl_reinj_id: int,
    cqn_id: int,
    batch_name: str,
    method_id: int,
    task_id: int,
    mass_tolerance_ppm: Optional[float] = None,
    in_competition: bool = True,
    is_beta_blocker: bool = True,
    validate_internal_standards: bool = True,
    year: Optional[int] = None,
    analysis_mdb=None,
    analysis_ch_client=None,
):
    """Entry point for processing ONE sample. Called in a loop by the
    orchestrator (once per regular sample in the batch).

    It is the direct parallel of the V1 project's `process_raw_file`, the only
    difference being where the chromatograms come from: here DbDataProvider
    reads the database; there MSFileReader opens the .raw file.

    PARAMETERS:
        src_mdb_conn / src_ch_client: connections to the source (the lake).
            Reused across every sample in the batch (opened by the caller). Do
            NOT confuse them with the analysis connection - these are for
            READING spectral points only.
        sample_id: id in the lake's `samples` table (the target sample).
        sample_name: readable name used in sample_results.sample_name and logs.
        cqcl_id, cqcl_reinj_id, cqn_id: ids of the batch's 3 controls in
            `samples`. The same for every sample in that batch.
        batch_name: the batch name (e.g. '20250828_0512_25'), used for the log
            and for the `pseudo_path` that ends up as sample_path in the logs.
        method_id: which method to use (configured in `methods`).
        task_id: id of the `sample_processing_tasks` row this run's
            sample_results belong to.
        mass_tolerance_ppm: see the DbDataProvider docstring. None for TVII.
        in_competition: comes from the lake's samples.competition. False -> skip
            substances and groups with competition=1 in the method. It also
            controls intra-group aggregation (ANY when True, ALL when False).
        is_beta_blocker: comes from the lake's samples.beta_block. False -> skip
            substances and groups with beta_blocker=1 in the method. Orthogonal
            to in_competition - the 2 filters add up.
        validate_internal_standards: True (the default) evaluates DTW/R2 for
            every internal standard in the method and aborts the sample if the
            IS gate fails (STEP B). Use False on validation/calibration batches
            where the check does not matter. The internal standard is still used
            normally by the cut-off in STEP D - only the initial validation is
            skipped.

    INTERNAL FLOW (summary - details in the module docstring):
        a) Creates 4 DbDataProviders (sample + 3 controls).
        b) Evaluates the internal standards (the IS gate blocks on failure).
        c) Reads the in_competition and is_beta_blocker flags of the sample.
        d) Loops over the individual substances -> process_single_substance
           (filtering competition and beta_blocker according to the flags).
        e) Loops over the groups -> process_single_substance + aggregates
           result_group (same filters, plus skipping competition-only/BB-only
           groups when the sample is not in the matching category).

    Returns a STATUS string describing the sample's outcome (stored by the
    caller in the task's `result` JSON - everything used to become "processed",
    including a sample the IS gate aborted silently, leaving a 'done' task with
    0 results and a result field that lied):
      - "processed"           -> full pipeline, substances analysed.
      - "blocked_by_is_gate"  -> the IS gate failed or had no data; NOTHING was
                                 analysed or written to sample_results.
      - "method_not_found"    -> the method does not exist in the analysis database.
    The side effects (INSERTs into SQLite + ClickHouse) are still the main
    product - the status is audit metadata.
    """
    # Pseudo-path: the V1 pipeline passed the real .raw path here for the logs.
    # Since there is no file, we generate a "db://" identifier just to keep it
    # traceable. It is NOT used to open anything - it is only a label.
    pseudo_path = f"db://{batch_name}/Data/{sample_name}.raw"

    # A DEDICATED connection to the ANALYSIS database - where
    # method_substances and cutoff_processing_groups live, and where
    # sample_results is written. Do NOT confuse it with src_mdb_conn (the lake -
    # only used to resolve channel_id inside the DbDataProviders).
    #
    # CONNECTION REUSE (see the WinError 10048/10055 issue):
    #   - `analysis_mdb` passed by the caller: reused across ALL the samples in
    #     the batch. The caller opens it once, reuses it for every sample and
    #     closes it at the end.
    #   - None: created and closed locally (legacy - ad-hoc use and scripts).
    own_conn = analysis_mdb is None
    conn = analysis_mdb if analysis_mdb is not None else maria_conn()
    cursor = conn.cursor(dictionary=True)

    # One log per BATCH (every sample of the same batch writes to one file).
    log_dir = os.path.join(_REPO_ROOT, "log")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{batch_name}.log")

    try:
        log_message(f"🧠 Starting sample: {sample_name} (batch {batch_name})", log_file)

        # ====================================================================
        # STEP A - instantiate the 4 "virtual MSFileReaders"
        # ====================================================================
        # Each DbDataProvider represents ONE data channel (one row in
        # `samples`). The first GetChroData call resolves its channel_id.
        # The 3 controls (cqcl, cqcl_reinj, cqn) are the SAME for every sample
        # in the batch - they arrive pre-resolved from the caller.
        raw_sample = DbDataProvider(src_mdb_conn, src_ch_client, sample_id,
                                    sample_label=sample_name,
                                    mass_tolerance_ppm=mass_tolerance_ppm,
                                    year=year)
        raw_cqcl = DbDataProvider(src_mdb_conn, src_ch_client, cqcl_id,
                                  sample_label="CQCL",
                                  mass_tolerance_ppm=mass_tolerance_ppm,
                                  year=year)
        raw_cqcl_reinj = DbDataProvider(src_mdb_conn, src_ch_client, cqcl_reinj_id,
                                        sample_label="CQCL_reinj",
                                        mass_tolerance_ppm=mass_tolerance_ppm,
                                        year=year)
        raw_cqn = DbDataProvider(src_mdb_conn, src_ch_client, cqn_id,
                                 sample_label="CQN",
                                 mass_tolerance_ppm=mass_tolerance_ppm,
                                 year=year)

        # ====================================================================
        # STEP B - fetch the method and identify the IS gate
        # ====================================================================
        cursor.execute("SELECT * FROM methods WHERE id = %s", (method_id,))
        method = cursor.fetchone()
        if not method:
            log_message(f"❌ Method {method_id} not found.", log_file)
            return "method_not_found"

        # The IS gate: the "primary" internal standard whose chromatogram MUST
        # be sound before the rest of the processing is allowed. If it fails we
        # block everything (most likely an injection or extraction failure).
        # Configured in cutoff_processing_groups.internal_standard.
        # The JOIN on method_substances plus deleted_at IS NULL: if the IS gate
        # was soft-deleted, it stops being the gate (it neither validates nor
        # blocks).
        cursor.execute("""
            SELECT cpg.internal_standard
            FROM cutoff_processing_groups cpg
            INNER JOIN method_substances ms ON ms.id = cpg.internal_standard
            WHERE cpg.method_id = %s
              AND ms.deleted_at IS NULL
            LIMIT 1
        """, (method_id,))
        row_gate = cursor.fetchone()
        gate_internal_standard_id = row_gate["internal_standard"] if row_gate else None

        # List every internal standard in the method. The gate is just one of
        # them; the others are informational - failures become alerts in
        # internal_standards_audit but do not block the analysis.
        cursor.execute("""
            SELECT method_substances.*, sf.name AS scan_filter
            FROM method_substances
            LEFT JOIN scan_filters sf ON sf.id = method_substances.scan_filter_id
            WHERE method_substances.method_id = %s
              AND method_substances.type = 'internal_standard'
              AND method_substances.deleted_at IS NULL
        """, (method_id,))
        internal_standards = [_attach_mass_range(r) for r in cursor.fetchall()]

        # Accumulator for the internal standards that failed (non-gate). Saved
        # into sample_processing_tasks.internal_standards_audit at the end.
        internal_standards_alert = {"Error": "Internal standard not identified", "samples": []}
        # Global flag: if it turns False, we abort processing for this sample.
        substances_present = True

        # ----- Internal standard evaluation -----
        # The same DTW + R2 logic used for the substances later, except the
        # window is centred on the peak of the internal standard's own CQCL.
        #
        # Skipped when validate_internal_standards=False. Useful for
        # validation/calibration batches - the internal standard is still used by
        # STEP D (cut-off) through cutoff_processing_groups; only the initial
        # validation and the IS gate are skipped.
        if not validate_internal_standards:
            log_message(
                "⏭️  Internal standard validation SKIPPED. The IS gate is not "
                "applied; the cut-off still uses the internal standard normally.",
                log_file,
            )
        for istd in (internal_standards if validate_internal_standards else []):
            log_message(f"\n🔬 Evaluating internal standard: {istd['name']}", log_file)
            dtw_limit = float(istd["dtw_limit"]) if istd.get("dtw_limit") is not None else 1.0
            r2_threshold = float(istd["r2_threshold"]) if istd.get("r2_threshold") is not None else -12.0

            data_sample = get_data_and_process(raw_sample, istd["start_time"], istd["end_time"],
                                               istd["mass_range"], istd["scan_filter"],
                                               istd["smoothing_type"], istd["smoothing_value"])
            data_cqcl = get_data_and_process(raw_cqcl, istd["start_time"], istd["end_time"],
                                             istd["mass_range"], istd["scan_filter"],
                                             istd["smoothing_type"], istd["smoothing_value"])

            if data_cqcl.empty or data_sample.empty:
                log_message(f"⚠️ Internal standard '{istd['name']}' has no valid chromatographic data", log_file)
                if gate_internal_standard_id is not None and istd["id"] == gate_internal_standard_id:
                    log_message("🚫 BLOCKED: the IS gate has no data.", log_file)
                    substances_present = False
                    break
                continue

            peak_idx = data_cqcl["intensity"].idxmax()
            rt_peak = data_cqcl.loc[peak_idx, "retention_time"]
            rt_window = 0.1
            sample_window = data_sample[(data_sample["retention_time"] >= rt_peak - rt_window) & (data_sample["retention_time"] <= rt_peak + rt_window)]
            cqcl_window = data_cqcl[(data_cqcl["retention_time"] >= rt_peak - rt_window) & (data_cqcl["retention_time"] <= rt_peak + rt_window)]

            if len(sample_window) < 3 or len(cqcl_window) < 3:
                log_message(f"⚠️ Internal standard '{istd['name']}' has too few points for DTW", log_file)
                if gate_internal_standard_id is not None and istd["id"] == gate_internal_standard_id:
                    log_message("🚫 BLOCKED: the IS gate has too little data.", log_file)
                    substances_present = False
                    break
                continue

            s_y = sample_window["intensity"].values
            c_y = cqcl_window["intensity"].values
            s_y = (s_y - np.min(s_y)) / (np.max(s_y) - np.min(s_y) + 1e-9)
            c_y = (c_y - np.min(c_y)) / (np.max(c_y) - np.min(c_y) + 1e-9)
            dist = dtw.distance(s_y, c_y)

            r_squared_float, r_squared_str = calculate_r_squared(sample_window, cqcl_window)
            r2_failed = r_squared_float is None or r_squared_float < r2_threshold
            is_gate = gate_internal_standard_id is not None and istd["id"] == gate_internal_standard_id

            if dist > dtw_limit or r2_failed:
                log_message(f"❌ IS '{istd['name']}' failed: DTW={dist:.4f} (limit {dtw_limit}), R2={r_squared_str} (min {r2_threshold})", log_file)
                if is_gate:
                    log_message("🚫 BLOCKED: the IS gate failed.", log_file)
                    substances_present = False
                    break
                else:
                    internal_standards_alert["samples"].append({"sample": sample_name, "internal_standard": istd["name"]})
            else:
                log_message(f"✅ IS '{istd['name']}' passed: DTW={dist:.4f} ≤ {dtw_limit}, R2={r_squared_str} ≥ {r2_threshold}", log_file)

        if internal_standards_alert["samples"]:
            log_message(f"⚠️ ALERT: internal standards with high variation -> {internal_standards_alert}", log_file)
            try:
                cur_audit = conn.cursor()
                cur_audit.execute("""
                    UPDATE sample_processing_tasks
                    SET internal_standards_audit = %s, updated_at = NOW()
                    WHERE id = %s
                """, (json.dumps(internal_standards_alert), task_id))
                conn.commit()
                cur_audit.close()
            except Exception as e:
                log_message(f"⚠️ Failed to save internal_standards_audit for task {task_id}: {e}", log_file)

        # ====================================================================
        # STEP C - the sample's flags (competition + beta_blocker)
        # ====================================================================
        # Both come from the lake's samples table (passed in as parameters); the
        # ingestion records them. NULL -> True in both cases (safe default - when
        # a sample has no flag set, we analyse EVERYTHING).
        #
        # Effect on the scan:
        #   competition=False  -> skip substances and groups with competition=1
        #   beta_blocker=False -> skip substances and groups with beta_blocker=1
        #
        # Group aggregation (ANY vs ALL members SUSPECT -> group SUSPECT) is
        # still governed ONLY by the competition flag (legacy behaviour).
        log_message(
            f"ℹ️ Sample '{sample_name}' "
            f"competition={'T' if in_competition else 'F'} "
            f"beta_blocker={'T' if is_beta_blocker else 'F'}",
            log_file,
        )

        # The IS gate blocked - we abort before spending time on the
        # substances. The "blocked_by_is_gate" status goes into the task's
        # `result` JSON; the caller used to mark it "processed" and the task
        # ended 'done' with 0 results and no trace of why (a silent failure).
        if not substances_present:
            log_message(
                "🚫 Sample blocked by the IS gate - nothing analysed "
                "(status: blocked_by_is_gate).",
                log_file,
            )
            return "blocked_by_is_gate"

        # ====================================================================
        # STEP D - the individual substances
        # ====================================================================
        # Loads EVERY substance in the method with type='substance'.
        # The JOIN with analysis_rules brings the `analysis_rule_code` that picks
        # which handler (standard/two_peaks/any_points) gets applied.
        # The individual-substance filter applies beta_blocker when the sample is
        # NOT in the BB category (skipping BB-only substances). The competition
        # filter is still handled in the else branch below (kept compatible with
        # the legacy design).
        bb_filter = "" if is_beta_blocker else " AND ms.beta_blocker = 0"

        cursor.execute(f"""
            SELECT ms.*, ar.code AS analysis_rule_code, sf.name AS scan_filter
            FROM method_substances ms
            LEFT JOIN analysis_rules ar ON ar.id = ms.analysis_rule_id
            LEFT JOIN scan_filters sf ON sf.id = ms.scan_filter_id
            WHERE ms.method_id = %s AND ms.type = 'substance'
              AND ms.deleted_at IS NULL{bb_filter}
        """, (method_id,))
        substances = [_attach_mass_range(r) for r in cursor.fetchall()]

        if in_competition:
            # In competition -> analyse EVERYTHING (already filtered by
            # beta_blocker above).
            for substance in substances:
                _, _, _, errors = process_single_substance(
                    cursor, task_id, pseudo_path,
                    raw_sample, raw_cqcl, raw_cqcl_reinj, raw_cqn, substance, log_file,
                    ch_client=analysis_ch_client,
                )
                if errors:
                    log_message(f"❌ Errors on {substance['name']}: {errors}", log_file)

            # ----------------------------------------------------------------
            # STEP E - substance groups (in_competition=True)
            # ----------------------------------------------------------------
            # Each group is a set of fragments of the SAME molecule. Every
            # member is processed individually (producing sample_results) and
            # then an update consolidates `result_group`:
            #   in_competition=True -> group SUSPECT if ANY member is SUSPECT
            #     (the more conservative rule: if one fragment matches, suspect).
            cursor.execute("SELECT * FROM substance_groups WHERE method_id = %s", (method_id,))
            for group in cursor.fetchall():
                # deleted_at IS NULL: soft-deleted members leave the group; the
                # rest keep being processed. A group with no active member falls
                # into the `if not members: continue` below.
                cursor.execute("""
                    SELECT method_substances.*, ar.code AS analysis_rule_code, sf.name AS scan_filter
                    FROM substance_group_memberships
                    INNER JOIN method_substances ON substance_group_memberships.method_substance_id = method_substances.id
                    LEFT JOIN analysis_rules ar ON ar.id = method_substances.analysis_rule_id
                    LEFT JOIN scan_filters sf ON sf.id = method_substances.scan_filter_id
                    WHERE substance_group_memberships.substance_group_id = %s
                      AND method_substances.deleted_at IS NULL
                """, (group["id"],))
                members = [_attach_mass_range(r) for r in cursor.fetchall()]
                if not members:
                    continue
                # Skip beta-blocker groups when the sample is not BB. The flag
                # is read off the first member (every member of a group is
                # expected to carry the same value - the configuration
                # guarantees it).
                if not is_beta_blocker and members[0].get("beta_blocker"):
                    continue
                group_results, group_ids = [], []
                for sub in members:
                    sample_result_id, result, _, _ = process_single_substance(
                        cursor, task_id, pseudo_path,
                        raw_sample, raw_cqcl, raw_cqcl_reinj, raw_cqn, sub, log_file,
                        ch_client=analysis_ch_client,
                    )
                    group_results.append(result)
                    if sample_result_id:
                        group_ids.append(sample_result_id)
                # GROUP AGGREGATION RULE (applies to both branches):
                # the group is SUSPECT only when ALL members are SUSPECT.
                # Any other case (NEGATIVE, INVALID, a mix) -> NEGATIVE.
                # INVALID is never stored in sample_results
                # (process_single_substance does not insert on INVALID), so
                # `all()` naturally returns False and the group becomes NEGATIVE
                # - without marking a row that does not exist.
                final_result = "SUSPECT" if all(r == "SUSPECT" for r in group_results) else "NEGATIVE"
                for sid in group_ids:
                    cursor.execute("UPDATE sample_results SET result_group = %s WHERE id = %s",
                                   (final_result, sid))
                conn.commit()
        else:
            # Out of competition -> only substances with competition=0
            # (substances banned in competition only drop off the list).
            # The beta_blocker filter applies here too, with the same semantics.
            cursor.execute(f"""
                SELECT ms.*, ar.code AS analysis_rule_code, sf.name AS scan_filter
                FROM method_substances ms
                LEFT JOIN analysis_rules ar ON ar.id = ms.analysis_rule_id
                LEFT JOIN scan_filters sf ON sf.id = ms.scan_filter_id
                WHERE ms.method_id = %s AND ms.type = 'substance' AND ms.competition = 0
                  AND ms.deleted_at IS NULL{bb_filter}
            """, (method_id,))
            for sub in (_attach_mass_range(r) for r in cursor.fetchall()):
                process_single_substance(cursor, task_id, pseudo_path,
                                         raw_sample, raw_cqcl, raw_cqcl_reinj, raw_cqn, sub, log_file,
                                         ch_client=analysis_ch_client)

            # Groups: the inverse rule of the in_competition=True branch.
            #   in_competition=False -> the group is SUSPECT only if ALL members
            #     are (stricter - out of competition it demands convergence).
            cursor.execute("SELECT * FROM substance_groups WHERE method_id = %s", (method_id,))
            for group in cursor.fetchall():
                # deleted_at IS NULL: soft-deleted members drop out; the group
                # carries on with the active ones (or disappears if none are
                # left).
                cursor.execute("""
                    SELECT method_substances.*, ar.code AS analysis_rule_code, sf.name AS scan_filter
                    FROM substance_group_memberships
                    INNER JOIN method_substances ON substance_group_memberships.method_substance_id = method_substances.id
                    LEFT JOIN analysis_rules ar ON ar.id = method_substances.analysis_rule_id
                    LEFT JOIN scan_filters sf ON sf.id = method_substances.scan_filter_id
                    WHERE substance_group_memberships.substance_group_id = %s
                      AND method_substances.deleted_at IS NULL
                """, (group["id"],))
                members = [_attach_mass_range(r) for r in cursor.fetchall()]
                if not members:
                    continue
                first_sub = members[0]
                # Skip competition-only (legacy) or BB-only groups when the
                # sample is not in that category. The flag is assumed uniform
                # across the group's members.
                if first_sub["competition"]:
                    continue
                if not is_beta_blocker and first_sub.get("beta_blocker"):
                    continue
                group_results, group_ids = [], []
                for sub in members:
                    sample_result_id, result, _, _ = process_single_substance(
                        cursor, task_id, pseudo_path,
                        raw_sample, raw_cqcl, raw_cqcl_reinj, raw_cqn, sub, log_file,
                        ch_client=analysis_ch_client,
                    )
                    group_results.append(result)
                    if sample_result_id:
                        group_ids.append(sample_result_id)
                # The same rule as the in_competition=True branch (see the
                # comment there). SUSPECT only when ALL members are SUSPECT;
                # otherwise NEGATIVE.
                final_result = "SUSPECT" if all(r == "SUSPECT" for r in group_results) else "NEGATIVE"
                for sid in group_ids:
                    cursor.execute("UPDATE sample_results SET result_group = %s WHERE id = %s",
                                   (final_result, sid))
                conn.commit()

        return "processed"

    except Exception as e:
        log_message(f"❌ General error while processing {sample_name}: {e}", log_file)
        traceback.print_exc()
        raise
    finally:
        # The cursor is ALWAYS ours - always close it.
        try:
            cursor.close()
        except Exception:
            pass
        # The analysis connection: only closed when WE created it locally
        # (legacy mode). When it came from the caller (the hot path), the caller
        # owns the lifecycle.
        if own_conn:
            try:
                conn.close()
            except Exception:
                pass

