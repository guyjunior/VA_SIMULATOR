"""
virtual_analyst.py - the orchestration of the processing.

THE INTERPRETATION STEPS FOLLOW THE ARTICLE. The decision for each substance
(Step 1 centroid protocols, Step 2 Mamdani inference over the Pseudo R2 with
centroid defuzzification, Step 3 DTW normalised by (n+m) as a rescue, Step 4
C_hat / cut-off bands for MRL substances) lives in
`app/services/article_pipeline.py`, which implements Sections 2.3 and 2.5-2.8
of the article, Table 1 and Eq. 4. This module extracts the chromatograms,
subtracts the NQC, resamples them onto the 3000-point grid, hands them to
article_pipeline and stores the verdicts.

Kept from the application this was ported from, because the article's
interpretation steps do not cover them: the internal-standard gate on the
injection (step 2 below), the fallback to the reinjected PQC when the PQC gives
NEGATIVE, the substance groups and the persistence.

The plumbing: the analysis database is SQLite (app/db/sqlite.py accepts the
`%s` dialect) and ClickHouse comes from this project's docker-compose. The
chromatograms are read through DbDataProvider, which answers `GetChroData`
like MSFileReader (same smoothing kernel) but reads the lake database instead
of opening .raw files.

PQC = CQCL (positive quality control), NQC = CQN (negative quality control).

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
       pass or fail one without affecting the other.

    4. PROCESS EACH INDIVIDUAL SUBSTANCE
       Loop over method_substances.type='substance' (filtered by competition and
       beta_blocker according to the sample's flags). For each one it calls
       process_single_substance, which:
         a. Extracts the XIC of sample, cqcl, cqcl_reinj, cqcn (4 chromatograms)
         b. Applies adjust_intensity (subtracts the CQN baseline)
         c. Resamples the corrected traces onto the 3000-point grid and keeps
            the substance's RT window (0.8 min in the screening method)
         d. Runs the article's Steps 1-4 (article_pipeline.evaluate_peak);
            analysis_rule_code selects the variant:
              "standard"  -> Steps 1-4
              "any_points"-> Steps 1-4 with Protocol #2 (substances with no MRL)
              "two_peaks" -> Steps 1-4 on each PQC peak, most severe kept
            If the PQC gives NEGATIVE, the reinjected PQC is tried as well.
         e. Substances with a cut_off get their level from Step 4 (C_hat
            against the internal standard closest in tR)
         f. INSERTs into sample_results + stores the XICs in the analysis
            ClickHouse

    5. PROCESS THE SUBSTANCE GROUPS (substance_groups)
       Grouped substances are fragments of the same molecule. Each member
       becomes a row in sample_results, and then an update applies the
       consolidated `result_group`:
         the group is SUSPECT only when ALL its members are suspect (at any
         level); otherwise NEGATIVE - see group_result(). Same rule in and out
         of competition.

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

    Pseudo R2 = r2_score(y_true=sample, y_pred=PQC) on max-normalised traces
          inside the RT window. The ONLY input of the Mamdani system.
    DTW = Dynamic Time Warping distance divided by (n + m) (Eq. 4). Used only
          when the Pseudo R2 fails (a rescue), with bands 0.002/0.008/0.120.

    Output levels: NEGATIVE, SUSPECT_VERY_LOW, SUSPECT_LOW, SUSPECT_MED,
    SUSPECT_HIGH, SUSPECT_VERY_HIGH.
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
from app.services import article_pipeline as ap

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

# Severity order of the verdicts (shared with the article module).
_SEVERITY_RANK = ap.SEVERITY


def is_suspect_label(label: str) -> bool:
    return isinstance(label, str) and "SUSPECT" in label.upper()


def group_result(member_results) -> str:
    """Consolidated `result_group`: SUSPECT only when EVERY member is suspect
    at any level (SUSPECT_VERY_LOW ... SUSPECT_VERY_HIGH). A member that is
    NEGATIVE or INVALID - or a group with no members - gives NEGATIVE."""
    results = list(member_results)
    return "SUSPECT" if results and all(is_suspect_label(r) for r in results) else "NEGATIVE"

def calculate_r_squared(sample, cqcl):
    """R2 used ONLY by the internal-standard gate (injection QC, outside the
    interpretation steps of the article). The substances use
    article_pipeline.pseudo_r2."""
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


def nearest_internal_standard(cursor, substance):
    """Step 4: the internal standard used for C_hat is the one CLOSEST IN tR to
    the substance (distance between the centres of the RT windows)."""
    cursor.execute("""
        SELECT ms.*, sf.name AS scan_filter
        FROM method_substances ms
        LEFT JOIN scan_filters sf ON sf.id = ms.scan_filter_id
        WHERE ms.method_id = %s
          AND ms.type = 'internal_standard'
          AND ms.deleted_at IS NULL
    """, (substance["method_id"],))
    rows = [_attach_mass_range(r) for r in cursor.fetchall()]
    if not rows:
        return None
    centre = (float(substance["start_time"]) + float(substance["end_time"])) / 2.0
    return min(rows, key=lambda r: abs((float(r["start_time"]) + float(r["end_time"])) / 2.0 - centre))


def make_cutoff_estimator(ctx, raw_reference, data_reference):
    """Build the Step 4 callable for substances with an MRL (cut_off set);
    None for the others. It returns (C_hat, cut_off).

    C_hat = (A_sample / A_IS,sample) / (A_PQC / A_IS,PQC) * C_PQC, with the
    peak areas from cut_off_verification and C_PQC = fortified_concentration.
    """
    substance = ctx["substance"]
    if substance.get("cut_off") is None:
        return None

    def estimate():
        cut_off = float(substance["cut_off"])
        c_pqc = substance.get("fortified_concentration")
        if not c_pqc:
            log_message(f"⚠️ [{ctx['s_name']}] MRL substance without fortified_concentration - C_hat not estimable", ctx["log_file"])
            return None, cut_off
        istd = nearest_internal_standard(ctx["cursor"], substance)
        if istd is None:
            log_message(f"⚠️ [{ctx['s_name']}] no internal standard in the method - C_hat not estimable", ctx["log_file"])
            return None, cut_off
        is_sample = get_data_and_process(ctx["raw_sample"], istd["start_time"], istd["end_time"],
                                         istd["mass_range"], istd["scan_filter"],
                                         istd["smoothing_type"], istd["smoothing_value"])
        is_reference = get_data_and_process(raw_reference, istd["start_time"], istd["end_time"],
                                            istd["mass_range"], istd["scan_filter"],
                                            istd["smoothing_type"], istd["smoothing_value"])
        c_hat, zero_area, _areas = cut_off_verification(
            ctx["data_sample"], data_reference, is_sample, is_reference, float(c_pqc)
        )
        if c_hat is None:
            log_message(f"⚠️ [{ctx['s_name']}] C_hat: zero area in '{zero_area}' (IS '{istd['name']}')", ctx["log_file"])
        return c_hat, cut_off

    return estimate


def evaluate_against_reference(ctx, reference: str):
    """Run Steps 1-4 of the article (article_pipeline) for one substance
    against one reference injection: reference = "cqcl" (the PQC) or
    "cqcl_reinj". Returns (verdict, concentration)."""
    sub = ctx["substance"]
    s_name = ctx["s_name"]
    sample = ctx["grid_sample"]
    pqc = ctx[f"grid_{reference}"]
    raw_reference = ctx["raw_cqcl"] if reference == "cqcl" else ctx["raw_cqcl_reinj"]
    data_reference = ctx["data_cqcl"] if reference == "cqcl" else ctx["data_cqcl_reinj"]

    r2_threshold = sub.get("r2_threshold")
    ratio_limit = sub.get("any_points_ratio_limit")
    kwargs = dict(
        r2_threshold=float(r2_threshold) if r2_threshold is not None else None,
        any_points=(ctx["rule_code"] == "any_points"),
        ratio_limit=float(ratio_limit) if ratio_limit is not None else ap.PROTOCOL2_RATIO_LIMIT,
        cutoff_estimator=make_cutoff_estimator(ctx, raw_reference, data_reference),
    )

    if ctx["rule_code"] == "two_peaks":
        segments = ap.split_two_peaks(ctx["grid_rt"], pqc)
        if segments:
            decisions = [ap.evaluate_peak(sample[lo:hi], pqc[lo:hi], **kwargs) for lo, hi in segments]
            for k, d in enumerate(decisions, 1):
                log_message(f"ℹ️ [{s_name}] vs {reference} two_peak#{k}: {d.describe()}", ctx["log_file"])
            best = max(decisions, key=lambda d: ap.SEVERITY.get(d.verdict, -2))
            return best.verdict, best.concentration

    d = ap.evaluate_peak(sample, pqc, **kwargs)
    log_message(f"ℹ️ [{s_name}] vs {reference}: {d.describe()}", ctx["log_file"])
    return d.verdict, d.concentration

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

def run_article_rules(ctx):
    """The decision for one substance: the article's Steps 1-4 against the PQC
    (CQCL). When that is NEGATIVE, the reinjected PQC is tried as well - the
    reinjection fallback is not part of the article's interpretation steps and
    is kept exactly as in the application it was ported from."""
    result, concentration = evaluate_against_reference(ctx, "cqcl")
    if result == "INVALID":
        print(f"ℹ️ [{ctx['s_name']}] INVALID on the CQCL - skipping.")
        return {"invalid": True, "errors": []}

    if result == "NEGATIVE":
        result2, concentration2 = evaluate_against_reference(ctx, "cqcl_reinj")
        if result2 == "INVALID":
            print(f"ℹ️ [{ctx['s_name']}] INVALID on the REINJ - skipping.")
            return {"invalid": True, "errors": []}
        if result2 != "NEGATIVE":
            result, concentration = result2, concentration2

    return {"result": result, "concentration": concentration, "errors": []}

# Rule codes understood by run_article_rules (method_substances.analysis_rule):
#   "standard"   -> Steps 1-4
#   "any_points" -> Steps 1-4 with Protocol #2 enabled (substances with no MRL)
#   "two_peaks"  -> Steps 1-4 on each of the two PQC peaks (co-eluting isomers)
ARTICLE_RULE_CODES = ("standard", "any_points", "two_peaks")

# The "do not interpret" rule (VISUALISATION-only substances): the XICs are
# stored, nothing is judged.
NO_ANALYZE_RULE_CODE = "no_analyze"

def process_single_substance(cursor, task_id, sample_path,raw_sample, raw_cqcl, raw_cqcl_reinj, raw_cqcn, substance, log_file, ch_client=None):
    """Process ONE substance for ONE sample. Called in a loop for every
    substance in the method.

    Steps:
      1. Extract the 4 raw chromatograms (sample, cqcl, cqcl_reinj, cqn) over
         the substance's window (RT range, m/z range, smoothing).
      2. Apply adjust_intensity (X - CQN) to the first 3 - that gives the
         "_corrected" ones.
      3. Check that the corrected CQCL has points (otherwise INVALID, early exit).
      4. Resample the corrected traces onto the 3000-point grid and keep the
         substance's RT window.
      5. Run the article's Steps 1-4 (run_article_rules -> article_pipeline).
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
        result: 'NEGATIVE' | 'SUSPECT_VERY_LOW/LOW/MED/HIGH/VERY_HIGH'
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
            # analyst too. That is why the branch comes before the early
            # returns of the interpretation path.
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

            # === Section 2.3: the 3000-point grid + the substance's RT window ===
            grid = ap.make_grid()
            mask = ap.window_mask(grid, substance["start_time"], substance["end_time"])
            if mask.sum() < 2:
                print(f"⚠️ [{s_name}] RT window holds fewer than 2 grid points - ignoring.")
                return None, "INVALID", "0.0", []

            def on_grid(df):
                if df is None or df.empty:
                    return np.zeros(int(mask.sum()))
                return ap.to_grid(df["retention_time"], df["intensity"], grid)[mask]

            if rule_code not in ARTICLE_RULE_CODES:
                print(f"⚠️ [{s_name}] analysis_rule '{rule_code}' unknown - using 'standard'")
                rule_code = "standard"

            ctx = {
                "cursor": cursor,
                "substance": substance,
                "s_name": s_name,
                "log_file": log_file,
                "rule_code": rule_code,
                "raw_sample": raw_sample,
                "raw_cqcl": raw_cqcl,
                "raw_cqcl_reinj": raw_cqcl_reinj,
                "data_sample": data_sample,
                "data_cqcl": data_cqcl,
                "data_cqcl_reinj": data_cqcl_reinj,
                "grid_rt": grid[mask],
                "grid_sample": on_grid(corrected_sample),
                "grid_cqcl": on_grid(corrected_cqcl),
                "grid_cqcl_reinj": on_grid(corrected_reinj),
            }

            handler_result = run_article_rules(ctx)

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
            substances and groups with competition=1 in the method.
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
            # then an update consolidates `result_group` (group_result).
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
                # the group is SUSPECT only when ALL members are suspect at any
                # level. Any other case (NEGATIVE, INVALID, a mix) -> NEGATIVE.
                # INVALID is never stored in sample_results
                # (process_single_substance does not insert on INVALID), so no
                # row is marked for it.
                final_result = group_result(group_results)
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

            # Groups: the same rule as the in_competition=True branch
            # (group_result); competition-only groups are skipped.
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
                # The same rule as the in_competition=True branch.
                final_result = group_result(group_results)
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

