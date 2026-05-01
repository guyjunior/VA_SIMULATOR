# Method Simulator — Virtual Analyst

> An interactive, single-file HTML application for simulating, validating and benchmarking a chromatographic decision pipeline used in **anti-doping analysis** (WADA prohibited-list substances). The Virtual Analyst combines **R² coefficient evaluation under fuzzy logic**, **Dynamic Time Warping (DTW)** shape-similarity, and **ITP-Cutoff threshold rules** to produce per-substance verdicts on whether a chromatographic peak should be flagged as suspect.

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status: Active](https://img.shields.io/badge/Status-Active-brightgreen.svg)
![Tech: Vanilla JS](https://img.shields.io/badge/Tech-Vanilla%20JS-yellow.svg)
![Domain: Anti-Doping](https://img.shields.io/badge/Domain-Anti--Doping-blue.svg)

---

## Table of Contents

1. [Abstract](#abstract)
2. [Motivation and Scope](#motivation-and-scope)
3. [Application Architecture](#application-architecture)
4. [Theoretical Background](#theoretical-background)
   - [R² Coefficient of Determination](#r-coefficient-of-determination)
   - [Dynamic Time Warping (DTW)](#dynamic-time-warping-dtw)
   - [Fuzzy Logic Classification](#fuzzy-logic-classification)
   - [ITP-Cutoff Concentration Rule](#itp-cutoff-concentration-rule)
5. [Pipeline Decision (4-Step Flow)](#pipeline-decision-4-step-flow)
6. [Method Configuration](#method-configuration)
7. [Worked Examples (Scenarios)](#worked-examples-scenarios)
8. [Interactive Analyst Controls](#interactive-analyst-controls)
9. [Experiment Metrics & Validation Results](#experiment-metrics--validation-results)
10. [User Interface Tabs](#user-interface-tabs)
11. [Tech Stack](#tech-stack)
12. [Getting Started](#getting-started)
13. [File Structure](#file-structure)
14. [Author](#author)
15. [License](#license)

---

## Abstract

The **Method Simulator — Virtual Analyst** is a self-contained, browser-based simulator that reproduces, in didactic and reproducible form, the decision logic of an automated chromatographic peak classifier deployed in routine anti-doping screening. The system addresses a recurring problem in liquid- and gas-chromatography mass-spectrometry (LC-MS / GC-MS) workflows: how to combine quantitative goodness-of-fit metrics (R²), shape-similarity measures (DTW), and concentration-based thresholds (ITP-Cutoff) into a coherent, interpretable, and tunable decision rule that minimises false negatives without overwhelming analysts with false positives.

The application provides three complementary views: (i) a **METHOD** tab for substance-by-substance parameter configuration; (ii) an **ANALYST** tab with three pre-loaded pedagogical scenarios drawn from the WADA prohibited list (S9 glucocorticoid Budesonide; S1 anabolic agent Ostarine; S4 hormone-and-metabolic-modulator Clomifene/Toremifene metabolites), each illustrating a distinct chromatographic challenge; and (iii) an **EXPERIMENT METRICS** tab presenting the empirical validation of the pipeline against an expert-reviewed reference set of 2,004 samples × 274 substances (549,096 total interpretations).

---

## Motivation and Scope

Routine anti-doping laboratories process thousands of samples per year, each screened against hundreds of prohibited substances. Fully manual peak review is neither scalable nor consistent across analysts. Conversely, naïve threshold-based automation tends to either (a) miss subtle peaks ("false negatives" — the most consequential error in this domain) or (b) flag too many irrelevant signals, increasing reviewer workload and degrading trust in the system.

This simulator was conceived to:

- **Make the decision logic transparent.** Each step of the pipeline is exposed in the UI with its inputs, thresholds, and intermediate values, allowing analysts and method developers to audit *why* a verdict was issued.
- **Permit substance-specific tuning.** Different chemical classes have different signal-to-noise characteristics, retention-time stability, and peak-shape morphology. The simulator allows per-substance adjustment of every threshold.
- **Provide reproducible benchmarks.** The bundled validation experiment (Section 9) reports a complete confusion matrix, allowing the community to assess the trade-offs of the chosen logic.
- **Support training.** The three pre-built scenarios are deliberately chosen to expose edge-case behaviours that frequently confuse less experienced reviewers.

---

## Application Architecture

The application is delivered as a **single HTML file** (`VIRTUAL_ANALYST.html`) with no build step, no external dependencies beyond a CDN-served charting library, and no server component. All state is maintained in JavaScript memory; all rendering is performed via the HTML5 Canvas API (chromatogram traces) and Chart.js (metrics dashboards).

The high-level structure is:

```
┌──────────────────────────────────────────────────────────────────┐
│  Top Bar (logo, title, theme toggle)                             │
├──────────────────────────────────────────────────────────────────┤
│  Tab Bar:  [ METHOD ] [ EXAMPLE 1 ] [ EXAMPLE 2 ] [ EXAMPLE 3 ]  │
│                       [ EXPERIMENT METRICS ]                     │
├──────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Active tab content                                              │
│   - METHOD:    per-substance parameter form                      │
│   - EXAMPLE *: live chromatogram, fuzzy panel, pipeline,         │
│                verdict box, control sliders                      │
│   - METRICS:   KPI cards, confusion matrix, charts, table        │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘
```

A central `substances` array holds, for each compound, both its **METHOD parameters** (retention-time window, mass range, scan filter, smoothing settings, R² and DTW thresholds, ITP-cutoff, integration factor) and its **per-scenario simulation seed** (deterministic noise vector, peak descriptors). Switching scenarios re-renders the chromatogram canvas and re-evaluates the pipeline against the currently configured thresholds in real time.

---

## Theoretical Background

### R² Coefficient of Determination

The coefficient of determination is computed between the *observed* chromatographic trace (after smoothing, within the substance-specific retention window) and a *reference template* representing the expected peak shape. Formally,

$$R^{2} = 1 - \frac{\sum_{i}(y_{i} - \hat{y}_{i})^{2}}{\sum_{i}(y_{i} - \bar{y})^{2}}$$

where $y_{i}$ is the observed signal at sampling point $i$, $\hat{y}_{i}$ the reference template value, and $\bar{y}$ the mean of the observed signal. Values approach **1.0** for a perfect fit and may be strongly negative when the observed window contains noise, baseline drift, or interfering co-elutions. In this simulator, R² is the **primary input to the fuzzy classifier** (Step 2 of the pipeline) and the only metric to which fuzzy logic is applied; all other steps use crisp thresholds.

### Dynamic Time Warping (DTW)

DTW measures the *minimum cumulative cost* of aligning two time-series of potentially different lengths or non-uniform local stretching. It is used here as a **shape-similarity fallback** when R² fails: a peak that is shifted in retention time, slightly broadened, or asymmetrically skewed may still match the expected morphology under a warped alignment, even if its point-wise R² collapses. The DTW score is normalised to a comparable scale and accepted when it falls **below** a configurable threshold (default range 0.020–0.120 depending on substance), denoting "shape close enough to remain suspect."

### Fuzzy Logic Classification

Fuzzy logic is applied **exclusively to the R² value** (a deliberate design choice highlighted in the UI's "fuzzy scope badge"). Five trapezoidal/triangular membership functions partition the R² domain:

| MF Name    | Linguistic Meaning           | Approx. R² Region |
|------------|------------------------------|-------------------|
| VERY LOW   | No fit — likely noise/blank  | R² < −2           |
| LOW        | Poor fit — interference      | −2 ≤ R² < 0       |
| MEDIUM     | Borderline                   | 0 ≤ R² < 0.6      |
| HIGH       | Good fit                     | 0.6 ≤ R² < 0.9    |
| VERY HIGH  | Excellent fit — suspect peak | R² ≥ 0.9          |

The classifier computes membership degrees, fires a small rule base, and produces a **defuzzified centroid** that drives the colour-coded R² bar in the Analyst panel. The rule firing strengths are exposed as horizontal bars so the user can audit the inference.

### ITP-Cutoff Concentration Rule

Even when R² and DTW agree that a peak is morphologically suspect, the substance-specific **ITP (International Threshold for Procedure) cutoff** acts as a final concentration guard: peaks below the cutoff (after integration) are reported as sub-threshold and not escalated, in line with WADA Technical Documents that define minimum reporting levels for each prohibited substance.

---

## Pipeline Decision (4-Step Flow)

Every interpretation flows through four sequential, auditable steps. The UI shows each step as a card that lights up green/red/skipped according to the live evaluation:

1. **Step 1 — Window Cut.**
   The signal is truncated to a substance-specific retention-time window centred on the expected PQC (Peak Quality Control) retention time. This eliminates upstream and downstream noise that would otherwise contaminate the R² calculation.

2. **Step 2 — R² Evaluation under Fuzzy Logic.**
   The cropped signal is correlated against the reference template; the resulting R² is fuzzified, and the defuzzified output determines whether the peak is provisionally accepted into the **fuzzy_r2 suspect zone**. If R² ≥ the configured threshold, the pipeline can short-circuit to Step 4. Otherwise, Step 3 is invoked as a fallback.

3. **Step 3 — DTW Threshold Check (Fallback).**
   When R² fails Step 2, the DTW distance between the observed and reference traces is computed. If DTW ≤ the configured threshold, the peak is rescued into the **fuzzy_dtw suspect zone** — recognising that shape similarity may persist where point-wise correlation does not. If both R² and DTW fail, the peak is classified as **negative** and the pipeline terminates.

4. **Step 4 — ITP-Cutoff Decision.**
   For peaks accepted by Step 2 or Step 3, the integrated concentration is compared against the substance's ITP cutoff. Above-cutoff peaks are flagged **POSITIVE (suspect)**; below-cutoff peaks are flagged **sub-threshold negative**, with the reason recorded in the verdict box.

The verdict box at the bottom of the Analyst panel displays both the final classification and a human-readable reason string (e.g. *"Reason: fuzzy_dtw rescue, above ITP-cutoff"*), enabling traceability for downstream review.

---

## Method Configuration

The METHOD tab exposes, per substance, the following configurable parameters:

| Parameter           | Description                                                              |
|---------------------|--------------------------------------------------------------------------|
| **Name**            | Substance identifier (matching the WADA category code).                  |
| **Start time**      | Lower bound of the retention-time window (minutes).                      |
| **End time**        | Upper bound of the retention-time window (minutes).                      |
| **Mass range**      | Quadrupole/Orbitrap m/z window to be extracted.                          |
| **Scan filter**     | Acquisition channel / SRM transition selector.                           |
| **Type**            | Substance class (glucocorticoid, SARM, SERM metabolite, etc.).           |
| **Smoothing type**  | Savitzky–Golay, moving average, etc.                                     |
| **Smoothing value** | Window size or polynomial order for the smoother.                        |
| **Window cut**      | Width of the Step 1 crop around the PQC retention time.                  |
| **Concentration**   | Reference / spike concentration used to scale the simulated peak.        |
| **ITP-Cutoff**      | Minimum reportable concentration (Step 4 threshold).                     |
| **Has ITP-Cutoff?** | Toggle to disable Step 4 for substances without a defined MRPL.          |
| **R² threshold**    | Slider (default range −12.00 to 1.00) controlling Step 2 acceptance.     |
| **DTW threshold**   | Slider (default range 0.020–0.120) controlling Step 3 fallback.          |
| **Integration factor** | Calibration constant translating peak area to concentration.          |

Configuration changes propagate immediately to the active scenario, allowing real-time exploration of parameter sensitivity.

---

## Worked Examples (Scenarios)

The Analyst tab ships with three deliberately chosen examples, each illustrating a distinct chromatographic phenomenon:

### Example 1 — S9 Budesonide (Glucocorticoid)

A textbook well-resolved peak with high signal-to-noise. Illustrates the **happy path**: R² ≥ threshold, fuzzy classifier outputs HIGH/VERY HIGH, Step 3 is skipped, and Step 4 issues a clean POSITIVE verdict. Reference concentration: 30.00 (arb. units) with ITP-cutoff at 15.00.

### Example 2 — S1 Ostarine (SARM, "Lynx Eyes")

A low-amplitude, noise-embedded peak that requires careful smoothing and a permissive R² threshold to be recovered. Demonstrates the **fuzzy_r2 borderline regime** — small parameter changes flip the verdict, and the fuzzy panel shows competing membership activations. The "Lynx Eyes" nickname refers to the analyst skill traditionally required to spot such peaks manually.

### Example 3 — S4 Clomifene (4-OH)* / Toremifene (4-OH)* (Double Peak)

Co-eluting hydroxylated metabolites of two SERMs producing a **double-peak signature**. This scenario stresses the R² calculation (which compares against a single-peak template and therefore degrades) and showcases the **DTW fallback** (Step 3): the warped alignment recognises the morphology as still consistent with a suspect signature, rescuing what would otherwise be a false negative.

---

## Interactive Analyst Controls

The Analyst panel exposes a control strip with the following live sliders and toggles, all of which redraw the chromatogram and re-evaluate the pipeline on every change:

- **Noise** — additive white-noise amplitude.
- **Sigma (σ)** — Gaussian peak width.
- **Centroid shift** — retention-time offset of the peak relative to its expected position.
- **Amplitude** — peak height scaling.
- **Skew** — peak asymmetry (positive = tailing; negative = fronting).
- **Smoothing** — runtime smoother strength (independent of method-defined smoothing).
- **Crop toggle** — enable/disable the Step 1 window cut for visualisation.

Two Canvas-rendered charts show the **input signal** (raw + smoothed overlay) and the **output / reference comparison**, with the active retention-time window highlighted.

A dedicated **Fuzzy Panel** displays:

- The defuzzified centroid value.
- The five MF activation levels as horizontal bars.
- The current R² value with its linguistic interpretation (e.g. "Perfect fit", "Borderline", "Poor fit").
- A scope badge clarifying that fuzzy logic is **applied exclusively to R² (Step 2)**.

The **Pipeline Panel** renders the four steps in compact card form, each annotated with its status flag (`triggered`, `skipped`, `positive`, `negative`) and current value, enabling at-a-glance audit of the decision flow.

---

## Experiment Metrics & Validation Results

The EXPERIMENT METRICS tab presents the empirical validation of the Virtual Analyst against an expert-reviewed reference dataset.

### Dataset

| Quantity                     | Value             |
|------------------------------|-------------------|
| Samples                      | 2,004             |
| Substances on prohibited list| 274               |
| **Total interpretations**    | **549,096**       |
| Evaluated (with signal)      | 234,735           |
| Blank windows (no signal, TN)| 314,361           |

Each interpretation is a sample × substance pair; "blank windows" are pairs where the instrument detected no chromatographic signal at all and which are therefore unambiguous true negatives. Both viewpoints (with and without blank windows) are reported to allow the reader to assess inflation of TN-derived metrics.

### Confusion Matrix (full, including blank windows)

|                         | Predicted Positive | Predicted Negative |
|-------------------------|-------------------:|-------------------:|
| **Actually Positive**   | TP = 2,110         | FN = 6             |
| **Actually Negative**   | FP = 749           | TN = 546,231       |

### Key Performance Indicators

| Metric         | Value     | Interpretation                                                                                               |
|----------------|-----------|--------------------------------------------------------------------------------------------------------------|
| Accuracy       | **99.86%**| Of all 549,096 interpretations, only 755 disagreed with the expert verdict.                                  |
| Precision      | **73.80%**| Of all flagged positives, ~74% are confirmed suspect; remainder are reviewer-resolvable false alarms.        |
| Recall         | **99.72%**| Of 2,116 truly suspect substances, 2,110 were detected — only 6 missed (FN).                                 |
| F1 Score       | **84.82%**| Harmonic mean of precision and recall.                                                                       |
| Specificity    | **99.86%**| 546,231 of 546,980 true negatives correctly classified.                                                      |
| FPR            | **0.14%** | Down from 0.32% when blank windows are included — false alarms are proportionally rare across the full set.  |
| AUC (ROC)      | **0.9979**| Near-ideal operating characteristic; rises from 0.9970 to 0.9979 once blank windows are factored in.         |

### Visualisations

The metrics tab includes:

- **Overview pie** — evaluated vs blank-window split of the 549,096 interpretations.
- **Comparison bar** — KPIs side-by-side under the two TN-inclusion regimes.
- **Composition pie** — TP/FP/TN/FN proportions.
- **TN breakdown bar** — substance-evaluated vs blank-window contributions to TN.
- **Top 20 TP chart** — substances most frequently and correctly flagged.
- **Top 20 FP chart** — substances driving the precision deficit (candidates for threshold re-tuning).
- **FN chart** — the six undetected suspect cases, by substance.
- **Bubble chart (Precision × Recall × FPR)** — per-substance operating points.
- **Per-substance metrics table** — full confusion matrix and KPIs for each of the 274 substances, with TN expanded by approximately 1,147 blank windows per substance.

---

## User Interface Tabs

| Tab                  | Purpose                                                                      |
|----------------------|------------------------------------------------------------------------------|
| **METHOD**           | Edit the parameters (thresholds, windows, smoothing, cutoffs) per substance. |
| **EXAMPLE 1 (S9)**   | Run the Virtual Analyst on the Budesonide reference scenario.                |
| **EXAMPLE 2 (S1)**   | Run the Virtual Analyst on the Ostarine "Lynx Eyes" scenario.                |
| **EXAMPLE 3 (S4)**   | Run the Virtual Analyst on the Clomifene/Toremifene double-peak scenario.    |
| **EXPERIMENT METRICS** | Inspect the validation dashboard against the 549,096-interpretation reference set. |

A dark/light theme toggle is available in the top bar.

---

## Tech Stack

- **HTML5 / CSS3** — semantic structure, CSS custom properties for theming.
- **Vanilla JavaScript (ES6+)** — no framework, no transpiler, no bundler.
- **Canvas 2D API** — chromatogram traces (input signal, smoothed overlay, output, reference template).
- **Chart.js** — KPI dashboards in the Experiment Metrics tab.
- **CSS Grid + Flexbox** — responsive panel layout.
- **Embedded base64 PNG** — favicon and top-bar logo (no external image dependencies).

The entire application is a **single self-contained HTML file** that runs from the local filesystem with no build, no server, and no installation.

---

## Getting Started

### Requirements

- Any modern browser (Chrome, Firefox, Edge, Safari) released within the past two years.
- No installation, no Node.js, no Python, no compiler.

### Run

1. Clone or download this repository.
2. Open `VIRTUAL_ANALYST.html` directly in your browser:
   - Double-click the file, **or**
   - Drag-and-drop it onto a browser window, **or**
   - Use `File → Open` from the browser menu.
3. The application loads instantly and is fully interactive — no internet connection is required after the initial load (Chart.js is fetched from a CDN on first use).

### Recommended Workflow

1. Begin in the **METHOD** tab and review the default parameters for the three example substances.
2. Switch to **EXAMPLE 1** to observe the happy-path behaviour. Adjust the noise and amplitude sliders to see how robust the verdict is.
3. Move to **EXAMPLE 2** and try lowering the R² threshold; observe how the fuzzy panel shifts membership activations.
4. Move to **EXAMPLE 3** and disable Step 3 (set DTW threshold to 0); observe the false-negative scenario the fallback was designed to prevent.
5. Open **EXPERIMENT METRICS** to interpret the empirical results in light of the parameter changes you have just explored.

---

## File Structure

```
.
├── VIRTUAL_ANALYST.html    # The entire application (HTML + CSS + JS in one file)
└── README.md               # This document
```

---

## Author

**Guy Junior**

If you use this simulator in teaching, internal training, or method-development discussions, attribution to the author is appreciated.

---

## License

MIT License

Copyright (c) 2026 Guy Junior

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

> **Disclaimer.** This simulator is a pedagogical and method-development tool. It is **not** a validated diagnostic instrument and must not be used as the sole basis for any anti-doping adverse analytical finding. All real-world casework must follow the applicable WADA International Standard for Laboratories (ISL), Technical Documents, and the laboratory's accredited standard operating procedures.
