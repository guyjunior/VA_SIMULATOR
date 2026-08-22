# Virtual Analyst

> Automated chromatographic peak classification for **anti-doping analysis**
> (WADA prohibited-list substances), combining **R² goodness-of-fit under fuzzy
> logic**, **Dynamic Time Warping (DTW)** shape-similarity, and
> **ITP-Cutoff threshold rules** to produce per-substance verdicts on whether a
> peak should be flagged as suspect.

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status: Active](https://img.shields.io/badge/Status-Active-brightgreen.svg)
![Domain: Anti-Doping](https://img.shields.io/badge/Domain-Anti--Doping-blue.svg)
![Platform: Windows](https://img.shields.io/badge/Platform-Windows-blue.svg)

This repository holds the two artefacts of the same work - the decision logic,
twice, at two altitudes:

| | What it is | Where |
|---|---|---|
| **Method Simulator** | A single-file, in-browser simulator of the decision pipeline: transparent, tunable, with three worked scenarios and the full validation experiment (2,004 samples × 274 substances). Open it and explore the logic - nothing to install. | [`VIRTUAL ANALYST.html`](VIRTUAL%20ANALYST.html) |
| **Virtual Analyst application** | The pipeline for real: upload `.raw` files, convert to mzML, store spectra in ClickHouse, run the analysis against the CQCL/CQCL_reinj/CQN controls, review the chromatograms. Ships its converters, its screening method and a control-only example batch. | [`virtual-analyst/`](virtual-analyst/) |

The simulator answers *"why does the pipeline decide this way?"*; the
application answers *"what does it decide on my data?"*. Read Part I for the
science, Part II to run it.

---

## Table of Contents

**Part I - Method Simulator (in-browser)**

1. [Abstract](#abstract)
2. [Motivation and Scope](#motivation-and-scope)
3. [Application Architecture](#application-architecture)
4. [Theoretical Background](#theoretical-background)
5. [Pipeline Decision (4-Step Flow)](#pipeline-decision-4-step-flow)
6. [Method Configuration](#method-configuration)
7. [Worked Examples (Scenarios)](#worked-examples-scenarios)
8. [Interactive Analyst Controls](#interactive-analyst-controls)
9. [Experiment Metrics & Validation Results](#experiment-metrics--validation-results)
10. [User Interface Tabs](#user-interface-tabs)
11. [Tech Stack](#tech-stack)
12. [Getting Started](#getting-started)

**Part II - Virtual Analyst application (`virtual-analyst/`)**

13. [Platform](#platform)
14. [Application requirements](#application-requirements)
15. [Installing from scratch](#installing-from-scratch)
16. [Getting it running](#getting-it-running)
17. [The flow](#the-flow)
18. [Where things live](#where-things-live)
19. [What was ported, and what was adapted](#what-was-ported-and-what-was-adapted)
20. [When something goes wrong](#when-something-goes-wrong)
21. [Credits](#credits)

**Repository**

22. [Author](#author) · [License](#license)

---

# Part I - Method Simulator (in-browser)

## Abstract

The **Method Simulator - Virtual Analyst** is a self-contained, browser-based simulator that reproduces, in didactic and reproducible form, the decision logic of an automated chromatographic peak classifier deployed in routine anti-doping screening. The system addresses a recurring problem in liquid- and gas-chromatography mass-spectrometry (LC-MS / GC-MS) workflows: how to combine quantitative goodness-of-fit metrics (R²), shape-similarity measures (DTW), and concentration-based thresholds (ITP-Cutoff) into a coherent, interpretable, and tunable decision rule that minimises false negatives without overwhelming analysts with false positives.

The application provides three complementary views: (i) a **METHOD** tab for substance-by-substance parameter configuration; (ii) an **ANALYST** tab with three pre-loaded pedagogical scenarios drawn from the WADA prohibited list (S9 glucocorticoid Budesonide; S1 anabolic agent Ostarine; S4 hormone-and-metabolic-modulator Clomifene/Toremifene metabolites), each illustrating a distinct chromatographic challenge; and (iii) an **EXPERIMENT METRICS** tab presenting the empirical validation of the pipeline against an expert-reviewed reference set of 2,004 samples × 274 substances (549,096 total interpretations).

---

## Motivation and Scope

Routine anti-doping laboratories process thousands of samples per year, each screened against hundreds of prohibited substances. Fully manual peak review is neither scalable nor consistent across analysts. Conversely, naïve threshold-based automation tends to either (a) miss subtle peaks ("false negatives" - the most consequential error in this domain) or (b) flag too many irrelevant signals, increasing reviewer workload and degrading trust in the system.

This simulator was conceived to:

- **Make the decision logic transparent.** Each step of the pipeline is exposed in the UI with its inputs, thresholds, and intermediate values, allowing analysts and method developers to audit *why* a verdict was issued.
- **Permit substance-specific tuning.** Different chemical classes have different signal-to-noise characteristics, retention-time stability, and peak-shape morphology. The simulator allows per-substance adjustment of every threshold.
- **Provide reproducible benchmarks.** The bundled validation experiment (Section 9) reports a complete confusion matrix, allowing the community to assess the trade-offs of the chosen logic.
- **Support training.** The three pre-built scenarios are deliberately chosen to expose edge-case behaviours that frequently confuse less experienced reviewers.

---

## Application Architecture

The application is delivered as a **single HTML file** (`VIRTUAL ANALYST.html`) with no build step, no external dependencies beyond a CDN-served charting library, and no server component. All state is maintained in JavaScript memory; all rendering is performed via the HTML5 Canvas API (chromatogram traces) and Chart.js (metrics dashboards).

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

DTW measures the *minimum cumulative cost* of aligning two time-series of potentially different lengths or non-uniform local stretching. It is used here as a **shape-similarity fallback** when R² fails: a peak that is shifted in retention time, slightly broadened, or asymmetrically skewed may still match the expected morphology under a warped alignment, even if its point-wise R² collapses. The DTW score is normalised to a comparable scale and accepted when it falls **below** a configurable threshold (default range 0.020-0.120 depending on substance), denoting "shape close enough to remain suspect."

### Fuzzy Logic Classification

Fuzzy logic is applied **exclusively to the R² value** (a deliberate design choice highlighted in the UI's "fuzzy scope badge"). Five trapezoidal/triangular membership functions partition the R² domain:

| MF Name    | Linguistic Meaning           | Approx. R² Region |
|------------|------------------------------|-------------------|
| VERY LOW   | No fit - likely noise/blank  | R² < −2           |
| LOW        | Poor fit - interference      | −2 ≤ R² < 0       |
| MEDIUM     | Borderline                   | 0 ≤ R² < 0.6      |
| HIGH       | Good fit                     | 0.6 ≤ R² < 0.9    |
| VERY HIGH  | Excellent fit - suspect peak | R² ≥ 0.9          |

The classifier computes membership degrees, fires a small rule base, and produces a **defuzzified centroid** that drives the colour-coded R² bar in the Analyst panel. The rule firing strengths are exposed as horizontal bars so the user can audit the inference.

### ITP-Cutoff Concentration Rule

Even when R² and DTW agree that a peak is morphologically suspect, the substance-specific **ITP (International Threshold for Procedure) cutoff** acts as a final concentration guard: peaks below the cutoff (after integration) are reported as sub-threshold and not escalated, in line with WADA Technical Documents that define minimum reporting levels for each prohibited substance.

---

## Pipeline Decision (4-Step Flow)

Every interpretation flows through four sequential, auditable steps. The UI shows each step as a card that lights up green/red/skipped according to the live evaluation:

1. **Step 1 - Window Cut.**
   The signal is truncated to a substance-specific retention-time window centred on the expected PQC (Peak Quality Control) retention time. This eliminates upstream and downstream noise that would otherwise contaminate the R² calculation.

2. **Step 2 - R² Evaluation under Fuzzy Logic.**
   The cropped signal is correlated against the reference template; the resulting R² is fuzzified, and the defuzzified output determines whether the peak is provisionally accepted into the **fuzzy_r2 suspect zone**. If R² ≥ the configured threshold, the pipeline can short-circuit to Step 4. Otherwise, Step 3 is invoked as a fallback.

3. **Step 3 - DTW Threshold Check (Fallback).**
   When R² fails Step 2, the DTW distance between the observed and reference traces is computed. If DTW ≤ the configured threshold, the peak is rescued into the **fuzzy_dtw suspect zone** - recognising that shape similarity may persist where point-wise correlation does not. If both R² and DTW fail, the peak is classified as **negative** and the pipeline terminates.

4. **Step 4 - ITP-Cutoff Decision.**
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
| **Smoothing type**  | Savitzky-Golay, moving average, etc.                                     |
| **Smoothing value** | Window size or polynomial order for the smoother.                        |
| **Window cut**      | Width of the Step 1 crop around the PQC retention time.                  |
| **Concentration**   | Reference / spike concentration used to scale the simulated peak.        |
| **ITP-Cutoff**      | Minimum reportable concentration (Step 4 threshold).                     |
| **Has ITP-Cutoff?** | Toggle to disable Step 4 for substances without a defined MRPL.          |
| **R² threshold**    | Slider (default range −12.00 to 1.00) controlling Step 2 acceptance.     |
| **DTW threshold**   | Slider (default range 0.020-0.120) controlling Step 3 fallback.          |
| **Integration factor** | Calibration constant translating peak area to concentration.          |

Configuration changes propagate immediately to the active scenario, allowing real-time exploration of parameter sensitivity.

---

## Worked Examples (Scenarios)

The Analyst tab ships with three deliberately chosen examples, each illustrating a distinct chromatographic phenomenon:

### Example 1 - S9 Budesonide (Glucocorticoid)

A textbook well-resolved peak with high signal-to-noise. Illustrates the **happy path**: R² ≥ threshold, fuzzy classifier outputs HIGH/VERY HIGH, Step 3 is skipped, and Step 4 issues a clean POSITIVE verdict. Reference concentration: 30.00 (arb. units) with ITP-cutoff at 15.00.

### Example 2 - S1 Ostarine (SARM, "Lynx Eyes")

A low-amplitude, noise-embedded peak that requires careful smoothing and a permissive R² threshold to be recovered. Demonstrates the **fuzzy_r2 borderline regime** - small parameter changes flip the verdict, and the fuzzy panel shows competing membership activations. The "Lynx Eyes" nickname refers to the analyst skill traditionally required to spot such peaks manually.

### Example 3 - S4 Clomifene (4-OH)* / Toremifene (4-OH)* (Double Peak)

Co-eluting hydroxylated metabolites of two SERMs producing a **double-peak signature**. This scenario stresses the R² calculation (which compares against a single-peak template and therefore degrades) and showcases the **DTW fallback** (Step 3): the warped alignment recognises the morphology as still consistent with a suspect signature, rescuing what would otherwise be a false negative.

---

## Interactive Analyst Controls

The Analyst panel exposes a control strip with the following live sliders and toggles, all of which redraw the chromatogram and re-evaluate the pipeline on every change:

- **Noise** - additive white-noise amplitude.
- **Sigma (σ)** - Gaussian peak width.
- **Centroid shift** - retention-time offset of the peak relative to its expected position.
- **Amplitude** - peak height scaling.
- **Skew** - peak asymmetry (positive = tailing; negative = fronting).
- **Smoothing** - runtime smoother strength (independent of method-defined smoothing).
- **Crop toggle** - enable/disable the Step 1 window cut for visualisation.

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
| Recall         | **99.72%**| Of 2,116 truly suspect substances, 2,110 were detected - only 6 missed (FN).                                 |
| F1 Score       | **84.82%**| Harmonic mean of precision and recall.                                                                       |
| Specificity    | **99.86%**| 546,231 of 546,980 true negatives correctly classified.                                                      |
| FPR            | **0.14%** | Down from 0.32% when blank windows are included - false alarms are proportionally rare across the full set.  |
| AUC (ROC)      | **0.9979**| Near-ideal operating characteristic; rises from 0.9970 to 0.9979 once blank windows are factored in.         |

### Visualisations

The metrics tab includes:

- **Overview pie** - evaluated vs blank-window split of the 549,096 interpretations.
- **Comparison bar** - KPIs side-by-side under the two TN-inclusion regimes.
- **Composition pie** - TP/FP/TN/FN proportions.
- **TN breakdown bar** - substance-evaluated vs blank-window contributions to TN.
- **Top 20 TP chart** - substances most frequently and correctly flagged.
- **Top 20 FP chart** - substances driving the precision deficit (candidates for threshold re-tuning).
- **FN chart** - the six undetected suspect cases, by substance.
- **Bubble chart (Precision × Recall × FPR)** - per-substance operating points.
- **Per-substance metrics table** - full confusion matrix and KPIs for each of the 274 substances, with TN expanded by approximately 1,147 blank windows per substance.

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

- **HTML5 / CSS3** - semantic structure, CSS custom properties for theming.
- **Vanilla JavaScript (ES6+)** - no framework, no transpiler, no bundler.
- **Canvas 2D API** - chromatogram traces (input signal, smoothed overlay, output, reference template).
- **Chart.js** - KPI dashboards in the Experiment Metrics tab.
- **CSS Grid + Flexbox** - responsive panel layout.
- **Embedded base64 PNG** - favicon and top-bar logo (no external image dependencies).

The entire application is a **single self-contained HTML file** that runs from the local filesystem with no build, no server, and no installation.

---

## Getting Started

### Requirements

- Any modern browser (Chrome, Firefox, Edge, Safari) released within the past two years.
- No installation, no Node.js, no Python, no compiler.

### Run

1. Clone or download this repository.
2. Open `VIRTUAL ANALYST.html` directly in your browser:
   - Double-click the file, **or**
   - Drag-and-drop it onto a browser window, **or**
   - Use `File → Open` from the browser menu.
3. The application loads instantly and is fully interactive - no internet connection is required after the initial load (Chart.js is fetched from a CDN on first use).

### Recommended Workflow

1. Begin in the **METHOD** tab and review the default parameters for the three example substances.
2. Switch to **EXAMPLE 1** to observe the happy-path behaviour. Adjust the noise and amplitude sliders to see how robust the verdict is.
3. Move to **EXAMPLE 2** and try lowering the R² threshold; observe how the fuzzy panel shifts membership activations.
4. Move to **EXAMPLE 3** and disable Step 3 (set DTW threshold to 0); observe the false-negative scenario the fallback was designed to prevent.
5. Open **EXPERIMENT METRICS** to interpret the empirical results in light of the parameter changes you have just explored.

---

# Part II - Virtual Analyst application (`virtual-analyst/`)

Upload `.raw` files, convert them to mzML, store the spectra in ClickHouse and run the
anti-doping analytical pipeline: every sample is compared against the **CQCL** and the
**CQCL_reinj**, and corrected by the **CQN**.

A standalone app - it does not depend on the production API, MariaDB or ClickHouse. What
came from there is the **code**: the conversion, the XIC extraction and the chemistry
pipeline are faithful ports, not reimplementations.

No access control, by design. This is a bench tool, running on the machine of whoever does
the analysis. **Do not expose it on a network without putting authentication in front.**

## Platform

**This is a Windows application.** It was built, tested and is meant to run on
Windows 10 or 11. That is not a preference; it follows from the file formats:

| Piece | Windows-bound? | Why |
|---|---|---|
| **msconvert** (Waters, Agilent, Sciex, Bruker, Shimadzu) | **Yes** | The vendor readers are Windows DLLs. ProteoWizard has no native Linux build - on Linux it becomes a Docker container or a Wine prefix. |
| **ThermoRawFileParser** | No | There is a self-contained Linux build; swapping the folder in `external_tools/` is enough. |
| The pipeline, the web app, the databases | No | Plain Python, and ClickHouse runs in a container either way. |

So a **Thermo-only** workflow could be moved to Linux by replacing one folder. A
workflow that has to read any other vendor cannot, not without turning msconvert
into a container or a Wine prefix - an infrastructure decision, not a code one.

Every command in this document is for the Windows shell (`cmd` or PowerShell),
run from the project folder.

---

## Application requirements

| | Version | Notes |
|---|---|---|
| **Windows** | 10 or 11 | Developed on Windows 11 (build 26200) |
| **Python** | **3.12** | Built and tested on **3.12.10**. See below. |
| **Docker Desktop** | 4.x with the **WSL 2** backend | Tested with Docker Engine 29.4.1 and Compose v5.1.3 |
| **Disk** | ~2.5 GB, plus the data | 582 MB repository + ~500 MB virtual environment + ~1 GB ClickHouse image + the data volume |
| **Memory** | 8 GB is comfortable | ClickHouse is capped at 4 GB in `docker/clickhouse/memory.xml` |

### Why Python 3.12 specifically

The scientific stack (NumPy, SciPy, pandas, scikit-learn, DTAIDistance,
Pyteomics) is installed from binary wheels, and those wheels lag behind the
newest interpreter. **3.12 is the version this was built and tested on.** 3.13
and 3.14 are untested here and may not have a wheel for every dependency -
`pip` then tries to compile from source, which needs a C toolchain you probably
do not have.

If several Python versions are installed, the launcher picks the one you ask for:

```bat
py -3.12 -m venv venv
```

### Why Docker is needed

ClickHouse stores the spectra - tens of millions of rows per batch - and there
is no embedded version of it. The container in `docker-compose.yml` is the whole
installation: no local service, no configuration beyond that file. SQLite, by
contrast, comes with Python (`sqlite3` in the standard library) and needs nothing.

The **WSL 2 backend** is what Docker Desktop uses on Windows to run Linux
containers; the installer enables it and asks to reboot if the Windows feature is
missing.

### Disk, in detail

The ClickHouse volume lives on the Docker data drive, which is **C:** by default -
not necessarily where this project is. Measured here: an ~87 MB screening `.raw`
becomes ~11 MiB in ClickHouse after ZSTD, so 14 samples take roughly 150 MB.
`py -m tools.check_env` reports the current usage. If C: is tight, move the
Docker data root in Docker Desktop → Settings → Resources → Advanced.

---

## Installing from scratch

### 1. Python 3.12

Download it from <https://www.python.org/downloads/release/python-31210/> (the
*Windows installer, 64-bit*). During the install, tick **"Add python.exe to
PATH"**; the `py` launcher comes with it either way.

```bat
py -3.12 --version
```

It should print `Python 3.12.x`.

### 2. Docker Desktop

Download it from <https://www.docker.com/products/docker-desktop/> and install
with the default options, which include the WSL 2 backend. Reboot if it asks.
**Docker Desktop must be running** - its whale icon in the tray - before any
`docker` command works.

```bat
docker --version
docker compose version
```

### 3. The project

Clone this repository - the application lives in `virtual-analyst/`, and every
command below runs from inside that folder. The converters and the example data
ship with it, so there is nothing else to download.

---

## Getting it running

With Docker Desktop running, from the repository root:

```bat
cd virtual-analyst
docker compose up -d
py -3.12 -m venv venv
venv\Scripts\pip install -r requirements.txt
copy .env.example .env
venv\Scripts\python -m tools.seed_example
venv\Scripts\python run.py
```

Open <http://127.0.0.1:5000>.

What each line does:

| Command | What it does | How long |
|---|---|---|
| `docker compose up -d` | Starts ClickHouse on ports 8124/9001 | seconds (the first run pulls the image: ~1 GB) |
| `py -3.12 -m venv venv` | Creates the virtual environment | seconds |
| `pip install -r requirements.txt` | Installs Flask and the scientific stack | 1-3 min |
| `copy .env.example .env` | The configuration. It works unedited. | instant |
| `python -m tools.seed_example` | Schema, example method, example samples | ~1 min |
| `python run.py` | Starts the server on 127.0.0.1:5000 | instant |

### Checking the installation

```bat
venv\Scripts\python -m tools.check_env
```

It verifies, in order: both converters, the connection to ClickHouse, the tables,
the SQLite files, and how much disk the data is taking. It ends with
`Environment ready.` when nothing is missing, and names what is missing when
something is. Run it before every long batch - finding out that ClickHouse is
down after converting twenty samples is a bad afternoon.

**`seed_example` is the one command that matters.** It creates the schema, loads
the example method (274 substances and 7 internal standards, the real screening
method) and ingests the four example samples - one urine plus the CQCL, the
CQCL_reinj and the CQN it is judged against. It takes about a minute and leaves
the app with everything in place except the analysis itself, which is the part
worth doing by hand the first time:

> On **New analysis**, pick `Method Example`, tick the sample `Sample`, and press
> **Guess from names** to fill the three controls. Then **Process**.

It is safe to run again: it re-does only what is missing. `--force` rebuilds the
method and re-ingests the samples from scratch.

### The example data

The four `.raw` files ship in `examples/raw/`, so a clone can run the example
without hunting for anything. Point `--raw-dir` somewhere else to seed from a
different set.

**They are all control injections - no athlete sample is published here.** Three
of them play their own part (`CQCL.raw`, `CQCL_reinj.raw`, `CQN.raw`), and the
fourth, `Sample.raw`, is a further CQCL injection standing in for the sample
under analysis. Its internal header still identifies it as a CQCL, and that is
deliberate: a real urine cannot be redistributed, and a fabricated one would not
exercise the pipeline honestly.

It also makes the demonstration a better one. A fortified control analysed
against a fortified control lights up for the compounds it was spiked with, so
the example produces a screen with actual findings on it, instead of the wall of
`NEGATIVE` a clean urine would give you. What it is **not** is a validation of
the method: a control judged against itself says nothing about sensitivity or
specificity.

### What a clone weighs

The repository carries its own binaries and its own example data:

| | Size | Why it is versioned |
|---|---|---|
| `external_tools/` | 258 MB | Converting a `.raw` is step one; a clone that cannot convert is not runnable |
| `examples/raw/` | 324 MB | The example is only reproducible if the data comes with it |

That is ~582 MB across 559 files, and every clone pays it. No single file reaches
GitHub's 100 MB hard limit, though the four `.raw` files (80-83 MB each) are past
the 50 MB mark where GitHub warns - [Git LFS](https://git-lfs.com/) is the usual
answer if the history starts to hurt.

### Running the pieces by hand

`seed_example` is those three steps in order, and each one still works alone:

```bat
venv\Scripts\python -m app.db.schema           rem create the schema
venv\Scripts\python -m tools.import_method data\method_article.json --name "Method Example"
venv\Scripts\python -m tools.check_env         rem what is missing, and why
```

`check_env` is the shortcut for finding out what is missing before you lose a batch
halfway: a missing converter, ClickHouse down, the schema not created.

### The converters

They ship with the repository, in `external_tools/`, as portable builds. Nothing is
installed system-wide and nothing needs configuring:

```
external_tools/
├── thermo_raw_file_parser/   ThermoRawFileParser (Compomics + Thermo RawFileReader)
└── msconvert/                ProteoWizard portable
```

The paths are resolved from the project root, and there is deliberately **no environment
override**: pinning an absolute path is exactly what masked a broken default for months in
the sibling project - it only surfaced when moving to another machine.

### Ports

This project's ClickHouse comes up on **8124** (HTTP) and **9001** (native), not on the
default ports. The reason: the machine this was built on already ran another
ClickHouse on 8123/9000. The two coexist - and the same will be true on any
machine that already has one.

---

## The flow

### 1. Samples

Upload the batch's `.raw` files and pick a converter. Every file becomes a row in `samples`
with status `pending`; conversion and ingestion run in a **queue with a single worker**.

Serial on purpose: the converter already uses the whole CPU, and two simultaneous ingests
fight over ClickHouse memory - which answers that with an error (`MEMORY_LIMIT_EXCEEDED`),
not with slowness. Serial is faster overall and far more predictable.

State lives in the database, not in the worker's memory. If Flask dies halfway, the sample
stays marked as it was - visible, not a ghost.

#### Why Thermo goes to ThermoRawFileParser

A Thermo `.raw` marks certain centroids as **exception data**: lock mass, AGC overflow and
internal calibration. They are not analyte.

msconvert **does not honour that flag** - not even with `--filter "peakPicking true 1-"`.
Those signals land in the mzML, become points in ClickHouse, and become **false positives**
in the pipeline.

TRFP with `-x` uses Thermo's own `RawFileReader.dll` and honours the flag - the same
behaviour as Xcalibur and FreeStyle. That is why automatic mode sends Thermo to TRFP and
everything else to msconvert; picking msconvert for a Thermo `.raw` is legitimate for
comparing the two, and wrong for producing results. The screen warns you when you do it.

#### A file locked at the source

Thermo's `RawfileDataService.exe` keeps `.raw` files open **with write access** after
TraceFinder opens the batch, and never lets go. That breaks the conversion, because TRFP
reopens the file with `File.OpenRead()` to compute the header's SHA-1 - a method that
demands "no writers".

The conversion tries directly first and, **only if** the failure is a lock, copies the file
to `data/tmp_raw/` and converts the copy. The order matters: in the survey that motivated
this code, 5,918 of 6,081 files were not locked - always copying would move hundreds of GB
for nothing, and the attempt that fails costs about 2 s (the converter never even reads a
spectrum).

Not every failure is an error. A file **being acquired** (or with an aborted acquisition)
and an **empty** `.raw` become `skipped` with the reason on screen, not `error`.

### 2. Method

Substances, internal standards, groups and the IS gate. This is what the pipeline reads for
every sample.

Three fields usually explain "nothing came out":

| Field | The trap |
|---|---|
| `scan_filter` | It matches the text seen in the mzML **string for string**. One extra space and the substance finds no channel -> `INVALID`. That is why the field offers the filters actually observed in the samples already ingested. |
| `start_time`/`end_time` | The RT window, in minutes. Outside it there is no chromatogram. |
| `analysis_rule` | With no rule, it falls back to the standard one (R² + DTW). |

**IS gate**: the internal standard that locks everything. If it fails on a sample, the
pipeline aborts that whole sample (`blocked_by_is_gate`) and stores no results at all -
there is no point looking for substances in a run that did not work. The other internal
standards only raise alerts.

**Groups**: fragments of the same molecule. In competition the group is suspect if *any*
member is; out of competition, only if *all* of them are.

### Loading the method from the production system

The method can be built from the production analysis database instead of typed in
by hand. It is a two-step import, on purpose: this app never holds production
credentials and does not depend on a MySQL driver.

```bat
rem 1. export - run with the SOURCE system's interpreter (it has mysql-connector)
<source-repo>\venv\Scripts\python tools\export_method.py ^
    --env <source-repo>\.env --id-min 255 --id-max 535 --out data\method_article.json

rem 2. import (run with this project's venv)
venv\Scripts\python -m tools.import_method data\method_article.json --name Method
```

The export is **read-only** - it issues SELECTs and nothing else. The JSON it
writes is the reproducible artefact: it records the method name, the id window
and every parameter, so the same file rebuilds the same method anywhere. It
deliberately does **not** record the host or the database it came from - the
artefact is published, and the address of an internal server should not travel
with it.

`--id-min/--id-max` select which substances to take. The window is recorded in
the file, so the artefact always says what it contains.

Three mappings do the real work, and each one fails silently if it is wrong:

| Mapping | Why it matters |
|---|---|
| `type` lowercased | The source stores `SUBSTANCE`; the pipeline queries `type = 'substance'`. MariaDB ignores case, SQLite does not - hence `COLLATE NOCASE` on the column, and the normalisation on import. Get this wrong and every task returns zero results with no error. |
| `analysis_rule` by **code** | Rule ids happen to line up today. The code is what the pipeline dispatches on. |
| `scan_filter` copied verbatim | It has to equal the filter string seen in the mzML, character for character, or the substance finds no channel and comes out `INVALID` - which produces no row. The importer cross-checks against the filters actually observed in the ingested samples and lists any that were never seen. |

**Groups split by the id window.** A substance typed `SUBSTANCE_GROUP` is reached
only through its group. If the window keeps one member and leaves the other out,
the survivor would be imported, counted, and never analysed. The importer
promotes it to an individually analysed substance and says so - or leaves it out
entirely with `--drop-orphans`. There is no silent third option, because either
choice changes how that compound is judged.

### 3. Processing

Pick N samples and the three controls. Every sample is processed against all three at once:

| Control | Role |
|---|---|
| **CQCL** | the reference - where and how the substance shows up |
| **CQCL_reinj** | the second reference, used when the CQCL comes out `NEGATIVE` |
| **CQN** | the baseline, subtracted from the others before comparing |

Beyond the result, the task stores a **snapshot** of the method as it was at run time. The
processing reads the method live; the snapshot exists so the result stays interpretable
after someone edits a threshold.

Possible values in `sample_results.result`:

| Value | Meaning |
|---|---|
| `NEGATIVE` | Not suspect. |
| `SUSPECT_LOW` / `_MED` / `_HIGH` | The level decided by the fuzzy step, after the R²+DTW gate passed. |
| `SUSPECT_VERY_HIGH` | The sample peak is taller than the control's - a shortcut that skips the fuzzy step. |
| `NOT_INTERPRETED` | A visualisation-only substance (the `no_analyze` rule): extracted and plotted, never classified. |
| `INVALID` | It could not be evaluated (no peak in the CQCL, a window with fewer than 2 points). **It produces no row** - the substance simply does not appear. |

That matters when reading the screen: **zero results is not the same as everything
negative.** Check the per-sample outcome on the status screen.

### 4. Results

A per-sample table, with the chromatogram behind each row. The traces are plotted **already
corrected by the CQN** (point-by-point subtraction, matching on the nearest RT, with the
negative clamped at zero) - that is how the pipeline saw them when it decided. The raw data
is still in ClickHouse: the "no correction" button shows it.

`result` is what the machine decided; `result_validation` is what the person decided. Both
sit side by side - that is what makes it possible to measure where the method gets it wrong.

---

## Where things live

```
data/
├── uploads/<batch>/    .raw exactly as the browser sent it
├── mzml/               converter output
├── tmp_raw/            local copy (for a file locked at the source)
├── lake.sqlite         samples, channels, observed scan filters
└── analysis.sqlite     method, substances, tasks, results
log/
├── runner.log          the conversion/ingestion queue and the processing
└── <batch>.log         the pipeline log, one per batch
```

### Why two SQLite files

They mirror the split the production system has between two databases. Not a flourish: the
two halves have tables with the **same name and different columns** - `scan_filters` in the
lake holds the text observed in the mzML (`scan_filter`), the one in the analysis database
holds the name configured in the method (`name`). Merging them into a single file would
force a table rename and a divergence from the production code, which is precisely what
this project is trying not to do.

### What goes into ClickHouse

The same as the production system, same engine and same `ORDER BY`:

| Database | Table | What it is |
|---|---|---|
| `lake` | `spectra_points` | the centroids. The big table. `ORDER BY (channel_id, mz, scan_index)` is what makes the XIC cheap: the m/z window becomes a contiguous range scan. |
| `lake` | `scans` | the run's RT grid. The `LEFT JOIN` against it is what turns a scan with no point in the window into a **zero** instead of a hole - a continuous XIC, the way Thermo does it. |
| `lake` | `chromatogram_points` | the TIC/BasePeak traces that already ship inside the mzML. |
| `analysis` | `processing_sampleparameter` | the 4 traces of every substance x sample, exactly as the pipeline saw them. |

Order of magnitude measured here: an ~87 MB screening `.raw` produced ~5,900 scans and ~1.4
million centroids, taking ~11 MiB in ClickHouse after ZSTD. The Docker volume lives on the
system drive - `py -m tools.check_env` shows how much is in use.

---

## What was ported, and what was adapted

| File | Origin | What changed |
|---|---|---|
| `app/services/virtual_analyst.py` | the production pipeline | Only the plumbing: MariaDB -> SQLite, this project's ClickHouse. **The chemistry was not touched.** |
| `app/services/data_provider.py` | the production `DbDataProvider` | Same. The calibration against MSFileReader (the RT edge tolerance, the smoothing sigma, the `LEFT JOIN` that zeroes an empty scan) was kept character for character. |
| `app/services/build_method_snapshot.py` | same | Only the connection. |
| `app/services/convert.py` | the production ingestion | It gained the on-screen converter choice (production decides by vendor, with no option). |
| `app/services/ingest.py` | same | Rewritten around the "uploaded file" path instead of a folder sweep; the writing itself is identical. |

`app/db/sqlite.py` is what makes that port possible: it teaches SQLite the dialect the
production code speaks (`cursor(dictionary=True)`, the `%s` placeholder, `NOW()`). The less
the pipeline changes, the smaller the chance its results diverge from the real system.

**A trap recorded there**: placeholder translation is a textual replace of `%s` with `?`.
That would break a literal containing `%s` - `LIKE '%s%'`, for instance. There is none today
(the pipeline's `LIKE` clauses use `'%SUSPECT%'`, with an uppercase S). The way out is to
pass the literal as a **parameter**, never inline in the SQL.

---

## When something goes wrong

| Symptom | Where to look |
|---|---|
| The task finished with 0 results | The task status screen. `blocked_by_is_gate` aborts the sample before any substance is analysed. |
| A substance never appears | The `scan_filter` differs from the one observed in the mzML, or the RT window is wrong -> `INVALID`, which produces no row. |
| `MEMORY_LIMIT_EXCEEDED` during ingestion | The ceilings live in `docker/clickhouse/memory.xml`. The insert already retries with backoff for that specific error. |
| Conversion fails with "being used by another process" | A handle left open by `RawfileDataService`. The code already retries from a local copy; if it persists, close the batch in TraceFinder. |
| A peak Xcalibur does not show | The sample was converted with msconvert. Reprocess it with ThermoRawFileParser (the screen has a per-sample button). |

### Installation and Windows

| Symptom | Cause and fix |
|---|---|
| `docker: error during connect` or `the system cannot find the file` | Docker Desktop is not running. Start it and wait for the whale icon to settle. |
| `Ports are not available: 8124` | Something else holds the port. `docker ps` shows what; change the host side of the mapping in `docker-compose.yml` and `CH_PORT` in the `.env` to match. |
| `pip` tries to compile and fails | The interpreter is not 3.12, so there is no matching wheel. Delete `venv\` and recreate it with `py -3.12 -m venv venv`. |
| `venv\Scripts\activate` is blocked in PowerShell | The execution policy. You do not need to activate anything - every command here calls `venv\Scripts\python` directly, on purpose. |
| `ModuleNotFoundError: No module named 'app'` | The command was run from another folder, or without `-m`. Run it from the project root, as `venv\Scripts\python -m tools.<name>`. |
| `ThermoRawFileParser not found` | The clone is missing `external_tools/` - it may have been cloned without the large files. |
| ClickHouse is up but `check_env` cannot reach it | `CH_PORT` in the `.env` must match the host side in `docker-compose.yml` (8124 by default, not 8123). |

---

## Credits

This project is a thin layer over other people's work. The scientific tools below
do the heavy lifting; if you publish results produced with it, cite them.

| Tool | What it does here | Licence |
|---|---|---|
| [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) | Converts Thermo `.raw` to mzML, honouring the exception-data flag | Apache-2.0 |
| [ProteoWizard / msconvert](https://proteowizard.sourceforge.io/) | Converts every non-Thermo vendor format to mzML | Apache-2.0 |
| [Pyteomics](https://github.com/levitsky/pyteomics) | Reads the mzML: spectra, scans, chromatograms | Apache-2.0 |
| [ClickHouse](https://clickhouse.com/) | Stores the centroids and serves the XIC queries | Apache-2.0 |
| [NumPy](https://numpy.org/) · [SciPy](https://scipy.org/) · [pandas](https://pandas.pydata.org/) | Numerical core: arrays, gaussian smoothing, dataframes | BSD-3-Clause |
| [scikit-learn](https://scikit-learn.org/) | R² between sample and control (`r2_score`) | BSD-3-Clause |
| [DTAIDistance](https://github.com/wannesm/dtaidistance) | Dynamic Time Warping between the normalised curves | Apache-2.0 |
| [Flask](https://flask.palletsprojects.com/) | The web layer | BSD-3-Clause |
| [Plotly.js](https://plotly.com/javascript/) | Chromatogram plots | MIT |
| [SQLite](https://sqlite.org/) | Method, tasks and results - via Python's stdlib `sqlite3` | Public domain |

### References

- **ThermoRawFileParser** - Hulstaert N., Shofstahl J., Sachsenberg T., Walzer M.,
  Barsnes H., Martens L., Perez-Riverol Y. *ThermoRawFileParser: Modular, Scalable,
  and Cross-Platform RAW File Conversion.* J. Proteome Res. 2020, 19(1), 537-542.
  <https://doi.org/10.1021/acs.jproteome.9b00328>
- **ProteoWizard** - Chambers M.C. *et al.* *A cross-platform toolkit for mass
  spectrometry and proteomics.* Nat. Biotechnol. 2012, 30, 918-920.
  <https://doi.org/10.1038/nbt.2377> · Kessner D., Chambers M., Burke R., Agus D.,
  Mallick P. *ProteoWizard: open source software for rapid proteomics tools
  development.* Bioinformatics 2008, 24(21), 2534-2536.
  <https://doi.org/10.1093/bioinformatics/btn323>
- **Pyteomics** - Goloborodko A.A., Levitsky L.I., Ivanov M.V., Gorshkov M.V.
  *Pyteomics - a Python framework for exploratory data analysis and rapid software
  prototyping in proteomics.* J. Am. Soc. Mass Spectrom. 2013, 24(2), 301-304.
  <https://doi.org/10.1007/s13361-012-0516-6> · Levitsky L.I., Klein J.A.,
  Ivanov M.V., Gorshkov M.V. *Pyteomics 4.0: Five Years of Development of a Python
  Proteomics Framework.* J. Proteome Res. 2019, 18(2), 709-714.
  <https://doi.org/10.1021/acs.jproteome.8b00717>
- **NumPy** - Harris C.R. *et al.* *Array programming with NumPy.* Nature 2020, 585,
  357-362. <https://doi.org/10.1038/s41586-020-2649-2>
- **SciPy** - Virtanen P. *et al.* *SciPy 1.0: fundamental algorithms for scientific
  computing in Python.* Nat. Methods 2020, 17, 261-272.
  <https://doi.org/10.1038/s41592-019-0686-2>
- **pandas** - McKinney W. *Data Structures for Statistical Computing in Python.*
  Proc. 9th Python in Science Conf. 2010, 56-61.
  <https://doi.org/10.25080/Majora-92bf1922-00a>
- **scikit-learn** - Pedregosa F. *et al.* *Scikit-learn: Machine Learning in
  Python.* J. Mach. Learn. Res. 2011, 12, 2825-2830.
- **DTAIDistance** - Meert W., Hendrickx K., Van Craenendonck T., Robberechts P.,
  Blockeel H., Davis J. *DTAIDistance.* Zenodo - cite the record matching the
  version you ran.

> Check every DOI and version against what is actually installed
> (`venv\Scripts\pip freeze`) before submitting. The pinned versions live in
> `requirements.txt`; the ones above are the canonical papers, not
> version-specific records.

### Vendor licences

`external_tools/` is **redistributed with this repository**, and those binaries
carry licences that are not this project's:

- **ThermoRawFileParser** is Apache-2.0, but it embeds Thermo Fisher's
  `RawFileReader` assemblies, which are covered by Thermo's own licence terms.
- **ProteoWizard** is Apache-2.0 at its core, but the vendor readers it bundles
  (Waters, Agilent, Sciex, Bruker, Shimadzu, Thermo) are proprietary libraries
  distributed under per-vendor agreements - the ones you accept when you download
  a build, which is why the official container is named
  `pwiz-skyline-i-agree-to-the-vendor-licenses`.

Confirm the redistribution terms of both before making the repository public or
attaching it to a publication. Shipping them is convenient; whether you may is a
question for the licences, not for this README.

---

## Author

**Guy Junior**

If you use this work - the simulator or the application - in teaching, internal training, or method-development discussions, attribution to the author is appreciated.

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

> **Disclaimer.** The simulator and the application are pedagogical and method-development tools. Neither is a validated diagnostic instrument, and neither must be used as the sole basis for any anti-doping adverse analytical finding. All real-world casework must follow the applicable WADA International Standard for Laboratories (ISL), Technical Documents, and the laboratory's accredited standard operating procedures.
