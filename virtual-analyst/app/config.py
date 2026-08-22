"""Single source of configuration - everything from the environment lands here.

Paths are anchored to the PROJECT ROOT, never to the process cwd. Pinning an
absolute converter path in the .env is exactly what masked a broken default for
months in the sibling project: it only broke when moving to another machine.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

from dotenv import load_dotenv

# The ported pipeline logs with emoji, and on Windows a redirected stdout
# assumes cp1252 - which cannot encode emoji and killed an analysis mid-batch
# with UnicodeEncodeError. This never showed up in production because the API
# runs on Linux, where UTF-8 is the default. Fix it here, in the one module
# everything imports.
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

ROOT = Path(__file__).resolve().parents[1]

load_dotenv(ROOT / ".env")


def _env(name: str, default: str = "") -> str:
    value = os.getenv(name)
    return default if value is None or value.strip() == "" else value.strip()


# --- Directories ------------------------------------------------------------
DATA_DIR = ROOT / "data"
UPLOAD_DIR = DATA_DIR / "uploads"     # .raw exactly as the browser sent it
MZML_DIR = DATA_DIR / "mzml"          # converter output
TMP_RAW_DIR = DATA_DIR / "tmp_raw"    # local copy, used when the original is locked
LOG_DIR = ROOT / "log"
EXTERNAL_TOOLS = ROOT / "external_tools"

for _directory in (UPLOAD_DIR, MZML_DIR, TMP_RAW_DIR, LOG_DIR):
    _directory.mkdir(parents=True, exist_ok=True)

# --- SQLite -----------------------------------------------------------------
# Two files, mirroring the two databases of the production system. Not a
# flourish: the two halves have tables with the SAME NAME and different columns
# (`scan_filters` in the lake holds the filter text seen in the mzML; the one in
# the analysis database holds the name configured in the method). Merging them
# into a single file would force a table rename and a divergence from the
# production code.
LAKE_DB = DATA_DIR / "lake.sqlite"
ANALYSIS_DB = DATA_DIR / "analysis.sqlite"

# --- ClickHouse -------------------------------------------------------------
CH_HOST = _env("CH_HOST", "localhost")
CH_PORT = int(_env("CH_PORT", "8124"))
CH_USER = _env("CH_USER", "va")
CH_PASSWORD = _env("CH_PASSWORD", "va")
CH_LAKE_DB = _env("CH_LAKE_DB", "lake")
CH_ANALYSIS_DB = _env("CH_ANALYSIS_DB", "analysis")

# --- Converters -------------------------------------------------------------
# The builds ship with the repository, so the path is resolved from the project
# root and nothing needs configuring. There is deliberately no environment
# override: pinning an absolute path is exactly what masked a broken default for
# months in the sibling project - it only broke when moving to another machine.
_TRFP_NAMES = ("ThermoRawFileParser.exe", "ThermoRawFileParser")


def _thermo_raw_parser() -> Path:
    """Windows ships an .exe, the Linux self-contained build has no suffix."""
    base = EXTERNAL_TOOLS / "thermo_raw_file_parser"
    for name in _TRFP_NAMES:
        if (base / name).exists():
            return base / name
    return base / _TRFP_NAMES[0]


THERMO_RAW_PARSER_PATH = _thermo_raw_parser()
MSCONVERT_PATH = EXTERNAL_TOOLS / "msconvert" / "msconvert.exe"

# --- Flask ------------------------------------------------------------------
SECRET_KEY = _env("FLASK_SECRET_KEY", "dev-do-not-use-in-production")
DEBUG = _env("FLASK_DEBUG", "1") == "1"
MAX_UPLOAD_MB = int(_env("MAX_UPLOAD_MB", "8192"))

# --- Pipeline ---------------------------------------------------------------
_tolerance = _env("MASS_TOLERANCE_PPM")
MASS_TOLERANCE_PPM = float(_tolerance) if _tolerance else None

# Batch sizes for the ClickHouse inserts. Inherited from the production
# ingestion: large enough for the insert to be efficient, small enough that the
# row list does not blow up the Python process memory.
BATCH_SCANS = 50_000
BATCH_POINTS = 500_000
