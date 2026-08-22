"""raw -> mzML.

A port of the conversion validated in production. What was kept word for word
are the details that already cost time:

  * TRFP v2.x requires `-b` with an ABSOLUTE output path. A relative path breaks
    when the subprocess cwd is not the expected one.
  * In that same version, `-p` means NO peak picking - the opposite of what it
    looks like. Peak picking is already the default; pass nothing.
  * `-x` (--excludeExceptionData) is the reason TRFP is here at all. See
    `pick_converter` below.
  * A failed conversion is not always an error. Three classes need different
    handling, and treating them as one fills the log with noise.

DIFFERENCE FROM PRODUCTION: there the converter is chosen by vendor, with no
option. Here the user picks it on screen, because this app also exists to
compare the two. The default is still the right one.
"""

from __future__ import annotations

import shutil
import subprocess
from contextlib import contextmanager
from pathlib import Path

from app import config

# --- failure classes --------------------------------------------------------
LOCKED = "locked"              # handle open at the source -> retry from a copy
ACQUIRING = "acquiring"        # run in progress OR aborted -> skip
EMPTY = "empty"                # .raw with header only -> skip
UNKNOWN = "unknown"            # a real error -> report it

# Signatures in the converter's stdout. The last two come from TRFP; msconvert
# does not expose them, but the lock signature comes from Windows itself and
# applies to both.
_SIGNATURES = (
    ("being used by another process", LOCKED),
    ("still being acquired", ACQUIRING),
    ("Empty RAW file", EMPTY),
)

SKIP_REASON = {
    ACQUIRING: (
        "File is being acquired, or its acquisition was aborted "
        "(the flag lives in the .raw header and is not cleared)"
    ),
    EMPTY: "The .raw file has no spectra (header only)",
}

# --- converters offered on screen -------------------------------------------
TRFP = "thermo_raw_file_parser"
MSCONVERT = "msconvert"
AUTO = "auto"

CONVERTERS = (AUTO, TRFP, MSCONVERT)


class ConversionFailure(RuntimeError):
    """The conversion failed, with the cause already classified."""

    def __init__(self, message: str, *, kind: str, output: str = ""):
        super().__init__(message)
        self.kind = kind
        self.output = output


def classify_failure(output: str) -> str:
    """Read the converter output and say what kind of failure it is.

    Substring matching, not a regex: the messages come from TRFP and from
    Windows, they are stable, and a regex here would only invent a new way to be
    wrong.
    """
    text = output or ""
    for marker, kind in _SIGNATURES:
        if marker in text:
            return kind
    return UNKNOWN


def detect_vendor(input_path: Path) -> str:
    """Identify the vendor from the acquisition file or folder.

    The rule that matters: a Thermo `.raw` is a FILE and a Waters `.raw` is a
    FOLDER - same extension, different vendors.
    """
    suffix = input_path.suffix.lower()
    is_file = input_path.is_file()
    is_dir = input_path.is_dir()

    if is_file and suffix == ".raw":
        return "thermo"
    if is_dir and suffix == ".raw":
        return "waters"
    if is_dir and suffix == ".d":
        if (input_path / "analysis.tdf").exists() or (input_path / "analysis.baf").exists():
            return "bruker"
        return "agilent"
    if is_file and suffix in (".wiff", ".wiff2"):
        return "sciex"
    if is_file and suffix == ".lcd":
        return "shimadzu"
    if is_file and suffix in (".mzml", ".mzxml", ".mzmlb", ".mgf"):
        return "mzml"
    return "unknown"


def pick_converter(vendor: str, choice: str) -> str:
    """Resolve `auto` and return either TRFP or MSCONVERT.

    Under `auto`: Thermo goes to TRFP, everything else to msconvert.

    WHY THERMO DOES NOT GO TO MSCONVERT
    A Thermo `.raw` marks certain centroids as EXCEPTION DATA - lock mass, AGC
    overflow, internal calibration. They are not analyte. msconvert does not
    honour that flag, not even with `--filter "peakPicking true 1-"`: those
    signals land in the mzML, become points in ClickHouse, and become FALSE
    POSITIVES in the pipeline. TRFP with `-x` uses Thermo's own
    RawFileReader.dll and honours the flag - the same behaviour as Xcalibur and
    FreeStyle.

    That is why the screen warns when the user picks msconvert for a Thermo
    .raw: it is a legitimate choice for comparing the two, and a wrong one for
    producing results.
    """
    if choice == TRFP:
        return TRFP
    if choice == MSCONVERT:
        return MSCONVERT
    return TRFP if vendor == "thermo" else MSCONVERT


def converter_warning(vendor: str, converter: str) -> str | None:
    if vendor == "thermo" and converter == MSCONVERT:
        return (
            "msconvert on a Thermo .raw keeps the exception data (lock mass, AGC "
            "overflow, internal calibration) in the mzML. That produces peaks "
            "Xcalibur does not show and can turn into false positives. Use "
            "ThermoRawFileParser."
        )
    if vendor != "thermo" and converter == TRFP:
        return (
            f"ThermoRawFileParser only reads Thermo .raw files; this one was detected "
            f"as '{vendor}'. The conversion will fail - use msconvert."
        )
    return None


def run_thermorawfileparser(
    input_path: Path,
    out_mzml_dir: Path,
    *,
    parser_path: Path | None = None,
    exclude_exception_data: bool = True,
) -> Path:
    """Convert a Thermo .raw to mzML with ThermoRawFileParser."""
    parser_path = parser_path or config.THERMO_RAW_PARSER_PATH
    if not Path(parser_path).exists():
        raise ConversionFailure(
            f"ThermoRawFileParser not found at {parser_path}. It ships with the "
            f"repository - a clone that is missing it was probably made without "
            f"Git LFS, or the file was excluded.",
            kind=UNKNOWN,
        )

    out_mzml_dir = out_mzml_dir.resolve()
    out_mzml_dir.mkdir(parents=True, exist_ok=True)
    mzml_path = (out_mzml_dir / f"{input_path.stem}.mzML").resolve()

    if mzml_path.exists():
        mzml_path.unlink()

    args = [
        str(parser_path),
        "-i", str(input_path.resolve()),
        "-b", str(mzml_path),   # ABSOLUTE - a relative path breaks
        "-f", "2",              # 2 = indexed mzML
    ]
    if exclude_exception_data:
        args.append("-x")

    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        output = f"{result.stdout}\n{result.stderr}"
        raise ConversionFailure(
            f"ThermoRawFileParser failed (exit {result.returncode}) on {input_path.name}\n"
            f"  args:   {args}\n"
            f"  STDERR:\n{result.stderr.strip()}\n"
            f"  STDOUT:\n{result.stdout.strip()}",
            kind=classify_failure(output),
            output=output,
        )
    return mzml_path


def run_msconvert(input_path: Path, out_mzml_dir: Path, *, msconvert_path: Path | None = None) -> Path:
    """Convert any supported format to mzML via msconvert."""
    msconvert_path = msconvert_path or config.MSCONVERT_PATH
    if not Path(msconvert_path).exists():
        raise ConversionFailure(
            f"msconvert not found at {msconvert_path}. It ships with the "
            f"repository - a clone that is missing it was probably made without "
            f"Git LFS, or the file was excluded.",
            kind=UNKNOWN,
        )

    out_mzml_dir = out_mzml_dir.resolve()
    out_mzml_dir.mkdir(parents=True, exist_ok=True)
    mzml_path = out_mzml_dir / f"{input_path.stem}.mzML"

    if mzml_path.exists():
        mzml_path.unlink()

    args = [
        str(msconvert_path),
        str(input_path.resolve()),
        "--mzML",
        "--filter", "peakPicking true 1-",
        "--outdir", str(out_mzml_dir),
        "--outfile", mzml_path.name,
    ]
    # No check=True on purpose: the output is what classifies the failure, and
    # CalledProcessError would throw it away.
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        output = f"{result.stdout}\n{result.stderr}"
        raise ConversionFailure(
            f"msconvert failed (exit {result.returncode}) on {input_path.name}\n"
            f"  args:   {args}\n"
            f"  STDERR:\n{result.stderr.strip()}\n"
            f"  STDOUT:\n{result.stdout.strip()}",
            kind=classify_failure(output),
            output=output,
        )
    return mzml_path


@contextmanager
def temporary_local_copy(source: Path, target_dir: Path | None = None):
    """Copy `source` to a local folder and delete it on the way out.

    This exists because of the lock at the source: Thermo's `RawfileDataService`
    keeps .raw files open WITH WRITE ACCESS after TraceFinder opens the batch,
    and never lets go. The converter uses `File.OpenRead()`, which demands "no
    writers", and is refused - even for a file untouched for weeks.

    There is no risk of partial data: if the .raw really is being written, the
    copy comes out truncated and the converter REFUSES it - the format's index
    lives at the end of the file, so its absence is detected immediately.
    """
    target_dir = target_dir or config.TMP_RAW_DIR
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / source.name
    try:
        if source.is_dir():
            if target.exists():
                shutil.rmtree(target)
            shutil.copytree(source, target)
        else:
            shutil.copy2(source, target)
        yield target
    finally:
        try:
            if target.is_dir():
                shutil.rmtree(target, ignore_errors=True)
            elif target.exists():
                target.unlink()
        except OSError:
            pass


def convert_to_mzml(input_path: Path, out_mzml_dir: Path, choice: str = AUTO) -> tuple[Path, str]:
    """Convert, returning (mzML path, converter used).

    It tries directly first and only copies to the local disk IF the failure is
    a lock. The order is deliberate: in the survey that motivated this code,
    5,918 of 6,081 files were not locked - always copying would move hundreds of
    GB for nothing. And the attempt that fails is cheap (median 2 s), because
    the converter never gets as far as reading a spectrum.
    """
    input_path = Path(input_path)
    vendor = detect_vendor(input_path)

    if vendor == "mzml":
        # Already converted: it just moves into the working folder.
        out_mzml_dir.mkdir(parents=True, exist_ok=True)
        target = out_mzml_dir / input_path.name
        if target.resolve() != input_path.resolve():
            shutil.copy2(input_path, target)
        return target, "copy"

    converter = pick_converter(vendor, choice)
    runner = run_thermorawfileparser if converter == TRFP else run_msconvert

    try:
        return runner(input_path, out_mzml_dir), converter
    except ConversionFailure as failure:
        if failure.kind != LOCKED:
            raise
        with temporary_local_copy(input_path) as copy:
            return runner(copy, out_mzml_dir), converter
