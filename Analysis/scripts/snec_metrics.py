#!/usr/bin/env python3
"""Shared utilities for SNEC speed/accuracy sweep analysis."""

from __future__ import annotations

from pathlib import Path
import atexit
import csv
import math
import re
import shutil
import subprocess
from typing import Dict, Iterable, List, Tuple

import numpy as np

DAY = 86400.0


def read_lum(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    rows: List[Tuple[float, float]] = []
    with path.open() as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            t, lum = float(parts[0]), float(parts[1])
            if lum > 0.0:
                rows.append((t, lum))
    arr = np.array(rows, dtype=float)
    if arr.size == 0:
        raise ValueError(f"No valid luminosity rows in {path}")
    return arr[:, 0], arr[:, 1]


def compute_mag_metrics(
    ref_lum_path: Path,
    run_lum_path: Path,
    ngrid: int = 2500,
) -> Dict[str, float]:
    t_ref, l_ref = read_lum(ref_lum_path)
    t_run, l_run = read_lum(run_lum_path)

    t_min = max(t_ref[0], t_run[0], 0.1 * DAY)
    t_max = min(t_ref[-1], t_run[-1])
    if t_max <= t_min:
        raise ValueError("No overlapping time domain between reference and run")

    t = np.linspace(t_min, t_max, ngrid)
    d = t / DAY
    l_ref_i = np.interp(t, t_ref, l_ref)
    l_run_i = np.interp(t, t_run, l_run)
    dmag = -2.5 * np.log10(l_run_i / l_ref_i)

    def rms(mask: np.ndarray) -> float:
        if not np.any(mask):
            return float("nan")
        x = dmag[mask]
        return float(np.sqrt(np.mean(x * x)))

    rms_0_5 = rms((d >= 0.0) & (d <= 5.0))
    rms_5_30 = rms((d > 5.0) & (d <= 30.0))
    rms_30_100 = rms((d > 30.0) & (d <= 100.0))
    rms_100_200 = rms((d > 100.0) & (d <= 200.0))
    rms_200_300 = rms((d > 200.0) & (d <= 300.0))
    rms_gt30 = rms(d > 30.0)
    rms_gt5 = rms(d > 5.0)
    rms_all = rms(d >= 0.0)
    max_abs = float(np.max(np.abs(dmag)))
    n_pass = int(rms_0_5 <= 0.10) + int(rms_5_30 <= 0.10) + int(rms_gt30 <= 0.10)

    return {
        "rms_0_5": rms_0_5,
        "rms_5_30": rms_5_30,
        "rms_30_100": rms_30_100,
        "rms_100_200": rms_100_200,
        "rms_200_300": rms_200_300,
        "rms_gt30": rms_gt30,
        "rms_gt5": rms_gt5,
        "rms_all": rms_all,
        "max_abs_mag": max_abs,
        "n_pass": float(n_pass),
    }


# ---------------------------------------------------------------------------
# Complete baseline parameter set — every sweep MUST start from this.
# Loaded from Analysis/scripts/baseline_parameters (single source of truth)
# so that sweep results never depend on the user's working parameters file.
# ---------------------------------------------------------------------------

_BASELINE_FILE = Path(__file__).resolve().parent / "baseline_parameters"


def _load_baseline(path: Path) -> Dict[str, str]:
    """Parse the canonical baseline_parameters file into a key→value dict."""
    params: Dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, _, value = line.partition("=")
        params[key.strip()] = value.strip()
    return params


BASELINE_PARAMETERS: Dict[str, str] = _load_baseline(_BASELINE_FILE)


_parameters_backup: Dict[str, Path] = {}


def write_parameters_file(
    path: Path,
    overrides: Dict[str, str] | None = None,
) -> None:
    """Write a complete, self-contained parameters file from the baseline template.

    Every parameter is written explicitly — nothing depends on the user's
    parameters file.  The original file is backed up on the first call and
    automatically restored on interpreter exit via ``restore_parameters_file``.
    """
    path = Path(path).resolve()
    key = str(path)
    if key not in _parameters_backup and path.exists():
        backup = path.with_suffix(".bak")
        shutil.copy2(path, backup)
        _parameters_backup[key] = backup
        atexit.register(restore_parameters_file, path)

    params = dict(BASELINE_PARAMETERS)
    if overrides:
        params.update(overrides)
    lines = [f"{key} = {value}" for key, value in params.items()]
    path.write_text("\n".join(lines) + "\n")


def restore_parameters_file(path: Path) -> None:
    """Restore the original parameters file from backup."""
    path = Path(path).resolve()
    key = str(path)
    backup = _parameters_backup.get(key)
    if backup is None:
        backup = path.with_suffix(".bak")
    if backup.exists():
        shutil.copy2(backup, path)
        backup.unlink()
        _parameters_backup.pop(key, None)


def run_cmd(
    cmd: List[str],
    cwd: Path,
    log_path: Path,
    timeout_s: float | None = None,
) -> Tuple[int, bool]:
    """Run a command, capturing output to log_path.  Returns (returncode, timed_out)."""
    timed_out = False
    rc = 0
    with log_path.open("w") as f:
        try:
            proc = subprocess.run(
                cmd,
                cwd=cwd,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=False,
                timeout=timeout_s,
            )
            rc = int(proc.returncode)
        except subprocess.TimeoutExpired:
            timed_out = True
            rc = -999
            f.write(f"\nTIMEOUT: command exceeded {timeout_s:.1f} s\n")
    return rc, timed_out


def parse_log_indicators(log_text: str) -> Dict[str, float]:
    """Extract diagnostic counters from an SNEC log file."""
    eos_problems = len(re.findall(r"EOS problem", log_text))
    scratch_steps = len(re.findall(r"Scratching entire step", log_text))
    saha_warnings = len(re.findall(r"convergence problem in saha solver", log_text))
    reached_tend = 1.0 if "Done! :-) tend reached" in log_text else 0.0
    ni_eval_total = 0.0
    ni_eval_remap = 0.0
    ni_match = re.search(r"Ni evals:\s+(\d+).*remap-forced resets:\s+(\d+)", log_text)
    if ni_match:
        ni_eval_total = float(ni_match.group(1))
        ni_eval_remap = float(ni_match.group(2))
    remap_count = 0.0
    remap_match = re.search(r"Grid remaps:\s+(\d+)", log_text)
    if remap_match:
        remap_count = float(remap_match.group(1))
    return {
        "eos_problems": float(eos_problems),
        "scratch_steps": float(scratch_steps),
        "saha_warnings": float(saha_warnings),
        "reached_tend": reached_tend,
        "ni_eval_total": ni_eval_total,
        "ni_eval_remap": ni_eval_remap,
        "remap_count": remap_count,
    }


def count_nan_in_lum(path: Path) -> int:
    """Count rows where luminosity is NaN in a lum_observed.dat file."""
    if not path.exists():
        return 0
    n = 0
    with path.open() as f:
        for line in f:
            s = line.split()
            if len(s) >= 2 and s[1].lower() == "nan":
                n += 1
    return n


def load_tsv(path: Path) -> List[Dict[str, str]]:
    """Read a TSV file written by write_tsv and return a list of row dicts."""
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# Alias for clarity in table-generation code.
read_tsv = load_tsv


def to_float(x: str) -> float:
    """Parse a string to float, stripping whitespace."""
    return float(x.strip())


def fortran_to_float(x: str) -> float:
    """Convert Fortran-style '1.0d-4' to Python float."""
    return float(x.lower().replace("d", "e"))


def tex_escape(text: str) -> str:
    """Escape special LaTeX characters in *text*."""
    out = text.replace("\\", r"\textbackslash{}")
    out = out.replace("_", r"\_")
    out = out.replace("%", r"\%")
    out = out.replace("&", r"\&")
    return out


def write_tex(path: Path, text: str) -> None:
    """Write a LaTeX snippet, creating parent dirs as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            out = {}
            for key in fieldnames:
                val = row.get(key, "")
                if isinstance(val, float):
                    if math.isnan(val):
                        out[key] = "nan"
                    else:
                        out[key] = f"{val:.5f}"
                else:
                    out[key] = str(val)
            w.writerow(out)
