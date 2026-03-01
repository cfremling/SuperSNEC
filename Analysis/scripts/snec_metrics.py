#!/usr/bin/env python3
"""Shared utilities for SNEC speed/accuracy sweep analysis."""

from __future__ import annotations

from pathlib import Path
import csv
import math
import re
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
    rms_gt30 = rms(d > 30.0)
    rms_gt5 = rms(d > 5.0)
    rms_all = rms(d >= 0.0)
    max_abs = float(np.max(np.abs(dmag)))
    n_pass = int(rms_0_5 <= 0.10) + int(rms_5_30 <= 0.10) + int(rms_gt30 <= 0.10)

    return {
        "rms_0_5": rms_0_5,
        "rms_5_30": rms_5_30,
        "rms_gt30": rms_gt30,
        "rms_gt5": rms_gt5,
        "rms_all": rms_all,
        "max_abs_mag": max_abs,
        "n_pass": float(n_pass),
    }


# ---------------------------------------------------------------------------
# Complete baseline parameter set — every sweep MUST start from this
# ---------------------------------------------------------------------------

BASELINE_PARAMETERS: Dict[str, str] = {
    # --- Profile and explosion ---
    "profile_name": '"profiles/stripped_star.short"',
    "comp_profile_name": '"profiles/stripped_star.iso.dat"',
    "initial_data": '"Thermal_Bomb"',
    "piston_vel": "5.0d9",
    "piston_tstart": "0.0d0",
    "piston_tend": "1.0d-2",
    "final_energy": "1.0d51",
    "bomb_tstart": "0.0d0",
    "bomb_tend": "0.1d0",
    "bomb_mass_spread": "0.1d0",
    "bomb_start_point": "1",
    # --- Grid ---
    "imax": "100",
    "grid_mode": '"adaptive_runtime"',
    "grid_update_interval_days": "1.0d0",
    "grid_pattern_file": '"tables/GridPattern.dat"',
    "grid_surface_alpha": "7.0",
    "grid_relax_days": "5.0",
    "grid_min_cell_frac": "1.0E-4",
    "grid_debug": "0",
    "grid_adaptive_interval": "0",
    # --- Excision ---
    "mass_excision": "1",
    "mass_excised": "1.4",
    # --- Physics ---
    "radiation": "1",
    "eoskey": "2",
    # --- Nickel ---
    "Ni_switch": "1",
    "Ni_mass": "0.03",
    "Ni_mix_fraction": "0.31d0",
    "Ni_mix_component2_fraction": "0.0d0",
    "Ni_mix_component2_extent": "0.31d0",
    "Ni_mix_kernel": "1",
    "Ni_period": "5.0E4",
    "Ni_period_max": "5.456E5",
    "Ni_fractional_change": "0.70d0",
    "Ni_quad_npoints": "70",
    # --- EOS / Ni smoothing ---
    "saha_ncomps": "3",
    "smooth_ni_luminosity": "0",
    # --- Boxcar ---
    "boxcar_smoothing": "1",
    "boxcar_smooth_ni": "1",
    "boxcar_nominal_mass_msun": "0.4d0",
    "boxcar_number_iterations": "4",
    "boxcar_ni_nominal_mass_msun": "0.4d0",
    "boxcar_ni_number_iterations": "4",
    # --- Opacity ---
    "opacity_floor_envelope": "0.01d0",
    "opacity_floor_core": "0.24d0",
    # --- Time limits ---
    "ntmax": "10000000000000",
    "tend": "5.0E6",
    # --- Output cadence (three-band) ---
    "dtout": "5.0d5",
    "dtout_mid": "1.0d5",
    "dtout_fast": "5.0d4",
    "dtout_scalar": "1.0d5",
    "dtout_scalar_mid": "5.0d4",
    "dtout_scalar_fast": "1.0d4",
    "output_mid_transition_days": "10.0d0",
    "output_late_transition_days": "2.0d1",
    "dtout_check": "1.7d9",
    "ntout": "-1",
    "ntout_scalar": "-1",
    "ntout_check": "-1",
    "ntinfo": "1000",
    "dtmin": "1.0d-10",
    "dtmax": "1.0E5",
    # --- Solver ---
    "epstol_hydro": "1.0d-4",
    "epstol_rad": "1.0d-4",
    "itmax_hydro": "100",
    "itmax_rad": "300",
    "ni_raytrace_opt": "1",
    "output_mode": "0",
    # --- Sedov test ---
    "sedov": "0",
    # --- Output directory (always overridden per sweep) ---
    "outdir": '"Data"',
}


def write_parameters_file(
    path: Path,
    overrides: Dict[str, str] | None = None,
) -> None:
    """Write a complete, self-contained parameters file from the baseline template.

    Every parameter is written explicitly — nothing depends on the user's
    parameters file.
    """
    params = dict(BASELINE_PARAMETERS)
    if overrides:
        params.update(overrides)
    lines = [f"{key} = {value}" for key, value in params.items()]
    path.write_text("\n".join(lines) + "\n")


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
