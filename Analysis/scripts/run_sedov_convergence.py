#!/usr/bin/env python3
"""Run Sedov blast-wave convergence tests for SNEC / SuperSNEC.

Configurations (spatial convergence, fixed time):
  Static grids  (legacy_pattern):
    1. GridPattern.dat          N=100
    2. GridPattern1000.dat      N=1000
    3. GridPattern_adaptive.dat N=100
  Adaptive grids (adaptive_runtime, surface-concentrated, no remapping):
    4. N=100
    5. N=1000

Time-step convergence (fixed N=1000, legacy grid):
    6-8. dtmax = 0.1, 0.01, 0.001

All runs use the Sedov unit-test setup: gamma=1.4, E0=1 erg, rho0=1 g/cm3,
no radiation, no gravity, tend=4 s.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SNEC_ORIG = Path("/Volumes/Data4TB/Projects/SNEC/SNEC-1.01_original")
TMPBASE = ROOT / "tmp" / "sedov_convergence"
BASE_PARAMS = ROOT / "parameters"  # full SuperSNEC parameter file
BASE_PARAMS_ORIG = SNEC_ORIG / "parameters_sedov"  # original SNEC Sedov params


def read_base_params() -> str:
    """Read the working SuperSNEC parameters as a template."""
    return BASE_PARAMS.read_text()


def patch_param(text: str, key: str, value: str) -> str:
    """Replace or append a parameter value in the text."""
    pattern = re.compile(rf'^(\s*{re.escape(key)}\s*=\s*).*$', re.MULTILINE)
    if pattern.search(text):
        return pattern.sub(rf'\g<1>{value}', text)
    else:
        return text + f"\n{key} = {value}\n"


def sedov_params(
    imax: int,
    grid_mode: str = "legacy_pattern",
    grid_pattern_file: str = "tables/GridPattern.dat",
    grid_surface_alpha: float = 5.0,
    dtmax: float = 100.0,
) -> str:
    text = read_base_params()
    # Grid settings
    text = patch_param(text, "imax", str(imax))
    text = patch_param(text, "grid_mode", f'"{grid_mode}"')
    text = patch_param(text, "grid_pattern_file", f'"{grid_pattern_file}"')
    text = patch_param(text, "grid_surface_alpha", f"{grid_surface_alpha:.1f}d0")
    text = patch_param(text, "grid_update_interval_days", "0.0d0")
    # Sedov-specific: uniform sphere, no stellar physics
    text = patch_param(text, "profile_name", '"profiles/sedov.short"')
    text = patch_param(text, "comp_profile_name", '"profiles/sedov.iso.dat"')
    text = patch_param(text, "final_energy", "1.0d0")
    text = patch_param(text, "bomb_tstart", "0.0d0")
    text = patch_param(text, "bomb_tend", "1.0d-4")
    text = patch_param(text, "bomb_mass_spread", "1.0d-36")
    text = patch_param(text, "mass_excision", "0")
    text = patch_param(text, "radiation", "0")
    text = patch_param(text, "eoskey", "1")
    text = patch_param(text, "Ni_switch", "0")
    text = patch_param(text, "Ni_mass", "0.0d0")
    text = patch_param(text, "saha_ncomps", "0")
    text = patch_param(text, "boxcar_smoothing", "0")
    # Timing
    text = patch_param(text, "tend", "4.0d0")
    text = patch_param(text, "dtout", "1.0d-2")
    text = patch_param(text, "dtout_mid", "1.0d-2")
    text = patch_param(text, "dtout_fast", "1.0d-2")
    text = patch_param(text, "dtout_scalar", "1.0d-2")
    text = patch_param(text, "dtout_scalar_mid", "1.0d-2")
    text = patch_param(text, "dtout_scalar_fast", "1.0d-2")
    text = patch_param(text, "output_mid_transition_days", "1.0d-6")
    text = patch_param(text, "output_late_transition_days", "1.0d-5")
    text = patch_param(text, "dtmax", f"{dtmax}d0")
    text = patch_param(text, "ntinfo", "10000")
    text = patch_param(text, "sedov", "1")
    text = patch_param(text, "output_mode", "0")
    return text


def snec_orig_sedov_params(imax: int, gridding: str = "from_file_by_mass",
                           dtmax: float = 100.0) -> str:
    """Generate parameters for original SNEC Sedov test."""
    text = BASE_PARAMS_ORIG.read_text()
    text = patch_param(text, "imax", str(imax))
    text = patch_param(text, "gridding", f'"{gridding}"')
    text = patch_param(text, "dtmax", f"{dtmax}d0")
    return text


def setup_orig_snec_run(name: str, params_text: str,
                        grid_pattern_src: str) -> Path:
    """Set up a run dir for original SNEC with a specific grid pattern file.

    Original SNEC hardcodes the grid path to tables/GridPattern.dat, so we
    create a tables/ dir and symlink the desired pattern file as GridPattern.dat.
    """
    rundir = TMPBASE / name
    if rundir.exists():
        shutil.rmtree(rundir)
    rundir.mkdir(parents=True)
    (rundir / "Data").mkdir()

    (rundir / "parameters").write_text(params_text)

    # Symlink binary and profiles from original SNEC
    os.symlink(SNEC_ORIG / "snec", rundir / "snec")
    os.symlink(SNEC_ORIG / "profiles", rundir / "profiles")

    # Build tables dir: symlink everything from original, but override
    # GridPattern.dat with the desired source file
    tables_dst = rundir / "tables"
    tables_dst.mkdir()
    tables_src = SNEC_ORIG / "tables"
    for f in tables_src.iterdir():
        if f.name == "GridPattern.dat":
            continue
        os.symlink(f, tables_dst / f.name)
    # Symlink the desired grid pattern as GridPattern.dat.
    # For adaptive pattern: use SNEC_ORIG's own 100-line version (compatible),
    # since SuperSNEC's version has 1000 lines.
    grid_src = Path(grid_pattern_src)
    if grid_src.name == "GridPattern_adaptive.dat":
        os.symlink(SNEC_ORIG / "tables" / grid_src.name,
                    tables_dst / "GridPattern.dat")
    else:
        os.symlink(ROOT / grid_pattern_src, tables_dst / "GridPattern.dat")

    return rundir


# ---- Original SNEC spatial convergence runs ----
#  grid_pattern_src is the file in SuperSNEC/tables/ to use as GridPattern.dat
ORIG_SNEC_CONFIGS = [
    # (name, imax, gridding, grid_pattern_src)
    ("orig_static_100z",          100,  "from_file_by_mass", "tables/GridPattern.dat"),
    ("orig_static_1000z",         1000, "from_file_by_mass", "tables/GridPattern1000.dat"),
    ("orig_static_adaptive_100z", 100,  "from_file_by_mass", "tables/GridPattern_adaptive.dat"),
]


# ---- Spatial convergence runs ----
SPATIAL_CONFIGS = [
    # (name, imax, grid_mode, pattern_file, alpha)
    ("static_100z",          100,  "legacy_pattern", "tables/GridPattern.dat",          5.0),
    ("static_1000z",         1000, "legacy_pattern", "tables/GridPattern1000.dat",      5.0),
    ("static_adaptive_100z", 100,  "legacy_pattern", "tables/GridPattern_adaptive.dat", 5.0),
    ("adaptive_100z",        100,  "adaptive_runtime", "tables/GridPattern.dat",        7.0),
    ("adaptive_1000z",       1000, "adaptive_runtime", "tables/GridPattern.dat",        7.0),
    ("adaptive_lo_100z",     100,  "adaptive_runtime", "tables/GridPattern.dat",        2.0),
    ("adaptive_lo_1000z",    1000, "adaptive_runtime", "tables/GridPattern.dat",        2.0),
]

# ---- Time-step convergence runs ----
TIMESTEP_CONFIGS = [
    # (name, imax, dtmax)
    ("dt_0.1",   1000, 0.1),
    ("dt_0.01",  1000, 0.01),
    ("dt_0.001", 1000, 0.001),
]


def setup_run(name: str, params_text: str) -> Path:
    rundir = TMPBASE / name
    if rundir.exists():
        shutil.rmtree(rundir)
    rundir.mkdir(parents=True)
    (rundir / "Data").mkdir()

    (rundir / "parameters").write_text(params_text)

    # Symlink binary, profiles, and tables
    os.symlink(ROOT / "snec", rundir / "snec")
    os.symlink(ROOT / "profiles", rundir / "profiles")
    os.symlink(ROOT / "tables", rundir / "tables")

    return rundir


def run_snec(rundir: Path) -> float:
    t0 = time.time()
    result = subprocess.run(
        [str(rundir / "snec")],
        cwd=str(rundir),
        capture_output=True, text=True, timeout=600,
    )
    elapsed = time.time() - t0
    if result.returncode != 0:
        last = result.stderr[-500:] if result.stderr else result.stdout[-500:]
        print(f"  FAILED ({elapsed:.1f}s): {last}")
    else:
        lines = result.stdout.strip().split("\n")
        for line in lines[-3:]:
            print(f"  {line}")
    return elapsed


def main():
    print("=" * 60)
    print("Sedov Blast-Wave Convergence Tests")
    print("=" * 60)

    if not (ROOT / "snec").exists():
        print("Building snec...")
        subprocess.run(["make", "-C", str(ROOT)], check=True,
                        capture_output=True)

    # --- Original SNEC baseline ---
    if SNEC_ORIG.exists() and (SNEC_ORIG / "snec").exists():
        print("\n--- Original SNEC Baseline ---")
        for name, imax, gridding, gpat_src in ORIG_SNEC_CONFIGS:
            label = f"sedov_{name}"
            print(f"\n{label}  (N={imax}, {gridding}, {gpat_src})")
            params = snec_orig_sedov_params(imax, gridding=gridding)
            rundir = setup_orig_snec_run(label, params, gpat_src)
            elapsed = run_snec(rundir)
            print(f"  Wall time: {elapsed:.1f}s")
    else:
        print(f"\nWARNING: Original SNEC not found at {SNEC_ORIG}, skipping")

    # --- Spatial convergence ---
    print("\n--- Spatial Convergence ---")
    for name, imax, gmode, gpat, alpha in SPATIAL_CONFIGS:
        label = f"sedov_{name}"
        print(f"\n{label}  (N={imax}, {gmode})")
        params = sedov_params(imax, grid_mode=gmode,
                              grid_pattern_file=gpat,
                              grid_surface_alpha=alpha)
        rundir = setup_run(label, params)
        elapsed = run_snec(rundir)
        print(f"  Wall time: {elapsed:.1f}s")

    # --- Time-step convergence ---
    print("\n--- Time-step Convergence ---")
    for name, imax, dtmax in TIMESTEP_CONFIGS:
        label = f"sedov_{name}"
        print(f"\n{label}  (N={imax}, dtmax={dtmax})")
        params = sedov_params(imax, grid_mode="legacy_pattern",
                              grid_pattern_file="tables/GridPattern1000.dat",
                              dtmax=dtmax)
        rundir = setup_run(label, params)
        elapsed = run_snec(rundir)
        print(f"  Wall time: {elapsed:.1f}s")

    print("\n" + "=" * 60)
    print("All Sedov runs complete.")
    print(f"Output in: {TMPBASE}")


if __name__ == "__main__":
    main()
