#!/usr/bin/env python3
"""Run energy conservation comparison: original SNEC vs SuperSNEC.

All runs use the stripped-star progenitor with SNEC-default physics settings.
SuperSNEC runs use legacy_pattern grid mode and SNEC-equivalent Ni settings
so the ONLY differences are the code changes, not parameter choices.

Configurations:
  1. Original SNEC, N=100   (EPSTOL=1e-7, nquad=150, Ni_period=1e4)
  2. Original SNEC, N=1000  (same)
  3. SuperSNEC SNEC-compat, N=100  (legacy_pattern, EPSTOL=1e-7, nquad=150, f_Ni conservative)
  4. SuperSNEC SNEC-compat, N=1000 (same)
  5. SuperSNEC baseline, N=100     (adaptive_runtime, EPSTOL=1e-4, nquad=70, f_Ni=0.20)
  6. SuperSNEC baseline, N=1000    (adaptive_runtime, EPSTOL=1e-4, nquad=150)
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SNEC_ORIG = Path("/Volumes/Data4TB/Projects/SNEC/SNEC-1.01_original")
SUPERSNEC = ROOT
TMPBASE = ROOT / "tmp" / "conservation_comparison"

# Stripped star profile (relative paths from run dir)
PROFILE = "profiles/stripped_star.short"
COMP_PROFILE = "profiles/stripped_star.iso.dat"


def write_snec_original_params(path: Path, imax: int) -> None:
    """Write parameter file for original SNEC."""
    # tend=5e6 s ~ 57.9 d to match SuperSNEC runs
    # dtout=1e4 s (~0.12 d) for dense .xg output
    path.write_text(f"""\
outdir              = "Data"
profile_name        = "{PROFILE}"
comp_profile_name   = "{COMP_PROFILE}"
initial_data        = "Thermal_Bomb"
final_energy        = 1.0d51
bomb_tstart         = 0.0d0
bomb_tend           = 0.1d0
bomb_mass_spread    = 0.1d0
bomb_start_point    = 1
imax                = {imax}
gridding            = "from_file_by_mass"
mass_excision       = 1
mass_excised        = 1.4
radiation           = 1
eoskey              = 2
Ni_switch           = 1
Ni_mass             = 0.03
Ni_boundary_mass    = 2.0d0
Ni_period           = 1.0d4
saha_ncomps         = 3
boxcar_smoothing    = 1
opacity_floor_envelope = 0.01d0
opacity_floor_core     = 0.24d0
ntmax               = 10000000000000
tend                = 5.0d6
dtout               = 1.0d4
dtout_scalar        = 1.0d4
dtout_check         = 1.0d4
ntout               = -1
ntout_scalar        = -1
ntout_check         = -1
ntinfo              = 1000
dtmin               = 1.0d-10
dtmax               = 1.0d4
sedov               = 0
""")


def write_supersnec_compat_params(path: Path, imax: int) -> None:
    """Write SuperSNEC params matching SNEC defaults exactly."""
    path.write_text(f"""\
outdir              = "Data"
profile_name        = "{PROFILE}"
comp_profile_name   = "{COMP_PROFILE}"
initial_data        = "Thermal_Bomb"
final_energy        = 1.0d51
bomb_tstart         = 0.0d0
bomb_tend           = 0.1d0
bomb_mass_spread    = 0.1d0
bomb_start_point    = 1
imax                = {imax}
grid_mode           = "legacy_pattern"
grid_pattern_file   = "tables/GridPattern.dat"
grid_update_interval_days = 0.0
mass_excision       = 1
mass_excised        = 1.4
radiation           = 1
eoskey              = 2
Ni_switch           = 1
Ni_mass             = 0.03
Ni_by_hand          = 1
Ni_mix_fraction     = 0.31
Ni_mix_kernel       = 1
Ni_period           = 1.0d4
Ni_period_max       = 1.0d4
Ni_fractional_change = 0.01
Ni_quad_npoints     = 150
saha_ncomps         = 3
boxcar_smoothing    = 1
smooth_ni_luminosity = 0
boxcar_smooth_ni    = 1
opacity_floor_envelope = 0.01d0
opacity_floor_core     = 0.24d0
ntmax               = 10000000000000
tend                = 5.0d6
dtout               = 1.0d4
dtout_scalar        = 1.0d4
dtout_check         = 1.0d4
ntout               = -1
ntout_scalar        = -1
ntout_check         = -1
ntinfo              = 1000
dtmin               = 1.0d-10
dtmax               = 1.0d4
sedov               = 0
epstol_hydro        = 1.0d-7
epstol_rad          = 1.0d-7
itmax_hydro         = 100
itmax_rad           = 300
ni_raytrace_opt     = 0
ni_ray_interp       = 0
output_mode         = 0
grid_remap_qphoto_stop = 0.50d0
dtout_mid           = 1.0d4
dtout_fast          = 1.0d4
dtout_scalar_mid    = 1.0d4
dtout_scalar_fast   = 1.0d4
output_mid_transition_days = 10.0d0
output_late_transition_days = 20.0d0
""")


def write_supersnec_baseline_params(path: Path, imax: int) -> None:
    """Write SuperSNEC adaptive_runtime with SNEC-identical solver/quad/Ni settings.

    The ONLY difference from compat is the grid mode: adaptive_runtime instead of
    legacy_pattern. Everything else (EPSTOL, nquad, Ni cadence, dtmax) matches
    the original SNEC defaults so we isolate the effect of adaptive gridding.
    """
    path.write_text(f"""\
outdir              = "Data"
profile_name        = "{PROFILE}"
comp_profile_name   = "{COMP_PROFILE}"
initial_data        = "Thermal_Bomb"
final_energy        = 1.0d51
bomb_tstart         = 0.0d0
bomb_tend           = 0.1d0
bomb_mass_spread    = 0.1d0
bomb_start_point    = 1
imax                = {imax}
grid_mode           = "adaptive_runtime"
grid_surface_alpha  = 7.0d0
grid_relax_days     = 5.0d0
grid_update_interval_days = 1.0d0
mass_excision       = 1
mass_excised        = 1.4
radiation           = 1
eoskey              = 2
Ni_switch           = 1
Ni_mass             = 0.03
Ni_by_hand          = 1
Ni_mix_fraction     = 0.31
Ni_mix_kernel       = 1
Ni_period           = 1.0d4
Ni_period_max       = 1.0d4
Ni_fractional_change = 0.01
Ni_quad_npoints     = 150
saha_ncomps         = 3
boxcar_smoothing    = 1
smooth_ni_luminosity = 0
boxcar_smooth_ni    = 1
opacity_floor_envelope = 0.01d0
opacity_floor_core     = 0.24d0
ntmax               = 10000000000000
tend                = 5.0d6
dtout               = 1.0d4
dtout_scalar        = 1.0d4
dtout_scalar_fast   = 1.0d4
dtout_scalar_mid    = 1.0d4
dtout_fast          = 1.0d4
dtout_mid           = 1.0d4
output_mid_transition_days = 10.0d0
output_late_transition_days = 20.0d0
dtout_check         = 1.0d4
ntout               = -1
ntout_scalar        = -1
ntout_check         = -1
ntinfo              = 1000
dtmin               = 1.0d-10
dtmax               = 1.0d4
sedov               = 0
epstol_hydro        = 1.0d-7
epstol_rad          = 1.0d-7
itmax_hydro         = 100
itmax_rad           = 300
ni_raytrace_opt     = 0
ni_ray_interp       = 0
output_mode         = 0
grid_remap_qphoto_stop = 0.50d0
""")


CONFIGS = [
    ("snec_orig_100z",       SNEC_ORIG,  write_snec_original_params,     100),
    ("snec_orig_1000z",      SNEC_ORIG,  write_snec_original_params,     1000),
    ("supersnec_compat_100z",  SUPERSNEC, write_supersnec_compat_params,  100),
    ("supersnec_compat_1000z", SUPERSNEC, write_supersnec_compat_params,  1000),
    ("supersnec_baseline_100z",  SUPERSNEC, write_supersnec_baseline_params, 100),
    ("supersnec_baseline_1000z", SUPERSNEC, write_supersnec_baseline_params, 1000),
]


def setup_run(name: str, code_root: Path, param_writer, imax: int) -> Path:
    """Set up a run directory with symlinks and parameters."""
    rundir = TMPBASE / name
    if rundir.exists():
        shutil.rmtree(rundir)
    rundir.mkdir(parents=True)
    (rundir / "Data").mkdir()

    # Write parameters
    param_writer(rundir / "parameters", imax)

    # Symlink code binary and profiles
    os.symlink(code_root / "snec", rundir / "snec")
    os.symlink(code_root / "profiles", rundir / "profiles")

    # Create tables dir: symlink all files, but provide the right GridPattern.dat
    tables_src = code_root / "tables"
    tables_dst = rundir / "tables"
    tables_dst.mkdir()
    for f in tables_src.iterdir():
        if f.name == "GridPattern.dat":
            continue  # we'll provide our own
        os.symlink(f, tables_dst / f.name)

    # For original SNEC: GridPattern.dat must match imax
    # SuperSNEC: uses grid_pattern_file parameter, but also needs GridPattern.dat present
    gridpat_src = SUPERSNEC / "tables"
    if imax <= 100:
        shutil.copy2(gridpat_src / "GridPattern.dat", tables_dst / "GridPattern.dat")
    else:
        shutil.copy2(gridpat_src / "GridPattern1000.dat", tables_dst / "GridPattern.dat")
    # Also copy GridPattern1000.dat if SuperSNEC needs it
    if not (tables_dst / "GridPattern1000.dat").exists():
        shutil.copy2(gridpat_src / "GridPattern1000.dat", tables_dst / "GridPattern1000.dat")

    return rundir


def run_snec(rundir: Path) -> float:
    """Run SNEC and return wall-clock time."""
    t0 = time.time()
    result = subprocess.run(
        [str(rundir / "snec")],
        cwd=str(rundir),
        capture_output=True, text=True, timeout=1800
    )
    elapsed = time.time() - t0
    if result.returncode != 0:
        print(f"  FAILED: {result.stderr[-300:]}")
    else:
        # Print last few lines
        lines = result.stdout.strip().split("\n")
        for line in lines[-4:]:
            print(f"  {line}")
    return elapsed


def load_conservation(rundir: Path):
    """Load conservation.dat and return key metrics."""
    path = rundir / "Data" / "conservation.dat"
    if not path.exists():
        return None
    data = np.loadtxt(path)
    return data


def main():
    print("=" * 72)
    print("Energy Conservation Comparison: Original SNEC vs SuperSNEC")
    print("=" * 72)

    results = {}

    for name, code_root, param_writer, imax in CONFIGS:
        print(f"\n--- {name} (N={imax}) ---")
        rundir = setup_run(name, code_root, param_writer, imax)
        elapsed = run_snec(rundir)
        print(f"  Wall time: {elapsed:.1f} s")

        data = load_conservation(rundir)
        if data is not None and data.ndim == 2 and data.shape[1] >= 11:
            time_d = data[:, 0] / 86400
            frac = data[:, 10]
            results[name] = (time_d, frac, elapsed)

            print(f"  Energy balance residual:")
            for t_target in [0.1, 1, 10, 30, 57]:
                idx = np.argmin(np.abs(time_d - t_target))
                print(f"    t={t_target:5.1f}d: eps_E = {frac[idx]:+.4%}")
        else:
            print("  No conservation data found or wrong format")
            results[name] = None

    # Summary table
    print("\n" + "=" * 72)
    print("SUMMARY TABLE")
    print("=" * 72)
    header = f"{'Configuration':<30s} {'Runtime':>8s}"
    for t in [0.1, 10, 30, 57]:
        header += f"  {'t='+str(t)+'d':>10s}"
    print(header)
    print("-" * len(header))
    for name, code_root, param_writer, imax in CONFIGS:
        r = results.get(name)
        if r is None:
            print(f"{name:<30s} {'FAILED':>8s}")
            continue
        time_d, frac, elapsed = r
        row = f"{name:<30s} {elapsed:7.1f}s"
        for t in [0.1, 10, 30, 57]:
            idx = np.argmin(np.abs(time_d - t))
            row += f"  {frac[idx]:+9.4%}"
        print(row)


if __name__ == "__main__":
    main()
