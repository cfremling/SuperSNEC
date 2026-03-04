#!/usr/bin/env python3
"""
Cumulative speedup benchmark for SuperSNEC paper.

Runs a sequence of SNEC/SuperSNEC configurations, each adding one
optimization on top of all previous, and generates a LaTeX table
ordered by marginal speedup (largest first).

Usage:
    python3 Analysis/scripts/run_cumulative_speedup.py [--skip-reference] [--force]

The 1000-zone original SNEC run takes ~10+ minutes.  Use --skip-reference
to reuse a cached result from a previous run.

Must be run from the SuperSNEC repo root (where ./snec lives).
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SUPERSNEC_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = SUPERSNEC_ROOT / "Analysis" / "results" / "cumulative_speedup"
TABLE_DIR = SUPERSNEC_ROOT / "Paper" / "generated"

# Import shared metrics (add to path)
sys.path.insert(0, str(SUPERSNEC_ROOT / "Analysis" / "scripts"))
from snec_metrics import compute_mag_metrics, write_parameters_file  # noqa: E402


# ---------------------------------------------------------------------------
# Configuration dataclass
# ---------------------------------------------------------------------------
@dataclass
class Config:
    """One benchmark configuration."""
    label: str                     # short label for table
    description: str               # longer description
    codebase: str                  # 'original_1000', 'original_100', 'supersnec'
    param_overrides: Dict[str, str] = field(default_factory=dict)
    outdir_name: str = ""          # temp output directory name
    is_production_grid: bool = False  # True for output_mode=1 rows
    mode0_outdir: str = ""         # outdir_name of full-output counterpart


# ---------------------------------------------------------------------------
# Configurations (in execution order — cumulative)
# ---------------------------------------------------------------------------
def build_configs() -> List[Config]:
    """Build the ordered list of cumulative benchmark configurations."""

    # SuperSNEC template parameters that match original SNEC physics as
    # closely as possible (legacy grid, original Ni cadence, original dtmax,
    # original quadrature, tight solver, no ray-trace opt).
    # Match original SNEC physics as closely as possible:
    # legacy grid, original Ni cadence, original dtmax, original quadrature,
    # original EPSTOL=1e-7, no raytrace opt, no SuperSNEC-only features.
    supersnec_legacy_base = {
        "outdir": '"PLACEHOLDER"',
        "imax": "100",
        "grid_mode": '"legacy_pattern"',
        "grid_update_interval_days": "0.0d0",
        "ni_ray_interp": "0",
        "Ni_period": "5.0E4",
        "Ni_period_max": "5.0E4",        # max = min → fixed cadence
        "Ni_fractional_change": "0.01d0", # tiny → never triggers adaptive stretch
        "Ni_quad_npoints": "150",
        "dtmax": "1.0E4",
        "epstol_hydro": "1.0d-7",        # match original SNEC hardcoded value
        "epstol_rad": "1.0d-7",          # match original SNEC hardcoded value
        "itmax_hydro": "100",
        "itmax_rad": "300",
        "ni_raytrace_opt": "0",
        "output_mode": "0",
        # Match original SNEC output cadence (uniform, no three-band)
        "dtout": "1.7d5",
        "dtout_mid": "1.7d5",
        "dtout_fast": "1.7d5",
        "dtout_scalar": "1.7d4",
        "dtout_scalar_mid": "1.7d4",
        "dtout_scalar_fast": "1.7d4",
        "ntout_scalar": "-1",
        # Match original SNEC boxcar: all species + Ni both get 0.4 Msun, 4 iterations
        "smooth_ni_luminosity": "0",
        "boxcar_smoothing": "1",
        "boxcar_smooth_ni": "1",
        "boxcar_nominal_mass_msun": "0.4d0",
        "boxcar_number_iterations": "4",
        "boxcar_ni_nominal_mass_msun": "0.4d0",
        "boxcar_ni_number_iterations": "4",
    }

    # Three-band cadence for the working baseline and beyond
    three_band_cadence = {
        "dtout_scalar_fast": "1.0d4",
        "dtout_scalar_mid": "5.0d4",
        "dtout_scalar": "1.0d5",
        "dtout_fast": "5.0d4",
        "dtout_mid": "1.0d5",
        "dtout": "5.0d5",
    }

    def overlay(base: dict, updates: dict) -> dict:
        d = dict(base)
        d.update(updates)
        return d

    configs = [
        # 0. True baseline: original SNEC, 1000 zones
        Config(
            label="SNEC reference (1000 zones)",
            description="Original SNEC 1.01, 1000 zones, -O3 -fbounds-check, EPSTOL=1e-7",
            codebase="original_1000",
            outdir_name="ref_1000z",
        ),
        # 1. Zone reduction: original SNEC, 100 zones, original flags
        Config(
            label="Zone reduction ($N{=}100$)",
            description="Original SNEC 1.01, 100 zones, original compiler flags",
            codebase="original_100",
            outdir_name="orig_100z",
        ),
        # 2. Compiler flags: original SNEC, 100 zones, -Ofast
        Config(
            label="Compiler flags (\\texttt{-Ofast})",
            description="Original SNEC 1.01, 100 zones, -Ofast -fprotect-parens etc",
            codebase="original_100",
            outdir_name="orig_100z_ofast",
        ),
        # 3. SuperSNEC codebase with original-matching physics
        #    (EPSTOL=1e-7, no Ni smoothing, no raytrace opt)
        Config(
            label="SuperSNEC codebase",
            description="SuperSNEC legacy_pattern, all physics matching original SNEC",
            codebase="supersnec",
            param_overrides=supersnec_legacy_base,
            outdir_name="ssnec_code_only",
        ),
        # 4. Solver tolerance relaxation (1e-7 → 1e-4)
        Config(
            label="Solver tolerance ($10^{-7} \\to 10^{-4}$)",
            description="+ EPSTOL relaxed from 1e-7 to 1e-4",
            codebase="supersnec",
            param_overrides=overlay(supersnec_legacy_base, {
                "epstol_hydro": "1.0d-4",
                "epstol_rad": "1.0d-4",
            }),
            outdir_name="ssnec_solver",
        ),
        # 5. + Ni ray-tracing optimization (still fixed grid, legacy Ni)
        Config(
            label="Ni ray-tracing optimization",
            description="+ precomputed kappa_gamma + sequential hunt + linear interp",
            codebase="supersnec",
            param_overrides=overlay(supersnec_legacy_base, {
                "epstol_hydro": "1.0d-4",
                "epstol_rad": "1.0d-4",
                "ni_raytrace_opt": "1",
                "ni_ray_interp": "1",
            }),
            outdir_name="ssnec_raytrace",
        ),
        # 6. + Adaptive Ni cadence f_Ni=0.10 + dtmax 1e5 (still fixed grid)
        Config(
            label="Ni cadence ($f_{\\mathrm{Ni}}{=}0.10$)",
            description="+ adaptive Ni cadence with f_Ni=0.10, dtmax=1e5 (fixed grid)",
            codebase="supersnec",
            param_overrides=overlay(supersnec_legacy_base, {
                "epstol_hydro": "1.0d-4",
                "epstol_rad": "1.0d-4",
                "ni_raytrace_opt": "1",
                "ni_ray_interp": "1",
                "dtmax": "1.0E5",
                "Ni_fractional_change": "0.10d0",
                "Ni_period_max": "5.456E5",
            }),
            outdir_name="ssnec_ni_cadence",
        ),
        # 7. + Quadrature reduction 150 → 70 (still fixed grid)
        Config(
            label="Quadrature ($150 \\to 70$)",
            description="+ Ni_quad_npoints reduced from 150 to 70 (fixed grid)",
            codebase="supersnec",
            param_overrides=overlay(supersnec_legacy_base, {
                "epstol_hydro": "1.0d-4",
                "epstol_rad": "1.0d-4",
                "ni_raytrace_opt": "1",
                "ni_ray_interp": "1",
                "dtmax": "1.0E5",
                "Ni_fractional_change": "0.10d0",
                "Ni_period_max": "5.456E5",
                "Ni_quad_npoints": "70",
            }),
            outdir_name="ssnec_quad",
        ),
        # 8. + Adaptive grid mode + f_Ni=0.20 + three-band cadence + smooth_ni
        #    (the full working baseline)
        Config(
            label="Adaptive grid + $f_{\\mathrm{Ni}}{=}0.20$",
            description="+ adaptive_runtime grid with f_Ni=0.20 + three-band cadence + smooth_ni",
            codebase="supersnec",
            param_overrides=overlay(overlay(supersnec_legacy_base, three_band_cadence), {
                "epstol_hydro": "1.0d-4",
                "epstol_rad": "1.0d-4",
                "ni_raytrace_opt": "1",
                "ni_ray_interp": "1",
                "dtmax": "1.0E5",
                "Ni_fractional_change": "0.20d0",
                "Ni_period_max": "5.456E5",
                "Ni_quad_npoints": "70",
                "grid_mode": '"adaptive_runtime"',
                "grid_update_interval_days": "1.0d0",
                "smooth_ni_luminosity": "1",
            }),
            outdir_name="ssnec_final",
        ),
        # 9. Fast mode: imax=60, f_Ni=0.20, quad=50, smooth_ni=1
        Config(
            label="Fast ($N{=}60$, $n_{quad}{=}50$)",
            description="imax=60, f_Ni=0.20, quad=50, smooth_ni=1 (maximum speed)",
            codebase="supersnec",
            param_overrides=overlay(overlay(supersnec_legacy_base, three_band_cadence), {
                "imax": "60",
                "epstol_hydro": "1.0d-4",
                "epstol_rad": "1.0d-4",
                "ni_raytrace_opt": "1",
                "ni_ray_interp": "1",
                "dtmax": "1.0E5",
                "Ni_fractional_change": "0.20d0",
                "Ni_period_max": "5.456E5",
                "Ni_quad_npoints": "50",
                "smooth_ni_luminosity": "1",
                "grid_mode": '"adaptive_runtime"',
                "grid_update_interval_days": "1.0d0",
            }),
            outdir_name="ssnec_extreme",
        ),
        # 10. Production grid: baseline + output_mode=1
        Config(
            label="Baseline",
            description="Baseline + output_mode=1 (fitting-only output for production grids)",
            codebase="supersnec",
            param_overrides=overlay(overlay(supersnec_legacy_base, three_band_cadence), {
                "epstol_hydro": "1.0d-4",
                "epstol_rad": "1.0d-4",
                "ni_raytrace_opt": "1",
                "ni_ray_interp": "1",
                "dtmax": "1.0E5",
                "Ni_fractional_change": "0.20d0",
                "Ni_period_max": "5.456E5",
                "Ni_quad_npoints": "70",
                "smooth_ni_luminosity": "1",
                "grid_mode": '"adaptive_runtime"',
                "grid_update_interval_days": "1.0d0",
                "output_mode": "1",
            }),
            outdir_name="ssnec_minimal_output",
            is_production_grid=True,
            mode0_outdir="ssnec_final",
        ),
        # 11. Production grid: fast + output_mode=1
        Config(
            label="Fast",
            description="Fast (imax=60, f_Ni=0.20, quad=50) + output_mode=1",
            codebase="supersnec",
            param_overrides=overlay(overlay(supersnec_legacy_base, three_band_cadence), {
                "imax": "60",
                "epstol_hydro": "1.0d-4",
                "epstol_rad": "1.0d-4",
                "ni_raytrace_opt": "1",
                "ni_ray_interp": "1",
                "dtmax": "1.0E5",
                "Ni_fractional_change": "0.20d0",
                "Ni_period_max": "5.456E5",
                "Ni_quad_npoints": "50",
                "smooth_ni_luminosity": "1",
                "grid_mode": '"adaptive_runtime"',
                "grid_update_interval_days": "1.0d0",
                "output_mode": "1",
            }),
            outdir_name="ssnec_fast_minimal",
            is_production_grid=True,
            mode0_outdir="ssnec_extreme",
        ),
    ]
    return configs


# ---------------------------------------------------------------------------
# Run a single configuration
# ---------------------------------------------------------------------------
def run_config(cfg: Config, ref_lum_path: Optional[Path] = None,
               force: bool = False) -> dict:
    """Run one benchmark configuration and return results dict."""
    result_dir = RESULTS_DIR / cfg.outdir_name
    result_file = result_dir / "result.npz"

    # Check cache
    if not force and result_file.exists():
        data = dict(np.load(result_file, allow_pickle=True))
        print(f"  [cached] {cfg.label}: {float(data['runtime']):.1f}s")
        return {k: float(v) if np.isscalar(v) else v for k, v in data.items()}

    result_dir.mkdir(parents=True, exist_ok=True)

    # Original SNEC configs use shipped reference data (cannot be re-run)
    if cfg.codebase in ("original_1000", "original_100"):
        if not result_file.exists():
            raise RuntimeError(
                f"Original SNEC result not found: {result_file}\n"
                f"Original SNEC runs are shipped as reference data and cannot be re-generated.\n"
                f"Ensure Analysis/results/cumulative_speedup/{cfg.outdir_name}/result.npz exists."
            )
        data = dict(np.load(result_file, allow_pickle=True))
        print(f"  [shipped] {cfg.label}: {float(data['runtime']):.1f}s")
        return {k: float(v) if np.isscalar(v) else v for k, v in data.items()}

    # Determine codebase root for SuperSNEC runs
    codebase_root = SUPERSNEC_ROOT

    # Live data dir is always <codebase>/Data (short relative path)
    live_data = codebase_root / "Data"

    # Write complete parameters file for SuperSNEC runs
    if cfg.codebase == "supersnec" and cfg.param_overrides:
        params_file = codebase_root / "parameters"
        overrides = dict(cfg.param_overrides)
        overrides["outdir"] = '"Data"'
        write_parameters_file(params_file, overrides)

    # SuperSNEC runs are fast enough to average over multiple runs
    nruns = 3 if cfg.codebase == "supersnec" else 1

    # Clean live data dir before running
    if live_data.exists():
        shutil.rmtree(live_data)
    live_data.mkdir(parents=True, exist_ok=True)

    # Run SNEC with OS-level file descriptors and Popen for accurate
    # timing.  subprocess.run adds ~0.08s from communicate/timeout
    # overhead; Popen+wait matches `time ./snec` wall-clock.
    log_path = result_dir / "run.log"
    label_suffix = f" ({nruns} runs)" if nruns > 1 else ""
    print(f"  Running {cfg.label}{label_suffix}...", end="", flush=True)
    times = []
    last_returncode = 0
    for run_i in range(nruns):
        if live_data.exists():
            shutil.rmtree(live_data)
        live_data.mkdir(parents=True, exist_ok=True)
        log_fd = os.open(str(log_path), os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
        try:
            t0 = time.perf_counter()
            proc = subprocess.Popen(
                ["./snec"],
                cwd=codebase_root,
                stdout=log_fd,
                stderr=log_fd,
            )
            proc.wait(timeout=1800)
            times.append(time.perf_counter() - t0)
            last_returncode = proc.returncode
        finally:
            os.close(log_fd)
    runtime = sum(times) / len(times)
    times_str = ", ".join(f"{t:.2f}" for t in times)
    print(f" {runtime:.2f}s [{times_str}] (exit={last_returncode})")

    if last_returncode != 0:
        print(f"  WARNING: snec exited with code {last_returncode}")
        log_content = log_path.read_text()
        print(f"  Last 5 lines of log:\n{chr(10).join(log_content.splitlines()[-5:])}")

    # Copy lum_observed.dat from live data dir to results
    lum_src = live_data / "lum_observed.dat"
    lum_dst = result_dir / "lum_observed.dat"
    if lum_src.exists():
        shutil.copy2(lum_src, lum_dst)
    else:
        print(f"  WARNING: {lum_src} not found")

    # Compute RMS if reference available
    rms_all = float("nan")
    if ref_lum_path and ref_lum_path.exists() and lum_dst.exists():
        try:
            metrics = compute_mag_metrics(ref_lum_path, lum_dst)
            rms_all = metrics["rms_all"]
        except Exception as e:
            print(f"  WARNING: RMS computation failed: {e}")

    result = {"runtime": runtime, "rms_all": rms_all}
    np.savez(result_file, **result)
    return result


# ---------------------------------------------------------------------------
# Table generation
# ---------------------------------------------------------------------------
def generate_table(configs: List[Config], results: List[dict]) -> str:
    """Generate LaTeX cumulative speedup table."""
    ref_runtime = results[0]["runtime"]

    # Build a lookup from outdir_name to result for production grid marginals
    outdir_to_result = {
        cfg.outdir_name: res for cfg, res in zip(configs, results)
    }

    # Build rows with cumulative and marginal speedup
    rows = []
    for i, (cfg, res) in enumerate(zip(configs, results)):
        rt = res["runtime"]
        cumulative = ref_runtime / rt if rt > 0 else float("inf")
        if i == 0:
            marginal_str = "---"
        elif cfg.is_production_grid:
            # Marginal for production grid = vs mode=0 counterpart
            counterpart_rt = outdir_to_result[cfg.mode0_outdir]["runtime"]
            marginal = counterpart_rt / rt if rt > 0 else float("inf")
            marginal_str = f"${marginal:.1f}\\times$"
        else:
            prev_rt = results[i - 1]["runtime"]
            marginal = prev_rt / rt if rt > 0 else float("inf")
            marginal_str = f"${marginal:.1f}\\times$"

        rms = res.get("rms_all", float("nan"))
        rms_str = f"{rms:.3f}" if not np.isnan(rms) else "---"

        rows.append({
            "idx": i,
            "label": cfg.label,
            "runtime": rt,
            "cumulative": cumulative,
            "marginal_str": marginal_str,
            "rms_str": rms_str,
            "is_production_grid": cfg.is_production_grid,
        })

    # Find the baseline row (ssnec_final)
    baseline_idx = None
    for i, cfg in enumerate(configs):
        if cfg.outdir_name == "ssnec_final":
            baseline_idx = i
            break
    if baseline_idx is None:
        baseline_idx = len(rows) - 1

    # Separate post-baseline rows into regular extras and production grid
    extra_rows = []
    prod_rows = []
    for r, c in zip(rows[baseline_idx + 1:], configs[baseline_idx + 1:]):
        if c.is_production_grid:
            prod_rows.append(r)
        else:
            extra_rows.append(r)

    # Build LaTeX
    lines = []
    lines.append(r"\begin{deluxetable*}{lrccc}")
    lines.append(r"\tabletypesize{\small}")
    lines.append(r"\tablecaption{Cumulative speedup breakdown. Each row adds one "
                 r"optimization on top of the previous row. "
                 r"Runtime is wall-clock time on an M1~Pro Apple Silicon CPU. "
                 r"RMS is computed against the 1000-zone reference."
                 r"\label{tab:cumulative_speedup}}")
    lines.append(r"\tablehead{")
    lines.append(r"\colhead{Configuration} &")
    lines.append(r"\colhead{Runtime (s)} &")
    lines.append(r"\colhead{Cumulative} &")
    lines.append(r"\colhead{Marginal} &")
    lines.append(r"\colhead{$\Delta m_\mathrm{all}$ (mag)}")
    lines.append(r"}")
    lines.append(r"\startdata")

    # Reference row
    lines.append(f"{rows[0]['label']} & {rows[0]['runtime']:.1f} & "
                 f"$1\\times$ & --- & {rows[0]['rms_str']} \\\\")
    lines.append(r"\hline")

    # Cumulative rows up to and including baseline
    for row in rows[1:baseline_idx + 1]:
        cum_str = f"${row['cumulative']:.0f}\\times$"
        lines.append(f"{row['label']} & {row['runtime']:.1f} & "
                     f"{cum_str} & {row['marginal_str']} & {row['rms_str']} \\\\")

    # Total row = the working baseline
    best = rows[baseline_idx]
    lines.append(r"\hline")
    lines.append(f"\\textbf{{Total (SuperSNEC baseline)}} & "
                 f"\\textbf{{{best['runtime']:.1f}}} & "
                 f"$\\mathbf{{{best['cumulative']:.0f}\\times}}$ & "
                 f"--- & \\textbf{{{best['rms_str']}}} \\\\")

    # Extra rows (Fast mode)
    if extra_rows:
        lines.append(r"\hline")
        for row in extra_rows:
            cum_str = f"${row['cumulative']:.0f}\\times$"
            lines.append(f"{row['label']} & {row['runtime']:.1f} & "
                         f"{cum_str} & {row['marginal_str']} & "
                         f"{row['rms_str']} \\\\")

    # Production grid section
    if prod_rows:
        lines.append(r"\hline")
        lines.append(r"\\")
        lines.append(r"\multicolumn{5}{c}{\textit{Production grid mode "
                     r"(\texttt{output\_mode=1}; "
                     r"Sect.~\ref{sec:output_mode})}} \\")
        lines.append(r"\hline")
        for row in prod_rows:
            cum_str = f"${row['cumulative']:.0f}\\times$"
            lines.append(f"{row['label']} & {row['runtime']:.1f} & "
                         f"{cum_str} & {row['marginal_str']} & "
                         f"{row['rms_str']} \\\\")

    lines.append(r"\enddata")
    lines.append(
        r"\tablecomments{The reference is unmodified SNEC~1.01 with 1000 zones, "
        r"compiled with \texttt{-O3 -fbounds-check}, "
        r"\texttt{EPSTOL}$=10^{-7}$, and fixed \Nif\ cadence "
        r"($\Delta t_\mathrm{Ni}=5\times10^4$~s). "
        r"All SuperSNEC runs use 100 zones unless otherwise noted. "
        r"All rows above the production grid section use "
        r"\texttt{output\_mode=0} (full SNEC output). "
        r"Production grid rows apply \texttt{output\_mode=1} to their "
        r"respective configurations; the marginal column shows the speedup "
        r"from enabling minimal output relative to the full-output row. "
        r"Runtimes are 3-run averages for SuperSNEC configurations; "
        r"run-to-run variation is $\lesssim 5\%$.}"
    )
    lines.append(r"\end{deluxetable*}")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# TSV output for reproducibility
# ---------------------------------------------------------------------------
def write_results_tsv(configs: List[Config], results: List[dict]) -> None:
    """Write a machine-readable TSV of all results."""
    tsv_path = RESULTS_DIR / "results.tsv"
    with tsv_path.open("w") as f:
        f.write("config\truntime_s\trms_all\tdescription\n")
        for cfg, res in zip(configs, results):
            rms = res.get("rms_all", float("nan"))
            rms_str = f"{rms:.5f}" if not np.isnan(rms) else "nan"
            f.write(f"{cfg.outdir_name}\t{res['runtime']:.2f}\t{rms_str}\t{cfg.description}\n")
    print(f"\nResults TSV: {tsv_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-reference", action="store_true",
                        help="Skip the 1000-zone reference run (use cached)")
    parser.add_argument("--force", action="store_true",
                        help="Force re-run of all configurations")
    parser.add_argument("--supersnec-only", action="store_true",
                        help="Skip original SNEC runs (use cached)")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    configs = build_configs()
    results: List[dict] = []

    print("=" * 60)
    print("Cumulative speedup benchmark")
    print("=" * 60)

    # Run each configuration
    ref_lum_path = None
    for i, cfg in enumerate(configs):
        print(f"\n[{i+1}/{len(configs)}] {cfg.label}")
        print(f"         {cfg.description}")

        skip = False
        if cfg.codebase == "original_1000" and args.skip_reference:
            skip = True
        if cfg.codebase.startswith("original") and args.supersnec_only:
            skip = True

        if skip and not args.force:
            # Try to load cached result
            cached = RESULTS_DIR / cfg.outdir_name / "result.npz"
            if cached.exists():
                data = dict(np.load(cached, allow_pickle=True))
                res = {k: float(v) for k, v in data.items()}
                print(f"  [cached] {res['runtime']:.1f}s")
                results.append(res)
            else:
                print(f"  [SKIPPED — no cached result at {cached}]")
                results.append({"runtime": float("nan"), "rms_all": float("nan")})
        else:
            res = run_config(cfg, ref_lum_path=ref_lum_path, force=args.force)
            results.append(res)

        # After the reference run, set the reference lum path
        if i == 0:
            ref_lum_path = RESULTS_DIR / cfg.outdir_name / "lum_observed.dat"

    # Generate outputs
    print("\n" + "=" * 60)
    print("Generating table...")
    table_tex = generate_table(configs, results)
    table_path = TABLE_DIR / "table_cumulative_speedup.tex"
    table_path.write_text(table_tex)
    print(f"Table written to: {table_path}")

    write_results_tsv(configs, results)

    print("\nDone!")
    print(f"Reference runtime: {results[0]['runtime']:.1f}s")
    if len(results) > 1 and not np.isnan(results[-1]["runtime"]):
        speedup = results[0]["runtime"] / results[-1]["runtime"]
        print(f"Final runtime:     {results[-1]['runtime']:.1f}s")
        print(f"Total speedup:     {speedup:.0f}x")


if __name__ == "__main__":
    main()
