#!/usr/bin/env python3
"""Generate adaptive-runtime validation LC figure (Figure 1).

Four-panel figure (2x2): top row = bolometric LC, bottom row = magnitude
residuals vs 1000-zone SNEC reference.  Left column = full LC (0-58 d),
right column = early zoom (0-3 d).

Models compared (all residuals measured against model 1):

  1. SNEC reference   -- original unmodified SNEC, N=1000, GridPattern.dat
  2. SNEC N=100       -- original unmodified SNEC, N=100,  GridPattern.dat
  3. SuperSNEC legacy -- SuperSNEC code, N=100, legacy_pattern grid
  4. SuperSNEC base   -- SuperSNEC code, N=100, adaptive_runtime, n_quad=70
  5. SuperSNEC fast   -- SuperSNEC code, N=60,  adaptive_runtime, n_quad=50

The script runs SuperSNEC to produce models 3-5 (setting Ni_quad_npoints
and smooth_ni_luminosity as runtime parameters).  Results are cached in
tmp/ and reused on subsequent calls.  Models 1-2 are pre-computed and
read from the original SNEC installation.

Usage:
  python3 Analysis/scripts/plot_adaptive_runtime_validation.py
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Tuple


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import (
    LogLocator, FuncFormatter, AutoMinorLocator, NullFormatter,
)
import numpy as np

# Publication-quality font settings
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

ROOT = Path(__file__).resolve().parents[2]
TMP = ROOT / "tmp"
DAY = 86400.0

PARAMS = ROOT / "parameters"

sys.path.insert(0, str(ROOT / "Analysis" / "scripts"))
from snec_metrics import write_parameters_file  # noqa: E402

# ---------------------------------------------------------------------------
# Model definitions and parameter documentation
# ---------------------------------------------------------------------------
# All models share:  profile = stripped_star.short, E_exp = 1e51 erg,
#                    M_Ni = 0.03 Msun, M_excised = 1.4 Msun
#
# Model 1: SNEC reference (1000 zones)
#   Code:       original unmodified SNEC-1.01
#   Grid:       from_file_by_mass (GridPattern.dat), N=1000
#   Ni mixing:  Ni_boundary_mass = 2.5 Msun (step function)
#   Ni heating: Ni_period = 5e4 s, n_quad = 150 (compile-time default)
#   Solver:     dtmax = 1e4 s
#   Runtime:    ~663 s
SNEC_REF_1000 = ROOT / "tmp" / "snec_original_1000z" / "lum_observed.dat"

# Model 2: SNEC N=100 (original code, 100 zones)
#   Code:       original unmodified SNEC-1.01
#   Grid:       from_file_by_mass (GridPattern.dat), N=100
#   Ni mixing:  Ni_boundary_mass = 2.5 Msun (step function)
#   Ni heating: Ni_period = 5e4 s, n_quad = 150 (compile-time default)
#   Solver:     dtmax = 1e4 s
SNEC_REF_100 = ROOT / "tmp" / "snec_original_100z" / "lum_observed.dat"

# Model 3: SuperSNEC legacy grid (100 zones)
#   Code:       SuperSNEC (this repo)
#   Grid:       legacy_pattern (GridPattern.dat), N=100
#   Ni mixing:  Ni_mix_fraction = 0.31, Ni_mix_kernel = 1 (box/step)
#   Ni heating: Ni_period = 5e4 s, Ni_period_max = 5.456e5 s,
#               Ni_fractional_change = 0.10, Ni_quad_npoints = 150
#   Ni smoothing: smooth_ni_luminosity = 1 (on)
#   Solver:     dtmax = 1e5 s
#   Grid params: (no adaptive grid — fixed legacy pattern)
SSNEC_LEGACY = dict(
    imax=100, nquad=150, label="legacy",
    extra_params={
        "grid_mode": '"legacy_pattern"',
        "Ni_fractional_change": "0.10d0",
        "smooth_ni_luminosity": "0",
    },
)

# Model 4: SuperSNEC baseline (100 zones, adaptive grid)
#   Code:       SuperSNEC (this repo)
#   Grid:       adaptive_runtime, N=100
#               grid_surface_alpha = 7.0, grid_relax_days = 5.0 (eff=5.0),
#               grid_update_interval_days = 1.0, grid_min_cell_frac = 1e-4
#   Ni mixing:  Ni_mix_fraction = 0.31, Ni_mix_kernel = 1 (box/step)
#   Ni heating: Ni_period = 5e4 s, Ni_period_max = 5.456e5 s,
#               Ni_fractional_change = 0.20, Ni_quad_npoints = 70
#   Ni smoothing: smooth_ni_luminosity = 1 (on)
#   Solver:     dtmax = 1e5 s
SSNEC_BASELINE = dict(imax=100, nquad=70, label="baseline", extra_params={
    "Ni_fractional_change": "0.20d0",
})

# Model 5: SuperSNEC fast (60 zones, adaptive grid)
#   Code:       SuperSNEC (this repo)
#   Grid:       adaptive_runtime, N=60
#               grid_surface_alpha = 7.0, grid_relax_days = 5.0 (eff=2.0),
#               grid_update_interval_days = 1.0, grid_min_cell_frac = 1e-4
#   Ni mixing:  Ni_mix_fraction = 0.31, Ni_mix_kernel = 1 (box/step)
#   Ni heating: Ni_period = 5e4 s, Ni_period_max = 5.456e5 s,
#               Ni_fractional_change = 0.20, Ni_quad_npoints = 50
#   Ni smoothing: smooth_ni_luminosity = 1 (on)
#   Solver:     dtmax = 1e5 s
SSNEC_FAST = dict(imax=60, nquad=50, label="fastest", extra_params={
    "Ni_fractional_change": "0.20d0",
})


# ---- Helpers ----

def read_lum(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Read lum_observed.dat, return (time_days, luminosity) with L>0."""
    t, l = [], []
    with path.open() as f:
        for line in f:
            s = line.split()
            if len(s) < 2:
                continue
            ti, li = float(s[0]), float(s[1])
            if li > 0.0:
                t.append(ti / DAY)
                l.append(li)
    return np.array(t), np.array(l)


def compute_residuals(
    t_ref: np.ndarray, l_ref: np.ndarray,
    t_run: np.ndarray, l_run: np.ndarray,
    ngrid: int = 2000, t_min: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray]:
    """Magnitude residual: dm = -2.5 * log10(L_run / L_ref)."""
    t_start = max(t_ref[0], t_run[0], t_min)
    t_end = min(t_ref[-1], t_run[-1])
    tg = np.linspace(t_start, t_end, ngrid)
    lr = np.interp(tg, t_ref, l_ref)
    lm = np.interp(tg, t_run, l_run)
    good = (lr > 0) & (lm > 0)
    dm = np.full_like(tg, np.nan)
    dm[good] = -2.5 * np.log10(lm[good] / lr[good])
    return tg, dm


def compute_rms(t: np.ndarray, dm: np.ndarray, t_lo: float, t_hi: float) -> float:
    """RMS of residuals in [t_lo, t_hi] days."""
    mask = (t >= t_lo) & (t <= t_hi) & ~np.isnan(dm)
    if not mask.any():
        return 0.0
    return float(np.sqrt(np.mean(dm[mask] ** 2)))


def run_model(
    imax: int, nquad: int, label: str,
    extra_params: dict[str, str] | None = None,
) -> Path:
    """Run SNEC with given imax and quadrature, return path to cached LC.

    extra_params: additional key=value overrides for the parameters file
                  (e.g. grid_mode, Ni_fractional_change).
    """
    cache = TMP / f"fig1_{label}.dat"
    if cache.exists():
        cache.unlink()
        print(f"  {label}: cleared stale cache, will re-run")

    TMP.mkdir(parents=True, exist_ok=True)
    outdir = TMP / f"fig1_{label}_run"
    outdir.mkdir(parents=True, exist_ok=True)

    # Build self-contained parameters file from baseline template.
    overrides = {
        "imax": str(imax),
        "Ni_quad_npoints": str(nquad),
        "smooth_ni_luminosity": "1",
        "outdir": f'"{outdir}/"',
        "grid_debug": "0",
    }
    if extra_params:
        overrides.update(extra_params)
    write_parameters_file(PARAMS, overrides)

    # Clean stale files
    for f in ROOT.glob("grid_diag_*.dat"):
        f.unlink()
    for f in ROOT.glob("fort.*"):
        f.unlink()

    # Run
    print(f"  {label}: running SNEC (N={imax}, n_quad={nquad})...")
    t0 = time.time()
    subprocess.run(["./snec"], capture_output=True, timeout=1200, cwd=ROOT)
    elapsed = time.time() - t0
    print(f"  {label}: {elapsed:.1f}s")

    # Copy LC
    lc_src = outdir / "lum_observed.dat"
    if lc_src.exists():
        shutil.copy2(lc_src, cache)
    else:
        raise FileNotFoundError(f"SNEC did not produce {lc_src}")

    return cache


# ---- Figure ----

def generate_figure(
    ref1000_path: Path,
    snec100_path: Path,
    legacy_path: Path,
    baseline_path: Path,
    fastest_path: Path,
    out_png: Path,
    tmax_early: float = 3.0,
    tmax_full: float = 58.0,
) -> None:
    t_ref, l_ref = read_lum(ref1000_path)
    t_s100, l_s100 = read_lum(snec100_path)
    t_leg, l_leg = read_lum(legacy_path)
    t_base, l_base = read_lum(baseline_path)
    t_fast, l_fast = read_lum(fastest_path)

    # Residuals (all vs 1000-zone reference)
    tg_s100, dm_s100 = compute_residuals(t_ref, l_ref, t_s100, l_s100)
    tg_leg, dm_leg = compute_residuals(t_ref, l_ref, t_leg, l_leg)
    tg_base, dm_base = compute_residuals(t_ref, l_ref, t_base, l_base)
    tg_fast, dm_fast = compute_residuals(t_ref, l_ref, t_fast, l_fast)

    # RMS for labels (t >= 0.1 d, matching snec_metrics rms_all cutoff)
    rms_s100 = compute_rms(tg_s100, dm_s100, 0.1, 60)
    rms_leg = compute_rms(tg_leg, dm_leg, 0.1, 60)
    rms_base = compute_rms(tg_base, dm_base, 0.1, 60)
    rms_fast = compute_rms(tg_fast, dm_fast, 0.1, 60)

    print(f"  RMS (all): SNEC-100={rms_s100:.3f}  legacy={rms_leg:.3f}"
          f"  baseline={rms_base:.3f}  fast={rms_fast:.3f}")

    # Styles
    ref_style = dict(color="#555555", lw=2.5, ls=(0, (5, 3)), zorder=1)
    s100_style = dict(color="#9467bd", lw=1.8, ls=(0, (1, 1)), zorder=2)
    leg_style = dict(color="#ff7f0e", lw=2.0, ls=(0, (3, 1.5, 1, 1.5)), zorder=3)
    base_style = dict(color="#1f77b4", lw=2.0, ls="-", zorder=5)
    fast_style = dict(color="#d62728", lw=2.0, ls="-", zorder=4)

    lbl_ref = "SNEC reference ($N$=1000)"
    lbl_s100 = f"SNEC ($N$=100, RMS={rms_s100:.3f})"
    lbl_leg = f"SuperSNEC fixed grid ($N$=100, RMS={rms_leg:.3f})"
    lbl_base = f"SuperSNEC baseline ($N$=100, $n_\\mathrm{{quad}}$=70, RMS={rms_base:.3f})"
    lbl_fast = f"SuperSNEC fast ($N$=60, $n_\\mathrm{{quad}}$=50, RMS={rms_fast:.3f})"

    FS_LABEL = 16
    FS_TICK = 14
    FS_LEGEND = 10.5

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[2, 1],
        hspace=0.0, wspace=0.15,
        left=0.09, right=0.98, top=0.97, bottom=0.09,
    )
    ax_lc_full = fig.add_subplot(gs[0, 0])
    ax_res_full = fig.add_subplot(gs[1, 0], sharex=ax_lc_full)
    ax_lc_early = fig.add_subplot(gs[0, 1])
    ax_res_early = fig.add_subplot(gs[1, 1], sharex=ax_lc_early)

    ylim_shared = (5e40, 5e42)

    def _scalar_log_fmt(x, pos):
        if x <= 0:
            return ""
        exp = int(np.floor(np.log10(x)))
        coeff = round(x / 10**exp, 1)
        if abs(coeff - 1.0) < 0.05:
            return r"$10^{%d}$" % exp
        if coeff == int(coeff):
            return r"$%d{\times}10^{%d}$" % (int(coeff), exp)
        return r"$%.1f{\times}10^{%d}$" % (coeff, exp)

    scalar_log_formatter = FuncFormatter(_scalar_log_fmt)

    # Plot LCs
    for ax, xlim in [(ax_lc_full, (0, tmax_full)), (ax_lc_early, (0, tmax_early))]:
        ax.plot(t_ref, l_ref, **ref_style)
        ax.plot(t_s100, l_s100, **s100_style)
        ax.plot(t_leg, l_leg, **leg_style)
        ax.plot(t_fast, l_fast, **fast_style)
        ax.plot(t_base, l_base, **base_style)
        ax.set_yscale("log")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim_shared)
        ax.set_ylabel(r"$L_{\mathrm{bol}}$ (erg s$^{-1}$)", fontsize=FS_LABEL)
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0), numticks=20))
        ax.yaxis.set_major_formatter(scalar_log_formatter)
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(1, 10), numticks=50))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="both", which="major", labelsize=FS_TICK, length=6)
        ax.tick_params(axis="both", which="minor", length=3)
        ax.tick_params(axis="x", labelbottom=False, which="both")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, alpha=0.3, which="major")

    # Plot residuals
    for ax, xlim in [(ax_res_full, (0, tmax_full)), (ax_res_early, (0, tmax_early))]:
        ax.axhspan(-0.1, 0.1, color="#aaddaa", alpha=0.35, zorder=0)
        ax.axhline(0, color="black", ls="--", lw=0.8, zorder=1)
        ax.plot(tg_s100, dm_s100, **s100_style)
        ax.plot(tg_leg, dm_leg, **leg_style)
        ax.plot(tg_fast, dm_fast, **fast_style)
        ax.plot(tg_base, dm_base, **base_style)
        ax.set_xlim(xlim)
        ax.set_ylim(-0.35, 0.35)
        ax.set_ylabel(r"$\Delta m$ (mag)", fontsize=FS_LABEL)
        ax.set_xlabel("Time since explosion (days)", fontsize=FS_LABEL)
        ax.tick_params(axis="both", which="major", labelsize=FS_TICK, length=6)
        ax.tick_params(axis="both", which="minor", length=3)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, alpha=0.3, which="major")

    ax_lc_early.set_ylabel("")
    ax_res_early.set_ylabel("")

    # Legend
    legend_handles = [
        Line2D([0], [0], color=ref_style["color"], lw=2.5,
               ls=(0, (4, 2)), label=lbl_ref),
        Line2D([0], [0], color=s100_style["color"], lw=2.0,
               ls=(0, (1, 1)), label=lbl_s100),
        Line2D([0], [0], color=leg_style["color"], lw=2.5,
               ls=(0, (3, 1.5, 1, 1.5)), label=lbl_leg),
        Line2D([0], [0], color=base_style["color"], lw=2.5,
               ls="-", label=lbl_base),
        Line2D([0], [0], color=fast_style["color"], lw=2.5,
               ls="-", label=lbl_fast),
    ]
    ax_lc_full.legend(
        handles=legend_handles, loc="upper right", fontsize=FS_LEGEND,
        frameon=False, handlelength=3.5,
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}")


def main() -> None:
    for p in [SNEC_REF_1000, SNEC_REF_100]:
        if not p.exists():
            print(f"ERROR: Reference LC not found: {p}", file=sys.stderr)
            sys.exit(1)

    # Run SuperSNEC models (cached in tmp/)
    legacy_lc = run_model(**SSNEC_LEGACY)
    baseline_lc = run_model(**SSNEC_BASELINE)
    fastest_lc = run_model(**SSNEC_FAST)

    # Generate figure
    out_png = ROOT / "Analysis" / "figures" / "adaptive_runtime_validation.png"
    generate_figure(
        SNEC_REF_1000, SNEC_REF_100,
        legacy_lc, baseline_lc, fastest_lc, out_png,
    )


if __name__ == "__main__":
    main()
