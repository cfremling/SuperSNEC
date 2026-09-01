#!/usr/bin/env python3
"""Analyse and plot Sedov convergence tests.

Reads SNEC .xg output from tmp/sedov_convergence/*, compares to the
exact Sedov solution, and produces:

  1. Error-vs-radius profiles for each variable (rho, v, p, eps)
     at the final time for every spatial-convergence run.

  2. L1-norm convergence plot: error vs number of zones N.

  3. Time-step convergence: error vs dtmax.

All figures are saved to Analysis/figures/.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Local module
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from sedov_exact import SedovSolution

ROOT = Path(__file__).resolve().parents[2]
TMPBASE = ROOT / "tmp" / "sedov_convergence"
FIGDIR = ROOT / "Analysis" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

# Sedov parameters (must match run_sedov_convergence.py)
GAMMA = 1.4
E0 = 1.0
RHO0 = 1.0
T_EVAL = 4.0    # final time of the Sedov runs
R_MIN = 1.0     # start of plot/error range (cf. Tasker et al. 2008 Fig. 6)
R_MAX = 2.0

# Variables to compare  (xg_file, label, log_scale)
VARIABLES = [
    ("rho",   r"$\rho$",  "density"),
    ("vel",   r"$v$",     "velocity"),
    ("press", r"$p$",     "pressure"),
    ("eps",   r"$\epsilon$", "spec. int. energy"),
]


# ------------------------------------------------------------------
# I/O helpers
# ------------------------------------------------------------------

def parse_xg(path: Path):
    """Parse a SNEC .xg file.  Return list of (time, col1[], col2[])."""
    snapshots = []
    current_time = None
    c1, c2 = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('"Time'):
                if current_time is not None:
                    snapshots.append((current_time, np.array(c1), np.array(c2)))
                current_time = float(line.split("=")[1])
                c1, c2 = [], []
            else:
                parts = line.split()
                if len(parts) >= 2:
                    c1.append(float(parts[0]))
                    c2.append(float(parts[1]))
    if current_time is not None:
        snapshots.append((current_time, np.array(c1), np.array(c2)))
    return snapshots


def last_snapshot(rundir: Path, var: str):
    """Return (mass[], value[]) for the last time-snapshot of *var*.xg."""
    snaps = parse_xg(rundir / "Data" / f"{var}.xg")
    if not snaps:
        return None, None
    _, mass, val = snaps[-1]
    return mass, val


def get_radii(rundir: Path):
    """Return (mass[], radius[]) from the last snapshot of radius.xg."""
    return last_snapshot(rundir, "radius")


# ------------------------------------------------------------------
# Exact solution
# ------------------------------------------------------------------

def exact_at_radii(r_arr):
    """Return exact (rho, vel, press, eps) at given radii, t=T_EVAL."""
    sol = SedovSolution(gamma=GAMMA, E0=E0, rho0=RHO0)
    return sol.evaluate(r_arr, T_EVAL)


# ------------------------------------------------------------------
# Error computation
# ------------------------------------------------------------------

def relative_error(exact, numerical):
    """| (exact - numerical) / exact |.

    Where exact = 0, falls back to |numerical| / max(|exact|) so the
    error is still plottable (e.g. numerical velocity beyond the shock).
    """
    peak = np.nanmax(np.abs(exact))
    if peak == 0:
        return np.full_like(exact, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.abs((exact - numerical) / exact)
    # Where exact is zero, use absolute error / peak
    zero_mask = np.abs(exact) < 1e-30
    err[zero_mask] = np.abs(numerical[zero_mask]) / peak
    err[~np.isfinite(err)] = np.nan
    return err


def interpolate_to_common_grid(r_num, val_num, n_common=3000):
    """Interpolate a numerical profile onto a fine uniform grid
    covering the valid analytical range [R_MIN, R_MAX].
    """
    r_common = np.linspace(R_MIN, R_MAX, n_common)
    val_interp = np.interp(r_common, r_num, val_num)
    return r_common, val_interp


def interpolate_run_to_grid(data, r_common):
    """Interpolate all variables from a run onto r_common.

    Computes eps = p/((gamma-1)*rho) from the interpolated p and rho
    instead of interpolating eps directly, to avoid artifacts where
    grids have large gaps (eps = p/rho amplifies interpolation errors).
    """
    result = {}
    for var in ("rho", "vel", "press"):
        result[var] = np.interp(r_common, data["r"], data[var])
    # Derive eps from interpolated p and rho
    eps = np.zeros_like(r_common)
    mask = result["rho"] > 1e-30
    eps[mask] = result["press"][mask] / ((GAMMA - 1.0) * result["rho"][mask])
    result["eps"] = eps
    return result


def l1_norm(exact, numerical, r_arr):
    """Mean relative error at cell positions (unweighted).

    Same error definition as relative_error() used in the 4-panel
    sub-panels: |(exact - numerical) / exact|, averaged over cells
    in [R_MIN, 0.95*r_shock].
    """
    sol = SedovSolution(gamma=GAMMA, E0=E0, rho0=RHO0)
    r_shock = sol.shock_radius(T_EVAL)
    err = relative_error(exact, numerical)
    mask = (np.isfinite(err)
            & (r_arr >= R_MIN)
            & (r_arr < 0.95 * r_shock))
    if mask.sum() == 0:
        return np.nan
    return np.mean(err[mask])


# ------------------------------------------------------------------
# Spatial convergence runs
# ------------------------------------------------------------------

# Draw order matters: thick Orig. SNEC first, then thin dashed SuperSNEC on top.
# Identical curves show as dashes riding on a thick solid line.
#   (name, label, N, color, linestyle, linewidth, zorder)
ALL_SPATIAL_RUNS = [
    # (name, label, N, color, linestyle, linewidth, zorder)
    # SNEC and SuperSNEC on the same grid produce identical results;
    # we plot only SuperSNEC and note the equivalence in the text.
    ("sedov_static_100z",        "SNEC/SuperSNEC N=100 (legacy)",   100,  "#1f77b4", "-",   1.0, 1),
    ("sedov_static_1000z",       "SNEC/SuperSNEC N=1000 (legacy)",  1000, "#ff7f0e", "-",   1.0, 1),
    ("sedov_adaptive_100z",      "SuperSNEC N=100 (adaptive)",      100,  "#2ca02c", "--",  1.0, 2),
    ("sedov_adaptive_1000z",     "SuperSNEC N=1000 (adaptive)",     1000, "#d62728", "--",  1.0, 2),
]

TIMESTEP_GROUPS = {
    # label: (run_prefix, marker, color, ls, dt_mult)
    "N=1000 legacy (SNEC=SuperSNEC)":
        ("sedov_dt", "o", "C0", "-", 0.95),
    "N=1000 surf-conc (SuperSNEC)":
        ("sedov_dt_sc", "D", "C2", "-.", 1.05),
}

DT_VALUES = [0.1, 0.01, 0.001]


def load_run(name):
    """Load the final-snapshot profiles for a run.

    Returns dict with keys 'r', 'rho', 'vel', 'press', 'eps' (numpy arrays),
    or None if the run directory doesn't exist.
    """
    rundir = TMPBASE / name
    if not rundir.exists():
        print(f"  WARNING: {rundir} not found, skipping")
        return None

    mass_r, r_arr = get_radii(rundir)
    if r_arr is None:
        print(f"  WARNING: no radius.xg in {rundir}")
        return None

    data = {"r": r_arr}
    for var, _, _ in VARIABLES:
        mass_v, val = last_snapshot(rundir, var)
        if val is not None:
            data[var] = val
        else:
            data[var] = np.full_like(r_arr, np.nan)

    return data


# ------------------------------------------------------------------
# Plotting: error vs radius
# ------------------------------------------------------------------

def plot_error_profiles():
    """One figure per variable: |relative error| vs radius for all spatial runs.

    Numerical profiles are interpolated onto a fine common grid so that
    the error is visible everywhere, even for grids with very few
    interior cells.
    """
    sol = SedovSolution(gamma=GAMMA, E0=E0, rho0=RHO0)
    r_shock = sol.shock_radius(T_EVAL)

    for var, tex, desc in VARIABLES:
        fig, ax = plt.subplots(figsize=(7, 4.5))

        for name, label, N, color, ls, lw, zo in ALL_SPATIAL_RUNS:
            data = load_run(name)
            if data is None:
                continue

            r_num = data["r"]
            val_num = data[var]
            mask = (r_num >= R_MIN) & (r_num <= R_MAX)
            r_m = r_num[mask]
            v_m = val_num[mask]

            # Exact solution at cell positions
            rho_ex, vel_ex, press_ex, eps_ex = exact_at_radii(r_m)
            exact_map = {"rho": rho_ex, "vel": vel_ex,
                         "press": press_ex, "eps": eps_ex}

            err = relative_error(exact_map[var], v_m)
            ax.semilogy(r_m, err, color=color, ls=ls, lw=lw,
                        label=label, alpha=0.9, zorder=zo)

        ax.axvline(r_shock, color="gray", ls=":", lw=0.7, label="shock")
        ax.set_xlabel("Radius [cm]")
        ax.set_ylabel(f"|Error| / max({tex})")
        ax.set_title(f"Sedov test: {desc} error vs radius  (t = {T_EVAL} s)")
        ax.legend(fontsize=7, loc="upper left", frameon=False)
        ax.set_xlim(R_MIN, R_MAX)
        ax.set_ylim(1e-6, 1e1)
        fig.tight_layout()
        outpath = FIGDIR / f"sedov_error_{var}.pdf"
        fig.savefig(outpath, dpi=150)
        fig.savefig(outpath.with_suffix(".png"), dpi=150)
        plt.close(fig)
        print(f"  Saved {outpath.name}")


# ------------------------------------------------------------------
# Plotting: L1-norm convergence
# ------------------------------------------------------------------

def plot_l1_convergence():
    """L1 relative error vs N for each variable, for static and adaptive grids."""
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
    axes = axes.ravel()

    # Distinct markers, colors, and small N-offsets so identical points
    # don't hide each other
    group_style = {
        "Orig. SNEC (legacy)":       dict(color="C0", marker="o", ms=8, ls="-",
                                           lw=1.5, N_mult=0.90),
        "SuperSNEC (legacy)":         dict(color="C1", marker="s", ms=7, ls="--",
                                           lw=1.5, N_mult=0.97),
        "SuperSNEC (adaptive)":    dict(color="C3", marker="D", ms=6, ls="-.",
                                           lw=1.5, N_mult=1.03),
    }

    for idx, (var, tex, desc) in enumerate(VARIABLES):
        ax = axes[idx]

        groups = {
            "Orig. SNEC (legacy)":  [("sedov_orig_static_100z", 100),
                                      ("sedov_orig_static_1000z", 1000)],
            "SuperSNEC (legacy)":    [("sedov_static_100z", 100),
                                      ("sedov_static_1000z", 1000)],
            "SuperSNEC (adaptive)": [("sedov_adaptive_100z", 100),
                                        ("sedov_adaptive_1000z", 1000)],
        }

        last_Ns, last_errs = None, None
        for glabel, runs in groups.items():
            sty = group_style[glabel]
            Ns, errs = [], []
            for rname, N in runs:
                data = load_run(rname)
                if data is None:
                    continue
                r_num = data["r"]
                val_num = data[var]
                rho_ex, vel_ex, press_ex, eps_ex = exact_at_radii(r_num)
                exact_map = {"rho": rho_ex, "vel": vel_ex,
                             "press": press_ex, "eps": eps_ex}
                e = l1_norm(exact_map[var], val_num, r_num)
                if not np.isnan(e):
                    Ns.append(N)
                    errs.append(e)

            if Ns:
                # Slight horizontal offset so identical points are visible
                Ns_plot = [n * sty["N_mult"] for n in Ns]
                ax.loglog(Ns_plot, errs, marker=sty["marker"],
                          color=sty["color"], ls=sty["ls"], lw=sty["lw"],
                          markersize=sty["ms"], label=glabel, zorder=3)
                last_Ns, last_errs = Ns, errs

        # Reference slopes anchored at a fixed point
        Nref = np.array([80, 1200])
        for order, ls_ref in [(1, ":"), (2, "--")]:
            # Anchor at N=100, y=0.1 so lines are visible in the plot range
            ref = 0.1 * (100.0 / Nref)**order
            ax.loglog(Nref, ref, color="gray", ls=ls_ref, lw=0.6,
                      label=f"$N^{{-{order}}}$" if idx == 0 else None)

        ax.set_ylabel(f"L1 error ({tex})")
        ax.set_ylim(1e-3, 2e-1)
        ax.set_title(desc, fontsize=10)
        if idx >= 2:
            ax.set_xlabel("Number of zones N")
        ax.legend(fontsize=7, frameon=False)

    fig.suptitle(f"Sedov spatial convergence (t = {T_EVAL} s)", fontsize=12)
    fig.tight_layout()
    outpath = FIGDIR / "sedov_convergence_spatial.pdf"
    fig.savefig(outpath, dpi=150)
    fig.savefig(outpath.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {outpath.name}")


# ------------------------------------------------------------------
# Plotting: time-step convergence
# ------------------------------------------------------------------

def plot_timestep_convergence():
    """L1 error vs dtmax for N=1000, multiple grid types."""
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
    axes = axes.ravel()

    for idx, (var, tex, desc) in enumerate(VARIABLES):
        ax = axes[idx]

        for glabel, (prefix, marker, color, ls, dt_mult) in TIMESTEP_GROUPS.items():
            dts, errs = [], []
            for dtmax in DT_VALUES:
                rname = f"{prefix}_{dtmax}"
                data = load_run(rname)
                if data is None:
                    continue
                r_num = data["r"]
                val_num = data[var]
                rho_ex, vel_ex, press_ex, eps_ex = exact_at_radii(r_num)
                exact_map = {"rho": rho_ex, "vel": vel_ex,
                             "press": press_ex, "eps": eps_ex}
                e = l1_norm(exact_map[var], val_num, r_num)
                if not np.isnan(e):
                    dts.append(dtmax)
                    errs.append(e)

            if dts:
                dts_plot = [d * dt_mult for d in dts]
                ax.loglog(dts_plot, errs, marker=marker, color=color,
                          ls=ls, markersize=6, lw=1.5, label=glabel)

        # Reference slopes
        dtref = np.array([5e-4, 0.2])
        for order, ls_ref in [(1, ":"), (2, "--")]:
            ref = 4e-2 * (0.1 / dtref)**order
            ax.loglog(dtref, ref, color="gray", ls=ls_ref, lw=0.6,
                      label=f"$\\Delta t^{{{order}}}$" if idx == 0 else None)

        ax.set_ylabel(f"L1 error ({tex})")
        ax.set_title(desc, fontsize=10)
        if idx >= 2:
            ax.set_xlabel(r"$\Delta t_{\max}$ [s]")
        ax.legend(fontsize=7, frameon=False)

    fig.suptitle(f"Sedov temporal convergence (N=1000, t = {T_EVAL} s)", fontsize=12)
    fig.tight_layout()
    outpath = FIGDIR / "sedov_convergence_temporal.pdf"
    fig.savefig(outpath, dpi=150)
    fig.savefig(outpath.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {outpath.name}")


# ------------------------------------------------------------------
# Plotting: exact + numerical overlay
# ------------------------------------------------------------------

def _get_exact_and_interp(sol, r_fine):
    """Compute exact profiles on a fine grid (cached-style helper)."""
    rho_ex, vel_ex, press_ex, eps_ex = sol.evaluate(r_fine, T_EVAL)
    return {"rho": rho_ex, "vel": vel_ex, "press": press_ex, "eps": eps_ex}


def _plot_one_profile(ax, var, tex, desc, r_fine, exact_profiles, r_shock,
                      use_log=False):
    """Plot one variable: exact (under) + all numerical runs (interpolated)."""
    # Analytical first (under everything)
    exact_vals = exact_profiles[var].copy()
    if use_log:
        exact_vals[exact_vals <= 0] = np.nan
    ax.plot(r_fine, exact_vals, "k-", lw=1.0,
            label="Analytical", zorder=1)

    for name, label, N, color, ls, lw, zo in ALL_SPATIAL_RUNS:
        data = load_run(name)
        if data is None:
            continue
        r_num = data["r"]
        val = data[var].copy()
        mask = (r_num >= R_MIN) & (r_num <= R_MAX)
        if use_log:
            val[val <= 0] = np.nan
        ax.plot(r_num[mask], val[mask], color=color, ls=ls, lw=0.7,
                label=label, alpha=0.9, zorder=zo + 1)

    # Shock position is obvious from the profiles; no marker needed
    ax.set_xlim(R_MIN, R_MAX)


def plot_profiles():
    """Individual profile figures (one per variable)."""
    sol = SedovSolution(gamma=GAMMA, E0=E0, rho0=RHO0)
    r_shock = sol.shock_radius(T_EVAL)
    r_fine = np.linspace(R_MIN, R_MAX, 5000)
    exact_profiles = _get_exact_and_interp(sol, r_fine)

    for var, tex, desc in VARIABLES:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        _plot_one_profile(ax, var, tex, desc, r_fine, exact_profiles, r_shock)
        ax.set_xlabel("Radius [cm]")
        ax.set_ylabel(f"{tex}  [{desc}]")
        ax.set_title(f"Sedov test: {desc}  (t = {T_EVAL} s)")
        ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        outpath = FIGDIR / f"sedov_profile_{var}.pdf"
        fig.savefig(outpath, dpi=150)
        fig.savefig(outpath.with_suffix(".png"), dpi=150)
        plt.close(fig)
        print(f"  Saved {outpath.name}")


def plot_profiles_4panel():
    """Combined 4-panel figure with error sub-panels (Tasker et al. 2008 style).

    Each column has a profile panel (top) and a relative-error panel (bottom),
    arranged as a 4x2 grid of subplots with shared x-axes.
    """
    sol = SedovSolution(gamma=GAMMA, E0=E0, rho0=RHO0)
    r_shock = sol.shock_radius(T_EVAL)
    r_fine = np.linspace(R_MIN, R_MAX, 5000)
    exact_profiles = _get_exact_and_interp(sol, r_fine)

    panel_order = [
        ("rho",   r"$\rho$",  "density",     False),
        ("vel",   r"$v$",     "velocity",    False),
        ("press", r"$p$",     "pressure",    False),
        ("eps",   r"$U$",     "int. energy", True),
    ]

    # Match style of plot_adaptive_runtime_validation.py
    FS_LABEL = 16
    FS_TICK = 14
    FS_LEGEND = 10.5

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

    # 2x2 outer grid, each cell contains a profile+error pair with hspace=0
    # Matches layout style of plot_adaptive_runtime_validation.py
    fig = plt.figure(figsize=(14, 12))
    outer = fig.add_gridspec(
        2, 2,
        hspace=0.12, wspace=0.15,
        left=0.09, right=0.98, top=0.97, bottom=0.07,
    )

    for idx, (var, tex, desc, use_log) in enumerate(panel_order):
        row, col = divmod(idx, 2)
        inner = outer[row, col].subgridspec(
            2, 1, height_ratios=[3, 1], hspace=0.0)
        ax_prof = fig.add_subplot(inner[0])
        ax_err = fig.add_subplot(inner[1], sharex=ax_prof)
        r_base = row * 2  # for conditional x-label below

        # --- Profile panel ---
        _plot_one_profile(ax_prof, var, tex, desc, r_fine,
                          exact_profiles, r_shock, use_log=use_log)
        ax_prof.set_ylabel(tex, fontsize=FS_LABEL)
        if use_log:
            ax_prof.set_yscale("log")
            ax_prof.set_ylim(7e-4, 5e0)
        plt.setp(ax_prof.get_xticklabels(), visible=False)
        ax_prof.tick_params(axis="both", which="major",
                            labelsize=FS_TICK, length=6)
        ax_prof.tick_params(axis="both", which="minor", length=3)
        ax_prof.grid(True, alpha=0.3, which="major")

        # --- Error sub-panel ---
        for name, label, N, color, ls, lw, zo in ALL_SPATIAL_RUNS:
            data = load_run(name)
            if data is None:
                continue
            r_num = data["r"]
            val_num = data[var]
            mask = (r_num >= R_MIN) & (r_num <= R_MAX)
            r_m = r_num[mask]
            v_m = val_num[mask]

            rho_ex, vel_ex, press_ex, eps_ex = exact_at_radii(r_m)
            exact_map = {"rho": rho_ex, "vel": vel_ex,
                         "press": press_ex, "eps": eps_ex}
            err = relative_error(exact_map[var], v_m)
            # Smooth with a running median to show trends
            if len(err) > 15:
                from scipy.ndimage import median_filter
                err_smooth = median_filter(err, size=max(5, len(err)//30))
            else:
                err_smooth = err
            ax_err.semilogy(r_m, err_smooth, color=color, ls=ls, lw=1.0,
                            alpha=0.9, zorder=zo + 1)

        ax_err.axvline(r_shock, color="black", ls="--", lw=1.0)
        ax_err.set_ylim(1e-4, 2e0)
        ax_err.set_xlim(R_MIN, R_MAX)
        ax_err.set_xlabel(r"$r$ (cm)", fontsize=FS_LABEL)
        ax_err.set_ylabel(r"$\delta$", fontsize=FS_LABEL)
        ax_err.tick_params(axis="both", which="major",
                           labelsize=FS_TICK, length=6)
        ax_err.tick_params(axis="both", which="minor", length=3)
        ax_err.grid(True, alpha=0.3, which="major")

        ax_prof.minorticks_on()
        ax_err.minorticks_on()

        # Only show x-label on bottom error panels (rows 1 and 3)
        if r_base < 2:
            ax_err.tick_params(axis="x", labelbottom=False)

    # Legend in the first profile panel
    fig.axes[0].legend(fontsize=FS_LEGEND, frameon=False, loc="upper left")

    outpath = FIGDIR / "sedov_profiles_4panel.pdf"
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath.name}")


# ------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------

def print_summary_table():
    """Print L1 errors for all runs and variables."""
    sol = SedovSolution(gamma=GAMMA, E0=E0, rho0=RHO0)

    header = f"{'Run':<35s}"
    for _, tex, desc in VARIABLES:
        header += f"  {desc:>18s}"
    print("\n" + header)
    print("-" * len(header))

    all_runs = [(n, l, N) for n, l, N, *_ in ALL_SPATIAL_RUNS]
    for glabel, (prefix, _, _, _, _) in TIMESTEP_GROUPS.items():
        for dtmax in DT_VALUES:
            rname = f"{prefix}_{dtmax}"
            all_runs.append((rname, f"{glabel} dt={dtmax}", 1000))

    for rname, label, N in all_runs:
        data = load_run(rname)
        if data is None:
            continue
        row = f"{label:<35s}"
        r_num = data["r"]
        rho_ex, vel_ex, press_ex, eps_ex = exact_at_radii(r_num)
        exact_map = {"rho": rho_ex, "vel": vel_ex,
                     "press": press_ex, "eps": eps_ex}
        for var, _, _ in VARIABLES:
            e = l1_norm(exact_map[var], data[var], r_num)
            row += f"  {e:18.4e}"
        print(row)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    print("=" * 60)
    print("Sedov Convergence Analysis")
    print("=" * 60)

    print("\nPlotting 4-panel profile figure...")
    plot_profiles_4panel()

    print("\nPlotting individual profile overlays...")
    plot_profiles()

    print("\nPlotting error profiles...")
    plot_error_profiles()

    print("\nPlotting spatial convergence...")
    plot_l1_convergence()

    print("\nPlotting temporal convergence...")
    plot_timestep_convergence()

    print_summary_table()


if __name__ == "__main__":
    main()
