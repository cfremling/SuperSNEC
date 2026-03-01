#!/usr/bin/env python3
"""Generate SN light-curve and velocity comparison figure (Figure 2).

Four-panel figure (2×2): top row = SN 2011dh, bottom row = SN 1993J.
Left column = bolometric luminosity (0–100 d), right column = photospheric
velocity (0–30 d).

Data sources:
  SN 2011dh  — Ergon et al. (2014, 2015): bolometric LC and Fe II 5169 velocity
  SN 1993J   — Richmond et al. (1994): bolometric LC; UNLV spectra: Fe II 5169

SuperSNEC models from Analysis/lc_comparisons/ (pre-computed GUI runs).

Usage:
  python3 Analysis/scripts/plot_lc_comparison_sne.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import (
    LogLocator, FuncFormatter, AutoMinorLocator, NullFormatter,
)
import numpy as np

# ---------------------------------------------------------------------------
# Publication-quality font settings (matching Figure 1)
# ---------------------------------------------------------------------------
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
LC_DIR = ROOT / "Analysis" / "lc_comparisons"
DAY = 86400.0

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------
# SN 2011dh
SN2011DH_DIR = LC_DIR / "sn2011dh"
SN2011DH_MODEL = (
    SN2011DH_DIR
    / "m16_m14_p160_profile50_R_100_Mex_1.98_Ek_1.00_Ni_0.074_20260224_165305"
    / "gui_run_2145" / "Data"
)
SN2011DH_OBS_LBOL = SN2011DH_DIR / "obs_lbol.dat"
SN2011DH_OBS_VEL = SN2011DH_DIR / "obs_vel.dat"

# SN 1993J
SN1993J_DIR = LC_DIR / "sn1993j"
SN1993J_MODEL = (
    SN1993J_DIR
    / "m15_m13_p240_profile52_R_400_Mex_1.63_Ek_1.35_Ni_0.100_20260223_170326"
    / "gui_run_2002" / "Data"
)
SN1993J_OBS_LBOL = (
    SN1993J_DIR / "SN1993J_bolometric_lc_uvoir_bb_EBV0.080_RV3.1.csv"
)
SN1993J_OBS_VEL = (
    SN1993J_DIR / "SN1993J_feii5169_velocity_from_unlv_dif_quality.csv"
)


# ---------------------------------------------------------------------------
# Data readers
# ---------------------------------------------------------------------------
def read_snec_lum(path: Path):
    """Read lum_observed.dat → (time_days, L_erg_s) with L > 0."""
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


def read_snec_vel(path: Path):
    """Read vel_photo.dat → (time_days, v_km_s) with v > 0."""
    t, v = [], []
    with path.open() as f:
        for line in f:
            s = line.split()
            if len(s) < 2:
                continue
            ti, vi = float(s[0]), float(s[1])
            if vi > 0.0:
                t.append(ti / DAY)
                v.append(vi / 1e5)  # cm/s → km/s
    return np.array(t), np.array(v)


def read_2011dh_obs_lbol(path: Path):
    """Read obs_lbol.dat → (t_days, L_erg_s, sigma_L_erg_s)."""
    t, l, sig = [], [], []
    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            s = line.split()
            if len(s) < 3:
                continue
            td = float(s[0])
            logL = float(s[1])
            sig_logL = float(s[2])
            L = 10**logL
            # Propagate log error: sigma_L = L * ln(10) * sigma_logL
            sL = L * np.log(10) * sig_logL
            t.append(td)
            l.append(L)
            sig.append(sL)
    return np.array(t), np.array(l), np.array(sig)


def read_2011dh_obs_vel(path: Path):
    """Read obs_vel.dat → (t_days, v_km_s, sigma_km_s)."""
    t, v, sig = [], [], []
    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            s = line.split()
            if len(s) < 3:
                continue
            t.append(float(s[0]))
            v.append(float(s[1]))
            sig.append(float(s[2]))
    return np.array(t), np.array(v), np.array(sig)


def read_1993j_obs_lbol(path: Path):
    """Read SN1993J CSV → (t_days, L_erg_s, sigma_L_erg_s)."""
    t, l, sig = [], [], []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            td = float(row["t_days"])
            logL = float(row["log10_L_uvoir_bb_erg_s"])
            sig_logL = float(row["log10_L_err_stat_dex"])
            if sig_logL > 0.15:
                continue
            L = 10**logL
            sL = L * np.log(10) * sig_logL
            t.append(td)
            l.append(L)
            sig.append(sL)
    return np.array(t), np.array(l), np.array(sig)


def read_1993j_obs_vel(path: Path):
    """Read SN1993J velocity CSV → (t_days, v_km_s,v_err_km_s)."""
    t, v, vs = [], [], []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            t.append(float(row["days_since_2449076.7"]))
            v.append(float(row["v_feii_km_s"]))
            vs.append(float(row["v_err_km_s"]))
    return np.array(t), np.array(v), np.array(vs)


# ---------------------------------------------------------------------------
# Figure generation
# ---------------------------------------------------------------------------
def generate_figure(out_path: Path) -> None:
    # Read all data
    t_mod_11dh, l_mod_11dh = read_snec_lum(SN2011DH_MODEL / "lum_observed.dat")
    t_vmod_11dh, v_mod_11dh = read_snec_vel(SN2011DH_MODEL / "vel_photo.dat")
    t_obs_11dh, l_obs_11dh, sl_obs_11dh = read_2011dh_obs_lbol(SN2011DH_OBS_LBOL)
    t_vobs_11dh, v_obs_11dh, sv_obs_11dh = read_2011dh_obs_vel(SN2011DH_OBS_VEL)

    t_mod_93j, l_mod_93j = read_snec_lum(SN1993J_MODEL / "lum_observed.dat")
    t_vmod_93j, v_mod_93j = read_snec_vel(SN1993J_MODEL / "vel_photo.dat")
    t_obs_93j, l_obs_93j, sl_obs_93j = read_1993j_obs_lbol(SN1993J_OBS_LBOL)
    t_vobs_93j, v_obs_93j, sv_obs_93j = read_1993j_obs_vel(SN1993J_OBS_VEL)

    # Explosion-time offsets: shift model curves relative to obs reference epoch
    # SN 2011dh: explosion 1 day later  → model shifts +1 d
    # SN 1993J:  explosion 4 days earlier → model shifts −4 d
    t_mod_11dh += 0.6
    t_vmod_11dh += 0.6
    t_mod_93j -= 4.0
    t_vmod_93j -= 4.0

    # ---- Styles ----
    model_style = dict(color="#1f77b4", lw=2.0, ls="-", zorder=3)
    obs_lc_style = dict(
        fmt="o", color="#333333", ms=3.5, mew=0.0, ecolor="#888888",
        elinewidth=0.8, capsize=0, zorder=5,
    )
    obs_vel_style = dict(
        fmt="s", color="#333333", ms=4.0, mew=0.0, ecolor="#888888",
        elinewidth=0.8, capsize=1.5, zorder=5,
    )

    FS_LABEL = 16
    FS_TICK = 14
    FS_LEGEND = 11
    FS_ANNOT = 15

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(
        2, 2,
        hspace=0.25, wspace=0.30,
        left=0.09, right=0.98, top=0.97, bottom=0.08,
    )
    ax_lc_11dh = fig.add_subplot(gs[0, 0])
    ax_vel_11dh = fig.add_subplot(gs[0, 1])
    ax_lc_93j = fig.add_subplot(gs[1, 0])
    ax_vel_93j = fig.add_subplot(gs[1, 1])

    # ---- Log-scale luminosity formatter (matching Figure 1) ----
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

    # ---- Plot LC panels ----
    for ax, t_mod, l_mod, t_obs, l_obs, sl_obs, sn_name in [
        (ax_lc_11dh, t_mod_11dh, l_mod_11dh,
         t_obs_11dh, l_obs_11dh, sl_obs_11dh, "SN 2011dh"),
        (ax_lc_93j, t_mod_93j, l_mod_93j,
         t_obs_93j, l_obs_93j, sl_obs_93j, "SN 1993J"),
    ]:
        ax.plot(t_mod, l_mod, **model_style)
        ax.errorbar(t_obs, l_obs, yerr=sl_obs, **obs_lc_style)
        ax.set_yscale("log")
        ax.set_xlim(0, 100)
        ax.set_ylabel(r"$L_{\mathrm{bol}}$ (erg s$^{-1}$)", fontsize=FS_LABEL)
        ax.set_xlabel("Time since explosion (days)", fontsize=FS_LABEL)
        ax.yaxis.set_major_locator(
            LogLocator(base=10, subs=(1.0, 2.0, 5.0), numticks=20)
        )
        ax.yaxis.set_major_formatter(scalar_log_formatter)
        ax.yaxis.set_minor_locator(
            LogLocator(base=10, subs=np.arange(1, 10), numticks=50)
        )
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis="both", which="major", labelsize=FS_TICK, length=6)
        ax.tick_params(axis="both", which="minor", length=3)
        ax.grid(True, alpha=0.3, which="major")
        ax.text(
            0.95, 0.95, sn_name, transform=ax.transAxes,
            fontsize=FS_ANNOT, fontweight="bold",
            ha="right", va="top",
        )

    # Set per-SN y-limits
    ax_lc_11dh.set_ylim(1e41, 4e42)
    ax_lc_93j.set_ylim(1e41, 1e43)

    # ---- Plot velocity panels ----
    for ax, t_mod, v_mod, t_obs, v_obs, sv_obs, sn_name in [
        (ax_vel_11dh, t_vmod_11dh, v_mod_11dh,
         t_vobs_11dh, v_obs_11dh, sv_obs_11dh, "SN 2011dh"),
        (ax_vel_93j, t_vmod_93j, v_mod_93j,
         t_vobs_93j, v_obs_93j, sv_obs_93j, "SN 1993J"),
    ]:
        ax.plot(t_mod, v_mod / 1e3, **model_style)  # km/s → 10³ km/s
        ax.errorbar(
            t_obs, v_obs / 1e3, yerr=sv_obs / 1e3, **obs_vel_style
        )
        ax.set_xlim(0, 30)
        ax.set_ylabel(
            r"$v_{\mathrm{ph}}$ ($10^3$ km s$^{-1}$)", fontsize=FS_LABEL
        )
        ax.set_xlabel("Time since explosion (days)", fontsize=FS_LABEL)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis="both", which="major", labelsize=FS_TICK, length=6)
        ax.tick_params(axis="both", which="minor", length=3)
        ax.grid(True, alpha=0.3, which="major")
        ax.text(
            0.95, 0.95, sn_name, transform=ax.transAxes,
            fontsize=FS_ANNOT, fontweight="bold",
            ha="right", va="top",
        )

    ax_vel_11dh.set_ylim(0, 15)
    ax_vel_93j.set_ylim(0, 25)

    # ---- Legend (shared, placed in top-left LC panel) ----
    legend_handles = [
        Line2D([0], [0], color=model_style["color"], lw=2.5, ls="-",
               label=r"SuperSNEC ($N$=100)"),
        Line2D([0], [0], color="#333333", marker="o", ms=4.5, lw=0,
               label="Observations"),
    ]
    ax_lc_11dh.legend(
        handles=legend_handles, loc="lower left", fontsize=FS_LEGEND,
        frameon=False, handlelength=2.5,
    )

    # ---- Save ----
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main() -> None:
    # Verify data exists
    for p in [
        SN2011DH_MODEL / "lum_observed.dat",
        SN2011DH_MODEL / "vel_photo.dat",
        SN2011DH_OBS_LBOL,
        SN2011DH_OBS_VEL,
        SN1993J_MODEL / "lum_observed.dat",
        SN1993J_MODEL / "vel_photo.dat",
        SN1993J_OBS_LBOL,
        SN1993J_OBS_VEL,
    ]:
        if not p.exists():
            raise FileNotFoundError(f"Missing data file: {p}")

    out_pdf = ROOT / "Analysis" / "figures" / "lc_comparison_sne.pdf"
    generate_figure(out_pdf)


if __name__ == "__main__":
    main()
