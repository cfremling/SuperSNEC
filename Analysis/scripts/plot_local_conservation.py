#!/usr/bin/env python3
"""Plot local energy conservation at mass boundaries for progenitor runs.

For each mass boundary m_b, at each snapshot we compute:
  E_enclosed = sum(eps * dm + 0.5 * v^2 * dm)  for zones inside m_b
  cumul_rad  = integral of L(m_b) * dt          (radiation leaving)
  cumul_PdV  = integral of P * 4*pi*r^2 * v * dt  (PdV work leaving)
  cumul_Ni   = integral of Ni_energy_rate * sum(NDF*dm)_inside * dt

Conservation residual:
  R(t) = E(t) - E(t0) + cumul_rad(t) + cumul_PdV(t) - cumul_Ni(t)
  Fractional = |R(t)| / |E(t0)| * 100  [percent]

This should be near zero if energy is locally conserved.
Budget starts at t0 ~ 0.5 d (after shock passage through all boundaries).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_sedov_convergence import parse_xg

ROOT = Path(__file__).resolve().parents[2]
TMPBASE = ROOT / "tmp" / "conservation_comparison"
FIGDIR = ROOT / "Analysis" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

MSUN = 1.989e33


def load_col2(p):
    d = np.loadtxt(p)
    return d[:, 1] if d.ndim == 2 else d


def compute_local_budget(rp, m_boundary, t_start_d=0.5):
    """Compute cumulative local energy conservation residual.

    Returns (times_d, frac_residual_pct).
    """
    eps_s = parse_xg(rp / "eps.xg")
    vel_s = parse_xg(rp / "vel.xg")
    lum_s = parse_xg(rp / "lum.xg")
    press_s = parse_xg(rp / "press.xg")
    radius_s = parse_xg(rp / "radius.xg")
    ni_dep_s = parse_xg(rp / "Ni_deposit_function.xg")
    ni_lum = np.loadtxt(rp / "Ni_total_luminosity.dat")

    t_start = t_start_d * 86400.0
    times, E_enc, L_bnd, PdV_bnd, Ni_rate = [], [], [], [], []

    for i in range(len(eps_s)):
        t = eps_s[i][0]
        if t < t_start:
            continue
        mc = eps_s[i][1]   # current mass grid (changes after remap)
        ev = eps_s[i][2]
        vv = vel_s[i][2]
        lv = lum_s[i][2]
        pv = press_s[i][2]
        rv = radius_s[i][2]
        nv = ni_dep_s[i][2]
        nz = len(ev)

        # --- Derive delta_mass from current snapshot mass grid ---
        dm_cur = np.diff(mc)
        dm_cur = np.append(dm_cur, dm_cur[-1])

        # --- Cell-center velocities (vel is at edges, eps at centers) ---
        v_c = np.zeros(nz)
        v_c[:nz - 1] = 0.5 * (vv[:nz - 1] + vv[1:nz])
        v_c[nz - 1] = vv[nz - 1]

        # --- Enclosed energy: sum full zones below boundary,
        #     plus fractional contribution of the straddling zone ---
        # mc[i] is the inner edge of zone i; zone i spans [mc[i], mc[i]+dm[i]]
        # Find the last zone fully inside the boundary
        nb = int(np.searchsorted(mc, m_boundary)) - 1
        nb = max(0, min(nb, nz - 2))

        # Full zones [0..nb-1]
        e_int = float(np.sum(ev[:nb] * dm_cur[:nb]))
        e_kin = float(np.sum(0.5 * v_c[:nb] ** 2 * dm_cur[:nb]))
        ni_enc = float(np.sum(nv[:nb] * dm_cur[:nb]))

        # Fractional zone nb: boundary falls within this zone
        if dm_cur[nb] > 0:
            frac_zone = min(1.0, max(0.0,
                (m_boundary - mc[nb]) / dm_cur[nb]))
        else:
            frac_zone = 0.0
        e_int += ev[nb] * dm_cur[nb] * frac_zone
        e_kin += 0.5 * v_c[nb] ** 2 * dm_cur[nb] * frac_zone
        ni_enc += nv[nb] * dm_cur[nb] * frac_zone

        # --- Boundary flux quantities: interpolate to exact m_boundary ---
        l_b = float(np.interp(m_boundary, mc, lv))
        p_b = float(np.interp(m_boundary, mc, pv))
        r_b = float(np.interp(m_boundary, mc, rv))
        v_b = float(np.interp(m_boundary, mc, vv))
        pdv = p_b * 4.0 * np.pi * r_b ** 2 * v_b

        # --- Ni heating inside boundary ---
        ni_total = float(np.sum(nv * dm_cur[:len(nv)]))
        frac_in = ni_enc / ni_total if ni_total > 0 else 0.0
        ni_lum_t = np.interp(t, ni_lum[:, 0], ni_lum[:, 1])

        times.append(t)
        E_enc.append(e_int + e_kin)
        L_bnd.append(l_b)
        PdV_bnd.append(pdv)
        Ni_rate.append(ni_lum_t * frac_in)

    times = np.array(times)
    E_enc = np.array(E_enc)
    L_bnd = np.array(L_bnd)
    PdV_bnd = np.array(PdV_bnd)
    Ni_rate = np.array(Ni_rate)

    E0 = E_enc[0]

    # Per-step residual: for each consecutive pair of snapshots,
    # the energy change should equal net flux + sources.
    # Accumulate |step_error| so the result is monotonically increasing.
    cumul_abs_err = np.zeros(len(times))
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        dE = E_enc[i] - E_enc[i - 1]
        L_avg = 0.5 * (L_bnd[i] + L_bnd[i - 1])
        PdV_avg = 0.5 * (PdV_bnd[i] + PdV_bnd[i - 1])
        Ni_avg = 0.5 * (Ni_rate[i] + Ni_rate[i - 1])
        step_err = abs(dE + (L_avg + PdV_avg - Ni_avg) * dt)
        cumul_abs_err[i] = cumul_abs_err[i - 1] + step_err

    frac_pct = cumul_abs_err / np.abs(E0) * 100.0

    return times / 86400.0, frac_pct


def main():
    plt.rcParams.update({"font.size": 12})

    # Find boundaries from a reference run
    # Use whichever SNEC 1000z run exists
    for ref_name in ["snec_orig_1000z", "snec_orig_1000z_adaptive"]:
        ref_p = TMPBASE / ref_name / "Data"
        if (ref_p / "He_init_frac.dat").exists():
            break

    mass_init = load_col2(ref_p / "mass_initial.dat")
    he_frac = np.loadtxt(ref_p / "He_init_frac.dat")[:, 1]
    idx_cohe = int(np.argmax(he_frac > 0.5))
    # Offset boundary by half a zone width so it falls MID-zone
    # in all grids, avoiding artifacts from landing on a zone edge
    m_cohe = mass_init[idx_cohe] + 0.5 * (mass_init[idx_cohe] - mass_init[idx_cohe - 1])

    mp = np.loadtxt(ref_p / "mass_photo.dat")
    m_photo_raw = mp[min(10, len(mp) - 1), 1]
    # Same offset for photosphere
    idx_ph = int(np.searchsorted(mass_init, m_photo_raw))
    m_photo = m_photo_raw + 0.5 * (mass_init[min(idx_ph, len(mass_init)-1)] - mass_init[max(idx_ph-1, 0)])

    runs = [
        ("snec_orig_1000z", "SNEC / SuperSNEC fixed, N=1000", "#1f77b4", "-"),
        ("snec_orig_100z", "SNEC / SuperSNEC fixed, N=100", "#2ca02c", "-"),
        ("supersnec_baseline_100z", "SuperSNEC adaptive, N=100", "#d62728", "--"),
        ("supersnec_baseline_200z", "SuperSNEC adaptive, N=200", "#9467bd", "--"),
        ("supersnec_baseline_500z", "SuperSNEC adaptive, N=500", "#8c564b", "--"),
        ("supersnec_baseline_1000z", "SuperSNEC adaptive, N=1000", "#ff7f0e", "--"),
    ]

    boundaries = [
        ("CO/He core boundary", m_cohe, f"m = {m_cohe / MSUN:.2f} $M_\\odot$"),
        ("Photosphere shell", m_photo, f"m = {m_photo / MSUN:.2f} $M_\\odot$"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(10, 5.2), sharey=True)

    for ax, (bname, m_b, mlabel) in zip(axes, boundaries):
        for rname, label, color, ls in runs:
            rp = TMPBASE / rname / "Data"
            if not (rp / "eps.xg").exists():
                continue
            t_d, frac = compute_local_budget(rp, m_b)
            ax.semilogy(t_d, frac, color=color, ls=ls, lw=1.4,
                        label=label, alpha=0.9)

        ax.set_xlabel("Time (days)")
        ax.set_title(f"{bname} ({mlabel})", fontsize=12)
        ax.set_xlim(0, 60)
        ax.set_ylim(None, 5.0)
        ax.minorticks_on()

    axes[0].set_ylabel("Energy conservation residual (%)")
    axes[0].legend(fontsize=8.5, frameon=False, loc="upper left")

    fig.tight_layout()
    outpath = FIGDIR / "local_energy_conservation.pdf"
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {outpath.name}")


if __name__ == "__main__":
    main()
