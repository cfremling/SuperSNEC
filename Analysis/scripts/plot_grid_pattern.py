#!/usr/bin/env python3
"""Plot grid pattern visualization: zone width vs mass coordinate at key epochs.

Shows how the adaptive grid redistributes zones to track the photosphere,
compared to the fixed GridPattern grid.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "tmp" / "conservation_comparison"
MSUN = 1.989e33
DAY = 86400.0


def parse_xg(path: Path) -> list[tuple[float, np.ndarray, np.ndarray]]:
    """Parse SNEC .xg file. Returns list of (time, col1, col2)."""
    snapshots = []
    current_time = None
    x_vals, y_vals = [], []

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if '"Time' in line:
            if current_time is not None and x_vals:
                snapshots.append((current_time, np.array(x_vals), np.array(y_vals)))
            current_time = float(line.split("=")[1])
            x_vals, y_vals = [], []
        else:
            parts = line.split()
            if len(parts) >= 2:
                x_vals.append(float(parts[0]))
                y_vals.append(float(parts[1]))

    if current_time is not None and x_vals:
        snapshots.append((current_time, np.array(x_vals), np.array(y_vals)))
    return snapshots


def find_snapshot(snapshots, target_time_s):
    """Find snapshot closest to target time."""
    times = [s[0] for s in snapshots]
    idx = np.argmin(np.abs(np.array(times) - target_time_s))
    return snapshots[idx]


def get_photosphere_mass(photo_path: Path, target_time_s: float) -> float:
    """Get photosphere mass coordinate at given time."""
    data = np.loadtxt(photo_path)
    idx = np.argmin(np.abs(data[:, 0] - target_time_s))
    return data[idx, 1]


def main():
    # Load adaptive and fixed grid mass profiles
    adaptive_mass = parse_xg(BASE / "supersnec_baseline_100z" / "Data" / "mass.xg")
    fixed_mass = parse_xg(BASE / "supersnec_compat_100z" / "Data" / "mass.xg")

    # Load photosphere mass data
    photo_m_adaptive = np.loadtxt(
        BASE / "supersnec_baseline_100z" / "Data" / "mass_photo.dat"
    )
    photo_m_fixed = np.loadtxt(
        BASE / "supersnec_compat_100z" / "Data" / "mass_photo.dat"
    )

    # Target epochs
    epochs = [
        (0.0, "Initial ($t = 0$)"),
        (10.0 * DAY, "$t \\approx 10\\,$d"),
        (25.0 * DAY, "$t \\approx 25\\,$d"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0), sharey=True)

    for col, (t_target, label) in enumerate(epochs):
        ax_m = axes[col]

        # --- Mass coordinate (top row) ---
        t_a, r_a, m_a = find_snapshot(adaptive_mass, t_target)
        t_f, r_f, m_f = find_snapshot(fixed_mass, t_target)

        dm_a = np.diff(m_a) / MSUN
        mc_a = 0.5 * (m_a[:-1] + m_a[1:]) / MSUN
        dm_f = np.diff(m_f) / MSUN
        mc_f = 0.5 * (m_f[:-1] + m_f[1:]) / MSUN

        ax_m.step(mc_f, dm_f, where="mid", color="C3", alpha=0.7,
                  lw=1.2, label="Fixed grid" if col == 1 else None)
        ax_m.step(mc_a, dm_a, where="mid", color="C0", alpha=0.9,
                  lw=1.5, label="Adaptive grid" if col == 1 else None)

        # Add zone boundary tick marks along the top of the panel
        for mi in m_a / MSUN:
            ax_m.axvline(mi, ymin=0.92, ymax=1.0, color="C0", lw=0.4, alpha=0.6)
        for mi in m_f / MSUN:
            ax_m.axvline(mi, ymin=0.0, ymax=0.08, color="C3", lw=0.4, alpha=0.6)

        if t_target > 0:
            idx = np.argmin(np.abs(photo_m_adaptive[:, 0] - t_target))
            ax_m.axvline(photo_m_adaptive[idx, 1] / MSUN, color="C0",
                         ls="--", lw=1.0, alpha=0.8, label="Photosphere" if col == 1 else None)
            idx = np.argmin(np.abs(photo_m_fixed[:, 0] - t_target))
            ax_m.axvline(photo_m_fixed[idx, 1] / MSUN, color="C3",
                         ls="--", lw=1.0, alpha=0.8)

        ax_m.set_yscale("log")
        ax_m.set_ylim(1e-4, 0.5)
        ax_m.set_xlabel("Mass coordinate ($M_\\odot$)")
        ax_m.set_title(label, fontsize=11)
        if col == 0:
            ax_m.set_ylabel("Zone width $\\Delta m$ ($M_\\odot$)")

    axes[1].legend(loc="lower left", fontsize=8, frameon=False,
                    bbox_to_anchor=(0.0, 0.10))

    plt.tight_layout()
    outpath = ROOT / "Analysis" / "figures" / "grid_pattern_visualization.pdf"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Saved: {outpath}")

    outpath_png = outpath.with_suffix(".png")
    fig.savefig(outpath_png, dpi=150, bbox_inches="tight")
    print(f"Saved: {outpath_png}")


if __name__ == "__main__":
    main()
