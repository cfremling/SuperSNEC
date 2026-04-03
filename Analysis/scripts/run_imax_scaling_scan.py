#!/usr/bin/env python3
"""Run imax scaling scan + quadrature scan and generate tables.

Runs SNEC at multiple zone counts with the current baseline parameters,
then runs a quadrature order scan at imax=100. Compares all light curves
against the original (unmodified) 1000-zone SNEC reference. Generates:
  - Paper/generated/table_imax_scaling.tex     (Table 1: zone-count scaling)
  - Paper/generated/table_quad_convergence.tex (Table 5: quadrature convergence)

Usage:
  python3 Analysis/scripts/run_imax_scaling_scan.py
  python3 Analysis/scripts/run_imax_scaling_scan.py --tables-only
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

ROOT = Path(__file__).resolve().parents[2]
TMP = ROOT / "tmp"
OUTDIR = ROOT / "Paper" / "generated"
PARAMS = ROOT / "parameters"

sys.path.insert(0, str(ROOT / "Analysis" / "scripts"))
from snec_metrics import write_parameters_file  # noqa: E402

# Original unmodified SNEC 1000-zone reference
SNEC_REF = TMP / "snec_original_1000z" / "lum_observed.dat"
SNEC_LEGACY_RUNTIME = 671  # seconds, measured externally

IMAX_VALUES = [60, 80, 100, 200, 500, 1000]
QUAD_VALUES = [30, 50, 70, 90, 110, 130, 150]
N_TIMING_RUNS = 3


def load_lc(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load luminosity file, return (time_days, mag)."""
    data = np.loadtxt(path)
    t_days = data[:, 0] / 86400.0
    with np.errstate(divide="ignore", invalid="ignore"):
        mag = -2.5 * np.log10(np.maximum(data[:, 1], 1e-90))
    return t_days, mag


def rms_in_range(
    t_ref: np.ndarray,
    m_ref: np.ndarray,
    t_test: np.ndarray,
    m_test: np.ndarray,
    t_lo: float,
    t_hi: float,
) -> float:
    """RMS magnitude residual in [t_lo, t_hi] days."""
    mask = (t_ref >= t_lo) & (t_ref <= t_hi)
    if not mask.any():
        return 0.0
    f = interp1d(t_test, m_test, bounds_error=False, fill_value=np.nan)
    m_interp = f(t_ref[mask])
    valid = ~np.isnan(m_interp)
    if not valid.any():
        return 0.0
    return float(np.sqrt(np.mean((m_interp[valid] - m_ref[mask][valid]) ** 2)))


def read_runtimes(path: Path) -> dict[int, float]:
    """Read a two-column (key, value) runtime TSV produced by a scan."""
    out: dict[int, float] = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        parts = line.strip().split()
        if len(parts) >= 2:
            out[int(parts[0])] = float(parts[1])
    return out


def median_runtime(n: int = N_TIMING_RUNS) -> float:
    """Run SNEC n times and return median wall-clock time."""
    times = []
    for _ in range(n):
        t0 = time.time()
        subprocess.run(
            ["./snec"],
            capture_output=True,
            timeout=1200,
            cwd=ROOT,
        )
        times.append(time.time() - t0)
    return float(np.median(times))


def run_imax_scan() -> dict[int, float]:
    """Run SNEC at each imax value, save LCs and runtimes."""
    runtimes: dict[int, float] = {}

    for imax in IMAX_VALUES:
        print(f"  imax={imax} ...", flush=True)
        outdir = TMP / f"imax_{imax}"
        outdir.mkdir(parents=True, exist_ok=True)

        write_parameters_file(PARAMS, {
            "imax": str(imax),
            "Ni_quad_npoints": "150",
            "smooth_ni_luminosity": "0",
            "outdir": f'"{outdir}/"',
            "grid_debug": "0",
        })

        # Clean stale diagnostic files
        for f in ROOT.glob("grid_diag_*.dat"):
            f.unlink()
        for f in ROOT.glob("fort.*"):
            f.unlink()

        t = median_runtime()
        runtimes[imax] = round(t, 1)
        print(f"    Runtime: {t:.1f}s")

        # Copy LC to tmp/
        lc_src = outdir / "lum_observed.dat"
        lc_dst = TMP / f"lum_baseline_imax{imax}.dat"
        if lc_src.exists():
            shutil.copy2(lc_src, lc_dst)

    # Write runtimes TSV
    tsv = TMP / "imax_scan_runtimes.tsv"
    with open(tsv, "w") as f:
        for imax in IMAX_VALUES:
            f.write(f"{imax}\t{runtimes[imax]}\n")

    return runtimes


def run_quad_scan() -> dict[int, float]:
    """Run quadrature scan at imax=100 with varying Ni_quad_npoints."""
    runtimes: dict[int, float] = {}

    for nq in QUAD_VALUES:
        print(f"  npoints={nq} ...", flush=True)
        outdir = TMP / f"quad_np{nq}"
        outdir.mkdir(parents=True, exist_ok=True)

        write_parameters_file(PARAMS, {
            "imax": "100",
            "Ni_quad_npoints": str(nq),
            "smooth_ni_luminosity": "0",
            "outdir": f'"{outdir}/"',
            "grid_debug": "0",
        })

        # Clean stale diagnostic files
        for f in ROOT.glob("grid_diag_*.dat"):
            f.unlink()
        for f in ROOT.glob("fort.*"):
            f.unlink()

        t = median_runtime()
        runtimes[nq] = round(t, 1)
        print(f"    Runtime: {t:.1f}s")

        # Copy LC
        lc_src = outdir / "lum_observed.dat"
        lc_dst = TMP / f"lum_quad_np{nq}.dat"
        if lc_src.exists():
            shutil.copy2(lc_src, lc_dst)

    # Write runtimes TSV
    tsv = TMP / "quad_scan_runtimes.tsv"
    with open(tsv, "w") as f:
        for nq in QUAD_VALUES:
            f.write(f"{nq}\t{runtimes[nq]}\n")

    return runtimes


def generate_tables(
    imax_runtimes: dict[int, float],
    quad_runtimes: dict[int, float],
) -> None:
    """Generate zone-count scaling and quadrature convergence tables."""
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Load reference LC
    t_ref, m_ref = load_lc(SNEC_REF)

    # ---- Table 1: Zone-count scaling ----
    lines: list[str] = []
    L = lines.append

    L(r"\begin{deluxetable}{rrrr}")
    L(r"\setlength{\tabcolsep}{6pt}")
    L(
        r"\tablecaption{Runtime and accuracy as a function of zone count $N$. "
        r"SuperSNEC rows use the \texttt{adaptive\_runtime} grid mode and "
        r"baseline settings (Sect.~\ref{sec:baseline}), with the SNEC default "
        r"gamma-ray deposition quadrature order "
        r"$n_\mathrm{quad} = 150$ (Sect.~\ref{sec:ni_quad}). "
        r"RMS is measured against "
        r"a 1000-zone reference run produced with the standard "
        r"(unmodified) SNEC code and default \texttt{GridPattern.dat}. "
        r"\label{tab:imax_scaling}}"
    )
    L(r"\tablehead{")
    L(
        r"\colhead{\hfill$N$} & \colhead{$t$ (s)} & "
        r"\colhead{$\Delta m_{5\text{--}30}$} & "
        r"\colhead{$\Delta m_\mathrm{all}$}"
    )
    L(r"}")
    L(r"\startdata")
    L(r"\cutinhead{SuperSNEC ($n_\mathrm{quad} = 150$)}")

    for imax in IMAX_VALUES:
        lc_path = TMP / f"lum_baseline_imax{imax}.dat"
        t, m = load_lc(lc_path)
        r1 = rms_in_range(t_ref, m_ref, t, m, 5, 30)
        r2 = rms_in_range(t_ref, m_ref, t, m, 0.5, 60)
        rt = imax_runtimes.get(imax, 0.0)
        L(f"{imax:4d} & {rt:5.1f} & {r1:.3f} & {r2:.3f} \\\\")

    L(r"\cutinhead{Reference (SNEC 1000-zone)}")
    L(f"1000 & {SNEC_LEGACY_RUNTIME} & \\nodata & \\nodata \\\\")

    L(r"\enddata")
    L(
        r"\tablecomments{$\Delta m_{5\text{--}30}$ and $\Delta m_\mathrm{all}$ "
        r"are RMS magnitude residuals in the 5--30\,d and $>$0.5\,d windows.}"
    )
    L(r"\end{deluxetable}")

    out_path = OUTDIR / "table_imax_scaling.tex"
    out_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out_path}")

    # ---- Table 5: Quadrature convergence ----
    lines = []
    L = lines.append

    L(r"\begin{deluxetable*}{rrrr}")
    L(r"\setlength{\tabcolsep}{12pt}")
    L(
        r"\tablecaption{Quadrature convergence at $N = 100$ zones. "
        r"$n_\mathrm{quad}$ is the number of integration points per "
        r"dimension in the gamma-ray deposition calculation. "
        r"All rows use the adaptive baseline (Sect.~\ref{sec:baseline}) "
        r"except for the varied $n_\mathrm{quad}$. "
        r"RMS is measured against a 1000-zone SNEC reference. "
        r"\label{tab:quad_convergence}}"
    )
    L(r"\tablehead{")
    L(
        r"\colhead{\hfill$n_\mathrm{quad}$} & \colhead{$t$ (s)} & "
        r"\colhead{$\Delta m_{5\text{--}30}$} & "
        r"\colhead{$\Delta m_\mathrm{all}$}"
    )
    L(r"}")
    L(r"\startdata")
    L(r"\cutinhead{SuperSNEC ($N = 100$)}")

    for nq in QUAD_VALUES:
        lc_path = TMP / f"lum_quad_np{nq}.dat"
        t, m = load_lc(lc_path)
        r1 = rms_in_range(t_ref, m_ref, t, m, 5, 30)
        r2 = rms_in_range(t_ref, m_ref, t, m, 0.5, 60)
        rt = quad_runtimes.get(nq, 0.0)
        L(f"{nq:4d} & {rt:5.1f} & {r1:.3f} & {r2:.3f} \\\\")

    L(r"\enddata")
    L(
        r"\tablecomments{$\Delta m_{5\text{--}30}$ and $\Delta m_\mathrm{all}$ "
        r"are RMS magnitude residuals in the 5--30\,d and $>$0.5\,d windows. "
        r"The RMS is essentially flat for $n_\mathrm{quad} \ge 50$.}"
    )
    L(r"\end{deluxetable*}")

    out_path = OUTDIR / "table_quad_convergence.tex"
    out_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tables-only",
        action="store_true",
        help="Skip SNEC runs; regenerate tables from cached tmp/ data.",
    )
    args = parser.parse_args()

    TMP.mkdir(parents=True, exist_ok=True)

    if not SNEC_REF.exists():
        print(f"ERROR: Reference LC not found: {SNEC_REF}", file=sys.stderr)
        sys.exit(1)

    if not args.tables_only:
        # Imax scan
        print("\n=== Imax scaling scan ===", flush=True)
        imax_runtimes = run_imax_scan()

        # Quadrature scan
        print("\n=== Quadrature scan ===", flush=True)
        quad_runtimes = run_quad_scan()
    else:
        # Load cached runtimes
        imax_runtimes = read_runtimes(TMP / "imax_scan_runtimes.tsv")
        quad_runtimes = read_runtimes(TMP / "quad_scan_runtimes.tsv")

    # Generate tables (always runs)
    print("\n=== Generating tables ===", flush=True)
    generate_tables(imax_runtimes, quad_runtimes)


if __name__ == "__main__":
    main()
