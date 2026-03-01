#!/usr/bin/env python3
"""Rebuild analysis artifacts and the paper PDF in one command.

Each sweep script is self-contained: it runs SNEC, writes TSV results, AND
generates its own LaTeX table(s).  With --skip-sweeps, every script is called
with --tables-only so tables are regenerated from cached data without re-running
SNEC.

Phase 1 (sweeps + tables):
  run_imax_scaling_scan.py             -> Tables 1, 5
  run_cumulative_speedup.py            -> Table 2
  run_adaptive_threshold_limits.py     -> Tables 3, 4
  run_ni_cadence_grid_mode_comparison.py -> Table 6
  run_remap_scheduling_sweep.py        -> Table 7
  run_legacy_ni_cadence_2d_sweep.py    -> Table 8

Phase 2 (figures):
  plot_adaptive_runtime_validation.py  -> Figure 1
  plot_lc_comparison_sne.py            -> Figure 2

Phase 3 (optional PDF):
  latexmk ...
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "Analysis" / "scripts"
PAPER_TEX = ROOT / "Paper" / "adaptive_gridding_for_snec.tex"
MAC_TEXBIN = Path("/Library/TeX/texbin")

# Each entry: (script, supports --tables-only)
SWEEP_TABLE_SCRIPTS = [
    (ANALYSIS / "run_imax_scaling_scan.py", True),
    (ANALYSIS / "run_cumulative_speedup.py", False),
    (ANALYSIS / "run_adaptive_threshold_limits.py", True),
    (ANALYSIS / "run_ni_cadence_grid_mode_comparison.py", True),
    (ANALYSIS / "run_remap_scheduling_sweep.py", True),
    (ANALYSIS / "run_legacy_ni_cadence_2d_sweep.py", True),
]

SOLVER_REPRO_SCRIPT = ANALYSIS / "reproduce_solver_analysis.py"


def run(cmd: list[str], cwd: Path = ROOT, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True, env=env)


def find_latexmk(explicit: str | None) -> str | None:
    if explicit:
        return explicit
    found = shutil.which("latexmk")
    if found:
        return found
    fallback = MAC_TEXBIN / "latexmk"
    if fallback.exists():
        return str(fallback)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-sweeps",
        action="store_true",
        help=(
            "Skip SNEC runs; call each sweep script with --tables-only "
            "to regenerate tables from cached results."
        ),
    )
    parser.add_argument(
        "--run-solver-analysis",
        action="store_true",
        help=(
            "Run the full solver-threshold reproduction pipeline "
            "(grid runs + repeats + solver figure/tables)."
        ),
    )
    parser.add_argument(
        "--skip-pdf",
        action="store_true",
        help="Do not build Paper/adaptive_gridding_for_snec.pdf.",
    )
    parser.add_argument(
        "--run-lc-cases",
        action="store_true",
        help="Re-run representative LC comparison cases before plotting overlays.",
    )
    parser.add_argument(
        "--latexmk",
        default=None,
        help="Path to latexmk executable (auto-detected by default).",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter used to run analysis scripts.",
    )
    args = parser.parse_args()

    # --- Phase 1: Sweeps + tables ---
    for script, has_tables_only in SWEEP_TABLE_SCRIPTS:
        if args.skip_sweeps:
            if has_tables_only:
                run([args.python, str(script), "--tables-only"])
            # Scripts without --tables-only are skipped entirely in skip mode
        else:
            run([args.python, str(script)])

    if not args.skip_sweeps and args.run_solver_analysis:
        if SOLVER_REPRO_SCRIPT.exists():
            run([args.python, str(SOLVER_REPRO_SCRIPT)])

    # --- Phase 2: Figures ---
    run([args.python, str(ANALYSIS / "plot_adaptive_runtime_validation.py")])
    run([args.python, str(ANALYSIS / "plot_lc_comparison_sne.py")])

    # --- Phase 3: PDF ---
    if not args.skip_pdf:
        latexmk = find_latexmk(args.latexmk)
        if not latexmk:
            print("latexmk not found. Install MacTeX (or make latexmk available) and rerun.")
            return 2

        env = os.environ.copy()
        if MAC_TEXBIN.exists():
            env["PATH"] = f"{MAC_TEXBIN}:{env.get('PATH', '')}"

        run(
            [
                latexmk,
                "-pdf",
                "-interaction=nonstopmode",
                "-file-line-error",
                PAPER_TEX.name,
            ],
            cwd=PAPER_TEX.parent,
            env=env,
        )
        print(f"Built {PAPER_TEX.with_suffix('.pdf')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
