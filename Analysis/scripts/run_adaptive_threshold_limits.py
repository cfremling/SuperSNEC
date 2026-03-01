#!/usr/bin/env python3
"""EPSTOL/ITMAX solver-threshold scan in adaptive_runtime mode.

Uses runtime parameters (no Fortran recompilation needed).
Generates Tables 3 & 4 (ITMAX-only and EPSTOL-only scans).

Usage:
  python3 Analysis/scripts/run_adaptive_threshold_limits.py
  python3 Analysis/scripts/run_adaptive_threshold_limits.py --tables-only
"""

from __future__ import annotations

import re
import shutil
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
import argparse
from typing import Dict, List

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
PARAMETERS = ROOT / "parameters"
DEFAULT_REF = ROOT / "tmp" / "snec_original_1000z"

# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------
sys.path.insert(0, str(ROOT / "Analysis" / "scripts"))
from snec_metrics import (  # noqa: E402
    compute_mag_metrics,
    count_nan_in_lum,
    parse_log_indicators,
    read_tsv,
    run_cmd,
    tex_escape,
    to_float,
    write_parameters_file,
    write_tex,
    write_tsv,
)

TABLE_OUTDIR = ROOT / "Paper" / "generated"


# ---------------------------------------------------------------------------
# Case definitions (from original run_adaptive_threshold_limits.py)
# ---------------------------------------------------------------------------

@dataclass
class ThresholdCase:
    case: str
    h_epstol: str
    h_itmax: str
    hr_epstol: str
    hr_itmax: str


def make_cases() -> List[ThresholdCase]:
    cases: List[ThresholdCase] = []

    # Baseline (tight): EPSTOL=1e-7, ITMAX=(100,300)
    cases.append(ThresholdCase("thrlim_base", "1.0d-7", "100", "1.0d-7", "300"))
    cases.append(ThresholdCase("thrlim_strict10x", "1.0d-8", "100", "1.0d-8", "300"))

    # EPSTOL relaxation sweep at fixed ITMAX=(100,300)
    for val, tag in [
        ("3.0d-7", "3e-7"),
        ("1.0d-6", "1e-6"),
        ("2.0d-6", "2e-6"),
        ("5.0d-6", "5e-6"),
        ("1.0d-5", "1e-5"),
        ("2.0d-5", "2e-5"),
        ("5.0d-5", "5e-5"),
        ("1.0d-4", "1e-4"),
        ("2.0d-4", "2e-4"),
        ("5.0d-4", "5e-4"),
    ]:
        cases.append(ThresholdCase(f"thrlim_eps_{tag}", val, "100", val, "300"))

    # ITMAX sweep at baseline EPSTOL=1e-7
    for h_it, hr_it in [
        ("80", "240"),
        ("60", "200"),
        ("40", "150"),
        ("30", "120"),
        ("20", "90"),
        ("15", "60"),
        ("10", "45"),
        ("8", "30"),
        ("6", "20"),
    ]:
        cases.append(ThresholdCase(f"thrlim_it_{h_it}_{hr_it}", "1.0d-7", h_it, "1.0d-7", hr_it))

    # Aggressive combined settings
    for eps, h_it, hr_it, tag in [
        ("1.0d-6", "40", "150", "a"),
        ("2.0d-6", "30", "120", "b"),
        ("5.0d-6", "20", "90", "c"),
        ("1.0d-5", "15", "60", "d"),
        ("2.0d-5", "10", "45", "e"),
        ("5.0d-5", "8", "30", "f"),
        ("1.0d-4", "6", "20", "g"),
    ]:
        cases.append(ThresholdCase(f"thrlim_combo_{tag}", eps, h_it, eps, hr_it))

    # Extended EPSTOL (probes degradation regime)
    for val, tag in [
        ("1.0d-3", "1e-3"),
        ("2.0d-3", "2e-3"),
    ]:
        cases.append(ThresholdCase(f"thrlim_eps_{tag}", val, "100", val, "300"))

    # Extreme ITMAX reduction
    for h_it, hr_it in [
        ("5", "15"),
        ("4", "12"),
        ("3", "9"),
        ("2", "6"),
    ]:
        cases.append(ThresholdCase(f"thrlim_it_{h_it}_{hr_it}", "1.0d-7", h_it, "1.0d-7", hr_it))

    return cases


# ---------------------------------------------------------------------------
# Table builders (Tables 3 & 4)
# ---------------------------------------------------------------------------

def _safe_status_diag(row: Dict[str, str]) -> str:
    status = str(row.get("status", "ok"))
    if status == "ok":
        return "ok"
    saha = int(float(row.get("saha_warnings", "0") or 0))
    lnan = int(float(row.get("lum_nan_count", "0") or 0))
    return f"{tex_escape(status)} (saha={saha}, NaN={lnan})"


def build_solver_itmax_only_table(rows_limits: List[Dict[str, str]]) -> str:
    rows = [
        r
        for r in rows_limits
        if r.get("status", "ok") == "ok" and (r["case"] == "thrlim_base" or r["case"].startswith("thrlim_it_"))
    ]
    if not rows:
        raise KeyError("No ITMAX-only rows found in threshold limits input.")

    base = next((r for r in rows if r["case"] == "thrlim_base"), None)
    if base is None:
        raise KeyError("Missing thrlim_base in ITMAX-only table input.")

    def _sort_key(r: Dict[str, str]) -> tuple[int, int]:
        if r["case"] == "thrlim_base":
            return (0, 10_000)
        return (1, -int(float(r["h_itmax"])))

    rows = sorted(rows, key=_sort_key)

    lines = [
        r"\begin{deluxetable*}{cccccc}",
        r"\tablecaption{ITMAX-only scan at fixed \texttt{EPSTOL}$=10^{-7}$ (adaptive-runtime setup).\label{tab:solver_itmax_only}}",
        r"\tablehead{",
        r"\colhead{$\mathrm{ITMAX}_{\mathrm{hyd}}$} &",
        r"\colhead{$\mathrm{ITMAX}_{\mathrm{rad}}$} &",
        r"\colhead{Runtime (s)} &",
        r"\colhead{Scratch} &",
        r"\colhead{RMS 5--30 d} &",
        r"\colhead{RMS $>30$ d}",
        r"}",
        r"\startdata",
    ]

    for r in rows:
        runtime = to_float(r["real_s"])
        lines.append(
            f"{int(float(r['h_itmax']))} & {int(float(r['hr_itmax']))} & "
            f"{runtime:.2f} & "
            f"{int(float(r.get('scratch_steps', '0') or 0))} & "
            f"{to_float(r['rms_5_30']):.4f} & {to_float(r['rms_gt30']):.4f} \\\\"
        )

    lines.extend(
        [
            r"\enddata",
            r"\tablecomments{Lowering \texttt{ITMAX} increases retry cost and does not produce reliable runtime gains in this fixed-\texttt{EPSTOL} test.}",
            r"\end{deluxetable*}",
        ]
    )
    return "\n".join(lines) + "\n"


def build_solver_epstol_only_table(
    rows_all: List[Dict[str, str]],
) -> str:
    rows = [
        r
        for r in rows_all
        if r["case"] == "thrlim_base" or r["case"].startswith("thrlim_eps_")
    ]
    if not rows:
        raise KeyError("No EPSTOL sweep rows found in threshold inputs.")

    def _eps_value(r: Dict[str, str]) -> float:
        return float(r["h_epstol"].replace("d", "e"))

    def _sort_key(r: Dict[str, str]) -> tuple[int, float]:
        if r["case"] == "thrlim_base":
            return (0, _eps_value(r))
        return (1, _eps_value(r))

    rows = sorted(rows, key=_sort_key)

    lines = [
        r"\begin{deluxetable*}{cccccc}",
        r"\tablecaption{EPSTOL sweep at fixed \texttt{ITMAX}$(\mathrm{hyd},\mathrm{rad})=(100,300)$ (adaptive-runtime setup).\label{tab:solver_epstol_only}}",
        r"\tablehead{",
        r"\colhead{$\mathrm{EPSTOL}$} &",
        r"\colhead{Runtime (s)} &",
        r"\colhead{RMS 5--30 d} &",
        r"\colhead{RMS $>30$ d} &",
        r"\colhead{Scratch} &",
        r"\colhead{Status / diagnostics}",
        r"}",
        r"\startdata",
    ]

    for r in rows:
        status = str(r.get("status", "ok"))
        rms_5_30 = "--" if status != "ok" else f"{to_float(r['rms_5_30']):.4f}"
        rms_gt30 = "--" if status != "ok" else f"{to_float(r['rms_gt30']):.4f}"
        lines.append(
            f"{tex_escape(r['h_epstol'])} & "
            f"{to_float(r['real_s']):.2f} & {rms_5_30} & {rms_gt30} & "
            f"{int(float(r.get('scratch_steps', '0') or 0))} & {_safe_status_diag(r)} \\\\"
        )

    lines.extend(
        [
            r"\enddata",
            r"\tablecomments{Fixed \texttt{ITMAX}$(\mathrm{hyd},\mathrm{rad})=(100,300)$. The dominant runtime lever is \texttt{EPSTOL}; values around $10^{-4}$ are fast, while $\gtrsim2\times10^{-4}$ gradually degrades quality (specially at $5-30$~d) or fails.}",
            r"\end{deluxetable*}",
        ]
    )
    return "\n".join(lines) + "\n"


def generate_tables(*results_tsvs: Path) -> None:
    """Generate Tables 3 & 4 from one or more results.tsv files."""
    TABLE_OUTDIR.mkdir(parents=True, exist_ok=True)

    all_rows: List[Dict[str, str]] = []
    for p in results_tsvs:
        if p.exists():
            all_rows.extend(read_tsv(p))
        else:
            print(f"WARNING: {p} not found, skipping.")

    if not all_rows:
        print("WARNING: no result rows found, skipping solver tables.")
        return

    write_tex(
        TABLE_OUTDIR / "table_solver_itmax_only.tex",
        build_solver_itmax_only_table(all_rows),
    )
    print(f"  Wrote {TABLE_OUTDIR / 'table_solver_itmax_only.tex'}")

    write_tex(
        TABLE_OUTDIR / "table_solver_epstol_only.tex",
        build_solver_epstol_only_table(all_rows),
    )
    print(f"  Wrote {TABLE_OUTDIR / 'table_solver_epstol_only.tex'}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", type=Path, default=DEFAULT_REF)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "Analysis" / "results" / "solver_threshold",
    )
    parser.add_argument("--case-filter", default="",
                        help="Optional regex filter on case labels.")
    parser.add_argument("--timeout-s", type=float, default=120.0,
                        help="Per-case timeout in seconds.")
    parser.add_argument(
        "--tables-only",
        action="store_true",
        help="Skip SNEC runs; regenerate tables from cached results.",
    )
    args = parser.parse_args()

    out_root = args.output_root

    if not args.tables_only:
        out_root.mkdir(parents=True, exist_ok=True)
        run_root = out_root / "runs"
        run_root.mkdir(parents=True, exist_ok=True)

        cases = make_cases()
        if args.case_filter:
            pat = re.compile(args.case_filter)
            cases = [c for c in cases if pat.search(c.case)]
            if not cases:
                raise ValueError(f"No cases matched --case-filter={args.case_filter!r}")
        write_tsv(out_root / "cases.tsv", [asdict(c) for c in cases], list(asdict(cases[0]).keys()))

        rows: List[Dict[str, object]] = []
        result_fields = [
            "case", "h_epstol", "h_itmax", "hr_epstol", "hr_itmax",
            "status", "real_s",
            "eos_problems", "scratch_steps", "saha_warnings",
            "lum_nan_count", "reached_tend",
            "rms_0_5", "rms_5_30", "rms_gt30", "rms_gt5", "rms_all",
            "max_abs_mag", "n_pass",
        ]

        live_outdir = ROOT / "Data"

        for i, c in enumerate(cases, 1):
            print(f"[{i}/{len(cases)}] {c.case}: epstol={c.h_epstol} itmax=({c.h_itmax},{c.hr_itmax})")

            row: Dict[str, object] = {
                "case": c.case,
                "h_epstol": c.h_epstol,
                "h_itmax": c.h_itmax,
                "hr_epstol": c.hr_epstol,
                "hr_itmax": c.hr_itmax,
                "status": "pending",
                "real_s": float("nan"),
                "eos_problems": float("nan"),
                "scratch_steps": float("nan"),
                "saha_warnings": float("nan"),
                "lum_nan_count": float("nan"),
                "reached_tend": 0.0,
                "rms_0_5": float("nan"),
                "rms_5_30": float("nan"),
                "rms_gt30": float("nan"),
                "rms_gt5": float("nan"),
                "rms_all": float("nan"),
                "max_abs_mag": float("nan"),
                "n_pass": float("nan"),
            }

            # Write self-contained parameters file with per-case solver settings.
            write_parameters_file(PARAMETERS, {
                "outdir": '"Data"',
                "epstol_hydro": c.h_epstol,
                "epstol_rad": c.hr_epstol,
                "itmax_hydro": c.h_itmax,
                "itmax_rad": c.hr_itmax,
            })

            t0 = time.perf_counter()
            run_rc, timed_out = run_cmd(
                ["./snec"], ROOT, out_root / f"{c.case}.log",
                timeout_s=args.timeout_s,
            )
            runtime = time.perf_counter() - t0
            row["real_s"] = runtime

            # Archive output directory.
            archived_outdir = run_root / f"Data_{c.case}"
            if archived_outdir.exists():
                shutil.rmtree(archived_outdir)
            if live_outdir.exists():
                shutil.copytree(live_outdir, archived_outdir)

            log_text = (out_root / f"{c.case}.log").read_text(errors="replace")
            row.update(parse_log_indicators(log_text))

            lum = live_outdir / "lum_observed.dat"
            row["lum_nan_count"] = float(count_nan_in_lum(lum))
            if timed_out:
                row["status"] = "timeout"
            elif run_rc != 0:
                row["status"] = "run_fail"
            elif not lum.exists():
                row["status"] = "no_output"
            else:
                try:
                    metrics = compute_mag_metrics(args.ref_dir / "lum_observed.dat", lum)
                    row.update(metrics)
                    row["status"] = "ok"
                except Exception:
                    row["status"] = "metric_fail"

            rows.append(row)
            write_tsv(out_root / "results.tsv", rows, result_fields)
            print(f"  -> status={row['status']} runtime={runtime:.1f}s")

        print(f"\nWrote {out_root / 'results.tsv'} ({len(rows)} cases)")

    # Generate tables (always runs)
    generate_tables(out_root / "results.tsv")


if __name__ == "__main__":
    main()
