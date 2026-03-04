#!/usr/bin/env python3
"""Cross-mode Ni cadence comparison: matched settings across adaptive_runtime and legacy_pattern.

Runs the same Ni cadence configurations in both grid modes to directly compare
Ni evaluation counts and accuracy. Generates Table 6 (Ni cadence comparison).

Usage:
  python3 Analysis/scripts/run_ni_cadence_grid_mode_comparison.py
  python3 Analysis/scripts/run_ni_cadence_grid_mode_comparison.py --tables-only
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import argparse
import shutil
import sys
import time
from typing import Dict, List


ROOT = Path(__file__).resolve().parents[2]
PARAMETERS = ROOT / "parameters"
DEFAULT_REF = ROOT / "tmp" / "snec_original_1000z"

sys.path.insert(0, str(ROOT / "Analysis" / "scripts"))
from snec_metrics import (  # noqa: E402
    compute_mag_metrics,
    count_nan_in_lum,
    fortran_to_float,
    parse_log_indicators,
    read_tsv,
    run_cmd,
    to_float,
    write_parameters_file,
    write_tex,
    write_tsv,
)

TABLE_OUTDIR = ROOT / "Paper" / "generated"


@dataclass
class CompCase:
    case: str
    grid_mode: str
    grid_update_interval_days: str
    ni_period: str
    ni_period_max: str
    ni_frac: str


def make_cases() -> List[CompCase]:
    """Generate matched Ni cadence configs across both grid modes."""
    base_period = "5.0d4"
    base_period_max = "5.456d5"

    # Key f_Ni values to compare across grid modes
    frac_vals = ["0.05d0", "0.10d0", "0.20d0", "0.70d0", "1.0d0"]

    cases: List[CompCase] = []
    for f in frac_vals:
        ftag_clean = f.lower().strip().replace("d-", "dm").replace("d+", "d").replace(".", "p").replace("-", "m")

        # adaptive_runtime with 1.0d interval
        cases.append(CompCase(
            case=f"comp_adapt_f{ftag_clean}",
            grid_mode="adaptive_runtime",
            grid_update_interval_days="1.0d0",
            ni_period=base_period,
            ni_period_max=base_period_max,
            ni_frac=f,
        ))
        # legacy_pattern (no grid updates)
        cases.append(CompCase(
            case=f"comp_legacy_f{ftag_clean}",
            grid_mode="legacy_pattern",
            grid_update_interval_days="0.0d0",
            ni_period=base_period,
            ni_period_max=base_period_max,
            ni_frac=f,
        ))

    return cases


def build_ni_cadence_comparison_table(rows: List[Dict[str, str]]) -> str:
    """Build a cross-mode Ni cadence comparison table."""
    by_frac: Dict[str, Dict[str, Dict[str, str]]] = {}
    for r in rows:
        if r.get("status", "") != "ok":
            continue
        fni_str = r["ni_fractional_change"]
        fni = fortran_to_float(fni_str)
        fni_key = f"{fni:.2f}"
        mode = r["grid_mode"]
        if fni_key not in by_frac:
            by_frac[fni_key] = {}
        by_frac[fni_key][mode] = r

    frac_keys = sorted(by_frac.keys(), key=float)

    lines = [
        r"\begin{deluxetable*}{lccccccc}",
        r"\tablecaption{Ni-cadence tolerance comparison across grid modes at matched settings "
        r"($\Delta t_{\mathrm{min}} = 5\times10^4$~s, $\Delta t_{\mathrm{max}} = 5.456\times10^5$~s, "
        r"100 zones). Each row shows a paired \texttt{adaptive\_runtime} and \texttt{legacy\_pattern} "
        r"run at the same $f_{\mathrm{Ni}}$.\label{tab:ni_cadence_comparison}}",
        r"\tablehead{",
        r"\colhead{$f_{\mathrm{Ni}}$} &",
        r"\colhead{Grid mode} &",
        r"\colhead{Runtime (s)} &",
        r"\colhead{Ni evals} &",
        r"\colhead{RMS 0--5\,d} &",
        r"\colhead{RMS 5--30\,d} &",
        r"\colhead{RMS $>30$\,d} &",
        r"\colhead{RMS all}",
        r"}",
        r"\startdata",
    ]

    for fk in frac_keys:
        modes = by_frac[fk]
        for mode_key, mode_label in [("adaptive_runtime", "adaptive"), ("legacy_pattern", "legacy")]:
            if mode_key not in modes:
                continue
            r = modes[mode_key]
            ni_evals = int(float(r.get("ni_eval_total", "0") or 0))
            lines.append(
                f"{fk} & "
                rf"\texttt{{{mode_label}}} & "
                f"{to_float(r['real_s']):.2f} & "
                f"{ni_evals} & "
                f"{to_float(r['rms_0_5']):.4f} & "
                f"{to_float(r['rms_5_30']):.4f} & "
                f"{to_float(r['rms_gt30']):.4f} & "
                f"{to_float(r['rms_all']):.4f} \\\\"
            )
        lines.append(r"\hline")

    if lines and lines[-1] == r"\hline":
        lines.pop()

    lines.extend(
        [
            r"\enddata",
            r"\tablecomments{In \texttt{adaptive\_runtime} mode, grid remaps reset the Ni cadence timer, "
            r"producing supplementary Ni evaluations that are included in the Ni evals count. "
            r"The \texttt{legacy\_pattern} mode has no remap-forced updates. "
            r"This mechanism allows \texttt{adaptive\_runtime} to tolerate aggressive "
            r"$f_{\mathrm{Ni}}$ values that would degrade accuracy in \texttt{legacy\_pattern}.}",
            r"\end{deluxetable*}",
        ]
    )
    return "\n".join(lines) + "\n"


def generate_tables(results_tsv: Path) -> None:
    """Generate Table 6 from results.tsv."""
    TABLE_OUTDIR.mkdir(parents=True, exist_ok=True)

    if not results_tsv.exists():
        print(f"WARNING: {results_tsv} not found, skipping Table 6.")
        return

    rows = read_tsv(results_tsv)
    write_tex(
        TABLE_OUTDIR / "table_ni_cadence_comparison.tex",
        build_ni_cadence_comparison_table(rows),
    )
    print(f"  Wrote {TABLE_OUTDIR / 'table_ni_cadence_comparison.tex'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", type=Path, default=DEFAULT_REF)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "Analysis" / "results" / "ni_cadence_grid_mode_comparison",
    )
    parser.add_argument("--timeout-s", type=float, default=120.0)
    parser.add_argument(
        "--tables-only",
        action="store_true",
        help="Skip SNEC runs; regenerate table from cached results.",
    )
    args = parser.parse_args()

    out_root = args.output_root

    if args.tables_only:
        generate_tables(out_root / "results.tsv")
        return

    out_root.mkdir(parents=True, exist_ok=True)
    run_root = out_root / "runs"
    run_root.mkdir(parents=True, exist_ok=True)

    cases = make_cases()
    write_tsv(out_root / "cases.tsv", [asdict(c) for c in cases], list(asdict(cases[0]).keys()))

    common_overrides = {
        "imax": "100",
        "dtmax": "1.0d5",
        "dtout_fast": "5.0d4",
        "dtout_mid": "1.0d5",
        "dtout_scalar_fast": "1.0d4",
        "dtout_scalar_mid": "5.0d4",
        "ntinfo": "1000",
        "outdir": '"Data"',
        "grid_debug": "0",
    }

    fields = [
        "case",
        "grid_mode",
        "grid_update_interval_days",
        "ni_period",
        "ni_period_max",
        "ni_fractional_change",
        "status",
        "real_s",
        "ni_eval_total",
        "ni_eval_remap",
        "eos_problems",
        "scratch_steps",
        "saha_warnings",
        "lum_nan_count",
        "reached_tend",
        "rms_0_5",
        "rms_5_30",
        "rms_gt30",
        "rms_gt5",
        "rms_all",
        "max_abs_mag",
        "n_pass",
    ]
    rows: List[Dict[str, object]] = []

    live_outdir = ROOT / "Data"

    for i, c in enumerate(cases, start=1):
            case_overrides = {
                **common_overrides,
                "grid_mode": f'"{c.grid_mode}"',
                "grid_update_interval_days": c.grid_update_interval_days,
                "Ni_period": c.ni_period,
                "Ni_period_max": c.ni_period_max,
                "Ni_fractional_change": c.ni_frac,
            }
            write_parameters_file(PARAMETERS, case_overrides)

            if live_outdir.exists():
                shutil.rmtree(live_outdir)
            live_outdir.mkdir(parents=True, exist_ok=True)

            log_path = out_root / f"{c.case}.log"
            t0 = time.perf_counter()
            rc, timed_out = run_cmd(["./snec"], ROOT, log_path, timeout_s=args.timeout_s)
            runtime = time.perf_counter() - t0

            row: Dict[str, object] = {
                "case": c.case,
                "grid_mode": c.grid_mode,
                "grid_update_interval_days": c.grid_update_interval_days,
                "ni_period": c.ni_period,
                "ni_period_max": c.ni_period_max,
                "ni_fractional_change": c.ni_frac,
                "status": "ok",
                "real_s": runtime,
                "ni_eval_total": 0.0,
                "ni_eval_remap": 0.0,
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

            log_text = log_path.read_text(errors="replace")
            row.update(parse_log_indicators(log_text))

            run_case_dir = run_root / c.case
            run_case_dir.mkdir(parents=True, exist_ok=True)
            (run_case_dir / "parameters_used").write_text(PARAMETERS.read_text())
            shutil.copy2(log_path, run_case_dir / "run.log")

            lum_path = live_outdir / "lum_observed.dat"
            if lum_path.exists():
                shutil.copy2(lum_path, run_case_dir / "lum_observed.dat")

            row["lum_nan_count"] = float(count_nan_in_lum(lum_path))

            if timed_out:
                row["status"] = "timeout"
            elif rc != 0:
                row["status"] = "run_fail"
            elif not lum_path.exists():
                row["status"] = "no_output"
            else:
                try:
                    metrics = compute_mag_metrics(args.ref_dir / "lum_observed.dat", lum_path)
                    row.update(metrics)
                    if row["reached_tend"] < 0.5:
                        row["status"] = "no_tend"
                    elif row["lum_nan_count"] > 0:
                        row["status"] = "lum_nan"
                    else:
                        row["status"] = "ok"
                except Exception:
                    row["status"] = "metric_fail"

            rows.append(row)
            write_tsv(out_root / "results.tsv", rows, fields)
            print(
                f"[{i:03d}/{len(cases):03d}] {c.case}: "
                f"mode={c.grid_mode} status={row['status']} "
                f"runtime={runtime:.2f}s ni_evals={row['ni_eval_total']:.0f}({row['ni_eval_remap']:.0f}) "
                f"RMS_all={row['rms_all']}"
            )

    print(f"\nWrote {out_root / 'results.tsv'} ({len(rows)} cases)")

    # Generate table (always runs after sweep)
    generate_tables(out_root / "results.tsv")


if __name__ == "__main__":
    main()
