#!/usr/bin/env python3
"""2D legacy-pattern Ni-cadence sweep: f_Ni x Ni_period_max.

Varies Ni_fractional_change and Ni_period_max jointly in legacy_pattern
mode at 100 zones. Generates Table 8.

Usage:
  python3 Analysis/scripts/run_legacy_ni_cadence_2d_sweep.py
  python3 Analysis/scripts/run_legacy_ni_cadence_2d_sweep.py --tables-only
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
class NiCase2D:
    case: str
    ni_frac: str
    ni_period_max: str


def fortran_to_float(x: str) -> float:
    return float(x.lower().replace("d", "e"))


def tag_from_fortran(x: str) -> str:
    tag = x.lower().strip()
    tag = tag.replace("d-", "dm")
    tag = tag.replace("d+", "d")
    tag = tag.replace(".", "p")
    tag = tag.replace("-", "m")
    return tag


# 2D grid
FRAC_VALS = ["0.05d0", "0.10d0", "0.20d0", "0.50d0", "0.70d0", "1.0d0"]
PMAX_VALS = ["5.0d4", "1.728d5", "3.456d5", "5.456d5", "8.64d5", "1.728d6"]


def make_cases() -> List[NiCase2D]:
    cases: List[NiCase2D] = []
    for fni in FRAC_VALS:
        for pm in PMAX_VALS:
            name = f"leg_2d_f{tag_from_fortran(fni)}_pm{tag_from_fortran(pm)}"
            cases.append(NiCase2D(case=name, ni_frac=fni, ni_period_max=pm))
    return cases


def build_legacy_ni_cadence_table(rows: List[Dict[str, str]]) -> str:
    """Build a 2D legacy-pattern Ni cadence sweep table (f_Ni x Ni_period_max)."""
    ok_rows = [r for r in rows if r.get("status", "") == "ok"]
    ok_rows.sort(key=lambda r: (fortran_to_float(r["ni_fractional_change"]),
                                 fortran_to_float(r["ni_period_max"])))

    lines = [
        r"\begin{deluxetable*}{cccccccc}",
        r"\tabletypesize{\scriptsize}",
        r"\tablecaption{Two-parameter $^{56}$Ni-cadence sweep in \texttt{legacy\_pattern} mode at 100 zones. "
        r"\texttt{Ni\_fractional\_change} ($f_{\mathrm{Ni}}$) and \texttt{Ni\_period\_max} "
        r"($\Delta t_{\mathrm{max}}$) are varied jointly; "
        r"\texttt{Ni\_period}$=5\times10^4$~s is fixed. "
        r"Solver baseline: \texttt{EPSTOL}$=10^{-4}$, \texttt{ITMAX}$(100,300)$."
        r"\label{tab:legacy_ni_cadence}}",
        r"\tablehead{",
        r"\colhead{$f_{\mathrm{Ni}}$} &",
        r"\colhead{$\Delta t_{\mathrm{max}}$ (s)} &",
        r"\colhead{Runtime (s)} &",
        r"\colhead{Ni evals} &",
        r"\colhead{RMS 0--5\,d} &",
        r"\colhead{RMS 5--30\,d} &",
        r"\colhead{RMS $>30$\,d} &",
        r"\colhead{RMS all}",
        r"}",
        r"\startdata",
    ]

    prev_fni = None
    for r in ok_rows:
        fni = fortran_to_float(r["ni_fractional_change"])
        pmax = fortran_to_float(r["ni_period_max"])
        ni_evals = int(float(r.get("ni_eval_total", "0") or 0))

        fni_key = f"{fni:.2f}"
        if prev_fni is not None and fni_key != prev_fni:
            lines.append(r"\hline")
        prev_fni = fni_key

        if pmax >= 1e6:
            pmax_str = f"${pmax/1e6:.3g}\\times10^6$"
        elif pmax >= 1e5:
            pmax_str = f"${pmax/1e5:.3g}\\times10^5$"
        else:
            pmax_str = f"${pmax/1e4:.3g}\\times10^4$"

        lines.append(
            f"{fni:.2f} & "
            f"{pmax_str} & "
            f"{to_float(r['real_s']):.2f} & "
            f"{ni_evals} & "
            f"{to_float(r['rms_0_5']):.4f} & "
            f"{to_float(r['rms_5_30']):.4f} & "
            f"{to_float(r['rms_gt30']):.4f} & "
            f"{to_float(r['rms_all']):.4f} \\\\"
        )

    n_ok = len(ok_rows)
    n_total = len(rows)
    lines.extend(
        [
            r"\enddata",
            rf"\tablecomments{{{n_ok} of {n_total} cases complete successfully. "
            r"The fixed grid has no remap-forced Ni recalculations. "
            r"Increasing $\Delta t_{\mathrm{max}}$ reduces Ni evaluations and runtime, "
            r"but at aggressive $f_{\mathrm{Ni}}$ the late-time RMS degrades sharply. "
            r"Conservative settings ($f_{\mathrm{Ni}} \leq 0.10$) are required for sub-0.1\,mag accuracy.}",
            r"\end{deluxetable*}",
        ]
    )
    return "\n".join(lines) + "\n"


def generate_tables(results_tsv: Path) -> None:
    """Generate Table 8 from results.tsv."""
    TABLE_OUTDIR.mkdir(parents=True, exist_ok=True)

    if not results_tsv.exists():
        print(f"WARNING: {results_tsv} not found, skipping Table 8.")
        return

    rows = read_tsv(results_tsv)
    write_tex(
        TABLE_OUTDIR / "table_legacy_ni_cadence.tex",
        build_legacy_ni_cadence_table(rows),
    )
    print(f"  Wrote {TABLE_OUTDIR / 'table_legacy_ni_cadence.tex'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", type=Path, default=DEFAULT_REF)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "Analysis" / "results" / "legacy_ni_cadence_2d",
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

    fixed_overrides = {
        "grid_mode": '"legacy_pattern"',
        "grid_update_interval_days": "0.0d0",
        "imax": "100",
        "Ni_period": "5.0d4",
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
        "ni_fractional_change",
        "ni_period_max",
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
            case_overrides = dict(fixed_overrides)
            case_overrides["Ni_fractional_change"] = c.ni_frac
            case_overrides["Ni_period_max"] = c.ni_period_max
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
                "ni_fractional_change": c.ni_frac,
                "ni_period_max": c.ni_period_max,
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
                f"[{i:03d}/{len(cases):03d}] {c.case}: status={row['status']} "
                f"runtime={runtime:.2f}s ni_evals={row['ni_eval_total']:.0f} "
                f"RMS(5-30,>30)=({row['rms_5_30']},{row['rms_gt30']})"
            )

    print(f"\nWrote {out_root / 'results.tsv'} ({len(rows)} cases)")

    # Generate table (always runs after sweep)
    generate_tables(out_root / "results.tsv")


if __name__ == "__main__":
    main()
