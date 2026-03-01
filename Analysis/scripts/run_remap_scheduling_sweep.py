#!/usr/bin/env python3
"""Sweep of adaptive remap scheduling modes (grid_adaptive_interval).

Runs 4 cases at 100 zones in adaptive_runtime mode and generates Table 7.

Usage:
  python3 Analysis/scripts/run_remap_scheduling_sweep.py
  python3 Analysis/scripts/run_remap_scheduling_sweep.py --tables-only
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
class SchedulingCase:
    case: str
    grid_adaptive_interval: str
    grid_update_interval_days: str
    description: str


# Sweep all 4 modes at multiple base intervals to show where
# physics-driven modes maintain accuracy while fixed degrades.
INTERVALS = ["0.50d0", "1.00d0", "2.00d0", "3.00d0"]
MODES = [
    ("0", "Fixed"),
    ("-1", "Photo+cap"),
    ("-2", "Photo"),
    ("-3", "Early-dense"),
]


def _itag(interval: str) -> str:
    return interval.replace("d0", "").replace(".", "p")


CASES: List[SchedulingCase] = [
    SchedulingCase(
        f"mode_{m.replace('-','m')}_{_itag(iv)}d",
        m, iv, f"{desc}, {iv.replace('d0','')}d base",
    )
    for iv in INTERVALS
    for m, desc in MODES
]


def build_remap_scheduling_table(rows: List[Dict[str, str]]) -> str:
    """Build a remap scheduling modes comparison table, grouped by base interval."""
    mode_order = {"0": 0, "-1": 1, "-2": 2, "-3": 3}
    mode_labels = {
        "0": "Fixed",
        "-1": "Photo+cap",
        "-2": "Photo",
        "-3": "Early-dense",
    }

    ok_rows = [r for r in rows if r.get("status", "") == "ok"]

    def _sort_key(r: Dict[str, str]) -> tuple:
        iv = fortran_to_float(r.get("grid_update_interval_days", "1.0d0"))
        mo = mode_order.get(r.get("grid_adaptive_interval", "0"), 99)
        return (iv, mo)

    ok_rows.sort(key=_sort_key)

    lines = [
        r"\begin{deluxetable*}{ccccccccc}",
        r"\tabletypesize{\scriptsize}",
        r"\tablecaption{Remap scheduling modes in \texttt{adaptive\_runtime} at 100 zones, "
        r"tested across four base intervals. "
        r"Mode 0 uses the base interval as a fixed remap period. "
        r"Mode $-1$: photosphere-driven with $t/20$ safety cap. "
        r"Mode $-2$: photosphere-driven (max 5\,d). "
        r"Mode $-3$: early-dense (0.25\,d for $t<3$\,d, 0.5\,d for 3--10\,d, then $t/20$). "
        r"All runs use the Ni-cadence and solver baselines from Sect.~\ref{sec:baseline}."
        r"\label{tab:remap_scheduling}}",
        r"\tablehead{",
        r"\colhead{Base (d)} &",
        r"\colhead{Mode} &",
        r"\colhead{Description} &",
        r"\colhead{Remaps} &",
        r"\colhead{Ni evals} &",
        r"\colhead{Runtime (s)} &",
        r"\colhead{RMS 5--30\,d} &",
        r"\colhead{RMS $>30$\,d} &",
        r"\colhead{RMS all}",
        r"}",
        r"\startdata",
    ]

    prev_iv = None
    for r in ok_rows:
        mode = r.get("grid_adaptive_interval", "0")
        desc = mode_labels.get(mode, "?")
        iv_raw = r.get("grid_update_interval_days", "1.0d0")
        iv = fortran_to_float(iv_raw)
        iv_str = f"{iv:.1f}" if iv == int(iv) else f"{iv:.2f}"
        remaps = int(float(r.get("remap_count", "0") or 0))
        ni_evals = int(float(r.get("ni_eval_total", "0") or 0))

        if prev_iv is not None and iv != prev_iv:
            lines.append(r"\hline")
        prev_iv = iv

        lines.append(
            f"{iv_str} & "
            f"{mode} & "
            f"{desc} & "
            f"{remaps} & "
            f"{ni_evals} & "
            f"{to_float(r['real_s']):.2f} & "
            f"{to_float(r['rms_5_30']):.4f} & "
            f"{to_float(r['rms_gt30']):.4f} & "
            f"{to_float(r['rms_all']):.4f} \\\\"
        )

    lines.extend(
        [
            r"\enddata",
            r"\tablecomments{At the 1\,d baseline interval all four modes perform comparably. "
            r"At longer intervals (2--3\,d), where the fixed mode degrades, "
            r"the photosphere-driven modes ($-1$, $-2$) converge to the same behavior as mode~0 "
            r"because the base interval dominates. "
            r"Mode $-3$ (early-dense) retains its dense early-time schedule regardless of the "
            r"base interval, maintaining lower RMS at 2--3\,d intervals.}",
            r"\end{deluxetable*}",
        ]
    )
    return "\n".join(lines) + "\n"


def generate_tables(results_tsv: Path) -> None:
    """Generate Table 7 from results.tsv."""
    TABLE_OUTDIR.mkdir(parents=True, exist_ok=True)

    if not results_tsv.exists():
        print(f"WARNING: {results_tsv} not found, skipping Table 7.")
        return

    rows = read_tsv(results_tsv)
    write_tex(
        TABLE_OUTDIR / "table_remap_scheduling.tex",
        build_remap_scheduling_table(rows),
    )
    print(f"  Wrote {TABLE_OUTDIR / 'table_remap_scheduling.tex'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", type=Path, default=DEFAULT_REF)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "Analysis" / "results" / "remap_scheduling",
    )
    parser.add_argument("--timeout-s", type=float, default=45.0)
    parser.add_argument("--nruns", type=int, default=3, help="Runs per case for median timing")
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

    fixed_overrides = {
        "outdir": '"Data"',
        "grid_mode": '"adaptive_runtime"',
        "imax": "100",
        "Ni_period": "5.0d4",
        "Ni_period_max": "5.456d5",
        "Ni_fractional_change": "0.70d0",
        "grid_debug": "0",
    }

    fields = [
        "case",
        "grid_adaptive_interval",
        "grid_update_interval_days",
        "description",
        "status",
        "real_s",
        "remap_count",
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

    nruns = args.nruns
    write_tsv(out_root / "cases.tsv", [asdict(c) for c in CASES], list(asdict(CASES[0]).keys()))

    for i, c in enumerate(CASES, start=1):
        overrides = dict(fixed_overrides)
        overrides["grid_adaptive_interval"] = c.grid_adaptive_interval
        overrides["grid_update_interval_days"] = c.grid_update_interval_days
        write_parameters_file(PARAMETERS, overrides)

        live_outdir = ROOT / "Data"
        runtimes: List[float] = []
        best_rc = 0
        best_timed_out = False

        for run_idx in range(nruns):
            if live_outdir.exists():
                shutil.rmtree(live_outdir)
            live_outdir.mkdir(parents=True, exist_ok=True)

            log_path = out_root / f"{c.case}_r{run_idx}.log"
            t0 = time.perf_counter()
            rc, timed_out = run_cmd(["./snec"], ROOT, log_path, timeout_s=args.timeout_s)
            runtime = time.perf_counter() - t0
            runtimes.append(runtime)
            if run_idx == 0:
                best_rc = rc
                best_timed_out = timed_out

        # Use median runtime
        runtimes.sort()
        median_runtime = runtimes[len(runtimes) // 2]
        median_idx = runtimes.index(median_runtime)

        # Re-run the median case to get its log and output (unless it was the last run)
        # Actually, the last run's output is still in live_outdir; use that for metrics
        # but report the median runtime. Re-read the median log for indicators.
        median_log = out_root / f"{c.case}_r{median_idx}.log"
        log_text = median_log.read_text(errors="replace")

        row: Dict[str, object] = {
            "case": c.case,
            "grid_adaptive_interval": c.grid_adaptive_interval,
            "grid_update_interval_days": c.grid_update_interval_days,
            "description": c.description,
            "status": "ok",
            "real_s": median_runtime,
            "remap_count": 0.0,
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

        row.update(parse_log_indicators(log_text))

        case_dir = run_root / c.case
        case_dir.mkdir(parents=True, exist_ok=True)
        (case_dir / "parameters_used").write_text(PARAMETERS.read_text())
        shutil.copy2(median_log, case_dir / "run.log")

        lum_path = live_outdir / "lum_observed.dat"
        if lum_path.exists():
            shutil.copy2(lum_path, case_dir / "lum_observed.dat")
        row["lum_nan_count"] = float(count_nan_in_lum(lum_path))

        if best_timed_out:
            row["status"] = "timeout"
        elif best_rc != 0:
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
        rt_str = ", ".join(f"{r:.2f}" for r in runtimes)
        print(
            f"[{i}/{len(CASES)}] {c.case}: status={row['status']} "
            f"runtimes=[{rt_str}] median={median_runtime:.2f}s "
            f"remaps={row['remap_count']:.0f} ni_evals={row['ni_eval_total']:.0f} "
            f"RMS(5-30,all)=({row['rms_5_30']},{row['rms_all']})"
        )

    print(f"\nWrote {out_root / 'results.tsv'} ({len(rows)} cases, {nruns} runs each)")

    # Generate table (always runs after sweep)
    generate_tables(out_root / "results.tsv")


if __name__ == "__main__":
    main()
