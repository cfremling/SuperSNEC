#!/usr/bin/env bash
set -euo pipefail

# Usage (run from repo root):
#   benchmark/benchmark_snec.sh [MAX_JOBS]
#
# Example:
#   benchmark/benchmark_snec.sh 8
#
# Environment:
#   SNEC_BIN  (optional) path to snec binary, default: ./snec

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

MAX_JOBS="${1:-8}"
TEMPLATE_DIR="$SCRIPT_DIR/template_run"

# Resolve SNEC binary to an absolute path
SNEC_BIN_REL="${SNEC_BIN:-./snec}"
SNEC_BIN="$(cd "$(dirname "$SNEC_BIN_REL")" && pwd)/$(basename "$SNEC_BIN_REL")"

if [[ ! -x "$SNEC_BIN" ]]; then
  echo "ERROR: SNEC binary '$SNEC_BIN' not found or not executable." >&2
  exit 1
fi

if [[ ! -d "$TEMPLATE_DIR" ]]; then
  echo "ERROR: Template directory '$TEMPLATE_DIR' does not exist." >&2
  exit 1
fi

echo "=== SNEC parallel benchmark ==="
echo "SNEC_BIN     = $SNEC_BIN"
echo "TEMPLATE_DIR = $TEMPLATE_DIR"
echo "MAX_JOBS     = $MAX_JOBS"
echo

BENCH_DIR="$SCRIPT_DIR/runs"
mkdir -p "$BENCH_DIR"

# Clear old runs
rm -rf "$BENCH_DIR"/job_* 2>/dev/null || true

for njobs in $(seq 1 "$MAX_JOBS"); do
  echo "-------------------------------------------"
  echo "Testing with $njobs parallel SNEC run(s)..."

  # Create per-job directories with copies of the template
  for j in $(seq 1 "$njobs"); do
    RUN_DIR="$BENCH_DIR/job_${njobs}_$j"
    mkdir -p "$RUN_DIR"
    rsync -a "$TEMPLATE_DIR"/ "$RUN_DIR"/
  done

  # Use python for subsecond timing
  start=$(python3 -c "import time; print(f'{time.time():.3f}')")

  pids=()
  for j in $(seq 1 "$njobs"); do
    RUN_DIR="$BENCH_DIR/job_${njobs}_$j"
    (
      cd "$RUN_DIR"
      "$SNEC_BIN" > snec.log 2>&1
    ) &
    pids+=("$!")
  done

  # Wait for all jobs in this batch
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
      echo "WARNING: a SNEC job (PID $pid) exited with non-zero status" >&2
    fi
  done

  end=$(python3 -c "import time; print(f'{time.time():.3f}')")
  elapsed=$(python3 -c "print(f'{$end - $start:.2f}')")

  echo "Finished $njobs run(s) in ${elapsed}s  (${elapsed}s wall)"
  throughput=$(python3 -c "print(f'{$njobs / ($end - $start):.3f}')")
  per_run=$(python3 -c "print(f'{($end - $start) / $njobs:.2f}')")
  echo "Throughput: ${throughput} runs/s  |  ${per_run}s per run"
  echo
done

echo "=== Benchmark complete ==="
echo "Optimal njobs: pick the row with highest throughput (runs/s)."
