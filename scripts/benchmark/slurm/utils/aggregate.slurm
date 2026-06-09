#!/bin/bash
#SBATCH --job-name=benchmark_aggregate
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00

# Aggregation job: collect per-replicate results → score.py format
# Runs after benchmark array job completes

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -z "${RUN_DIR:-}" ]]; then
    echo "ERROR: RUN_DIR must be exported by orchestrator" >&2
    exit 1
fi

echo "Aggregating results from $RUN_DIR..."
python3 "$SCRIPT_DIR/utils/aggregate_results.py" "$RUN_DIR" > "$RUN_DIR/results.json"

if [[ -f "$RUN_DIR/results.json" ]]; then
    echo "Results aggregated to $RUN_DIR/results.json"
    echo
    echo "To score results:"
    echo "  python ~/data-commons/test/benchmarking/alignment/score.py $RUN_DIR/results.json"
else
    echo "ERROR: Failed to create results.json" >&2
    exit 1
fi
