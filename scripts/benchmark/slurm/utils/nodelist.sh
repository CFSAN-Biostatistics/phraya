#!/bin/bash
# Query sinfo for available nodes and output comma-separated nodelist
# Usage: nodelist.sh [num_nodes]

set -euo pipefail

NUM_NODES="${1:-3}"  # Default: 3 nodes for 3 replicates
PARTITION="${SLURM_PARTITION:-batch}"  # Default partition

# Query sinfo for nodes in IDLE or MIX state
AVAILABLE_NODES=$(sinfo -p "$PARTITION" -t idle,mix -h -o "%N" | head -n "$NUM_NODES")

if [[ -z "$AVAILABLE_NODES" ]]; then
    echo "ERROR: No available nodes in partition $PARTITION" >&2
    exit 1
fi

# Expand node ranges (e.g., "node[01-03]" → "node01,node02,node03")
NODELIST=$(scontrol show hostnames "$AVAILABLE_NODES" | paste -sd,)

if [[ -z "$NODELIST" ]]; then
    echo "ERROR: Failed to expand nodelist from: $AVAILABLE_NODES" >&2
    exit 1
fi

echo "$NODELIST"
