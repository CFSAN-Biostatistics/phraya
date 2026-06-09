#!/bin/bash
#SBATCH --job-name=stream_triad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0:30:00

# STREAM Triad platform characterization
# Measures sustained main-memory bandwidth for BNT normalization
# Usage: sbatch stream.slurm <run_dir>

set -euo pipefail

RUN_DIR=$1

if [[ -z "$RUN_DIR" ]]; then
    echo "ERROR: RUN_DIR not specified" >&2
    exit 1
fi

# Download and build STREAM if not present
STREAM_DIR="$HOME/.cache/stream"
mkdir -p "$STREAM_DIR"
cd "$STREAM_DIR"

if [[ ! -f stream ]]; then
    echo "Downloading STREAM benchmark..."
    wget -q https://www.cs.virginia.edu/stream/FTP/Code/stream.c

    echo "Compiling STREAM..."
    gcc -O3 -march=native -fopenmp -DSTREAM_ARRAY_SIZE=400000000 -o stream stream.c

    if [[ ! -f stream ]]; then
        echo "ERROR: STREAM compilation failed" >&2
        exit 1
    fi
fi

# Run STREAM, capture output
echo "Running STREAM Triad on node $SLURM_NODEID..."
./stream | tee "$RUN_DIR/stream_${SLURM_NODEID}.txt"

# Extract Triad GB/s (line format: "Triad:    80234.5   ...")
TRIAD=$(grep "Triad:" "$RUN_DIR/stream_${SLURM_NODEID}.txt" | awk '{print $2}')

if [[ -z "$TRIAD" ]]; then
    echo "ERROR: Failed to extract Triad result" >&2
    exit 1
fi

echo "STREAM Triad result: $TRIAD GB/s"

# Append to shared results file (node_id triad_gbps)
echo "$SLURM_NODEID $TRIAD" >> "$RUN_DIR/stream_triad.txt"

echo "STREAM characterization complete for node $SLURM_NODEID"
