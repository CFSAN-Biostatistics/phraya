#!/bin/bash
# Phraya --strategy fast wrapper (seed-vote subsampling, ±150bp, divergence cutoff)
# Usage: phraya-fast.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/_phraya_run.sh" fast "$@"
