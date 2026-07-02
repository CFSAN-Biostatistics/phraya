#!/bin/bash
# Phraya --strategy balanced wrapper (Myers fitting ≤500bp, WFA fallback, ±50bp)
# Usage: phraya.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/_phraya_run.sh" balanced "$@"
