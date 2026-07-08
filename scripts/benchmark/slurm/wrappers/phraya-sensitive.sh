#!/bin/bash
# Phraya --strategy sensitive wrapper (WFA all-anchors K=∞, ±25bp coverage)
# Usage: phraya-sensitive.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/_phraya_run.sh" sensitive "$@"
