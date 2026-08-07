#!/bin/bash
# Phraya protein-mode wrapper, --strategy balanced (ADR-0013)
# Usage: phraya-protein.sh <proteome.fasta> <query.fasta> <out_dir> <threads>
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/_phraya_run_protein.sh" balanced "$@"
