#!/bin/bash
# Phraya --strategy sensitive --gap-model linear wrapper (ADR-0014).
#
# Once ADR-0014 ships, `sensitive` defaults to gap-affine (see
# phraya-sensitive.sh, unchanged — omitting --gap-model already gets the new
# default). This variant exists to keep the old linear-cost WFA path directly
# comparable in the same run: it's the exact-equivalence-with-Myers guarantee
# some users may still want, and the "before" side of any linear-vs-affine
# perf/IEC comparison.
#
# Usage: phraya-sensitive-linear.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/_phraya_run.sh" sensitive "$1" "$2" "$3" "$4" "$5" linear
