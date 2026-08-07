#!/usr/bin/env python3
"""Compute Indel Event Concordance (IEC) — see BENCHMARK_EXPANSION.md.

For each simulated read/query with N true indel events (from a `gen_synthetic*.py`
`.truth.tsv` sidecar), a read is *concordant* if its reported alignment's CIGAR
contains exactly N indel (I/D) operations. IEC = concordant / (reads with N > 0).

This is deliberately an op-*count* comparison, not an op-*position* comparison:
`phraya filter --format tsv` emits one row per (variant position, supporting read)
with no explicit per-read alignment-start column, so reconstructing each indel
op's absolute reference position from that sidecar is more machinery than the
signal is worth. Op-count concordance needs only `provenance` (read/query ID,
joins directly to `truth.tsv`'s ID column — `phraya` sets provenance to the raw
query ID) and `cigar` (identical on every TSV row for a given read, since CIGAR
is a whole-read property), both already present in `phraya filter --format tsv`
output.

Usage:
    phraya filter run.phraya --format tsv | \\
        compute_indel_recovery.py --truth reads.fq.truth.tsv

    compute_indel_recovery.py --truth reads.fq.truth.tsv --tsv variants.tsv
"""
import argparse
import re
import sys
from collections import defaultdict

# phraya's CIGAR alphabet is M/X/I/D (see wfa_simd.rs); count I/D op *runs*, not
# lengths — one 5D run is one event, matching how truth.tsv records one indel event
# regardless of its length.
_CIGAR_OP = re.compile(r"(\d+)([MXID])")


def count_indel_ops(cigar: str) -> int:
    return sum(1 for _count, op in _CIGAR_OP.findall(cigar) if op in ("I", "D"))


def load_truth(path: str) -> dict[str, int]:
    """Returns {read_id: n_true_indel_events}."""
    truth = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        id_col = header.index("read_id") if "read_id" in header else header.index("query_id")
        events_col = header.index("indel_events")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            events = fields[events_col]
            n_events = 0 if not events else len(events.split(";"))
            truth[fields[id_col]] = n_events
    return truth


def load_cigars_by_provenance(fh) -> dict[str, str]:
    """Returns {read_id: cigar} from `phraya filter --format tsv` output
    (header: position, ref_base, all_alleles, mapq, confidence, cigar,
    edit_distance, coverage, avg_base_quality, provenance — see
    `phraya-cli/src/main.rs::output_tsv`). All rows sharing a read carry the
    same CIGAR (a whole-read property), so the first row wins.

    `provenance` is `query.id()`, which for FASTQ input is the *whole* header
    line after `@` (e.g. `read_0 pos=17358` for a gen_synthetic.py read, not
    just `read_0`) — split on the first whitespace and key on the leading
    token, matching `truth.tsv`'s bare `read_id`/`query_id` column.
    """
    header = fh.readline().rstrip("\n").split("\t")
    cigar_col = header.index("cigar")
    prov_col = header.index("provenance")
    cigars: dict[str, str] = {}
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        if len(fields) <= max(cigar_col, prov_col):
            continue
        read_id = fields[prov_col].split(" ", 1)[0]
        cigars.setdefault(read_id, fields[cigar_col])
    return cigars


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True, help="gen_synthetic*.py .truth.tsv sidecar")
    ap.add_argument("--tsv", help="phraya filter --format tsv output (default: stdin)")
    args = ap.parse_args()

    truth = load_truth(args.truth)
    with (open(args.tsv) if args.tsv else sys.stdin) as fh:
        cigars = load_cigars_by_provenance(fh)

    reads_with_indels = {rid: n for rid, n in truth.items() if n > 0}
    concordant = 0
    unplaced = 0
    mismatched: list[tuple[str, int, int]] = []

    for read_id, n_true in reads_with_indels.items():
        cigar = cigars.get(read_id)
        if cigar is None:
            unplaced += 1
            continue
        n_reported = count_indel_ops(cigar)
        if n_reported == n_true:
            concordant += 1
        else:
            mismatched.append((read_id, n_true, n_reported))

    total = len(reads_with_indels)
    iec = concordant / total if total > 0 else float("nan")

    print(f"reads_with_indel_events\t{total}")
    print(f"concordant\t{concordant}")
    print(f"unplaced\t{unplaced}")
    print(f"discordant\t{len(mismatched)}")
    print(f"iec\t{iec:.4f}" if total > 0 else "iec\tnan")

    if mismatched:
        print("\n# discordant reads (read_id, n_true_events, n_reported_ops)", file=sys.stderr)
        for read_id, n_true, n_reported in mismatched[:20]:
            print(f"#   {read_id}\t{n_true}\t{n_reported}", file=sys.stderr)
        if len(mismatched) > 20:
            print(f"#   ... and {len(mismatched) - 20} more", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
