#!/usr/bin/env python3
"""Placement accuracy (PA) for long-read alignment, scored against ground truth.

Truth TSV (from gen_synthetic_long.py) columns:
    read_id  chrom  true_start  strand  true_end  n_subs  n_ins  n_del  indel_events
`true_start` is the leftmost reference coordinate (0-based) of the simulated read.

Two input kinds, chosen by file extension (no content sniffing, no silent fallback):
    *.queries        phraya cross-space queries sidecar: zstd(msgpack({read_id: [[space, pos, identity], ...]}))
    *.sam            SAM text (e.g. minimap2 -a); only primary records are scored
                     (flags 0x4 unmapped, 0x100 secondary, 0x800 supplementary are ignored)

Scoring, per truth read:
    mapped        >= 1 placement/record (on any reference)
    correct       the read's single *best* placement is on the true chrom and within
                  --tolerance bp of true_start.
                  phraya: best = highest-identity placement (earliest position on ties).
                  SAM:    the primary record; leading soft/hard clip is subtracted so the
                          position is where the read's first base would sit.
    correct_any   (phraya only; equals `correct` for SAM) at least one retained placement
                  is on the true chrom within tolerance. phraya retains many per-locus
                  placements, so this is reported separately and never used as PA.

    PA        = correct / n_truth       (unmapped reads count against the aligner)
    precision = correct / mapped

Output: one TSV line (see FIELDS), or `--header` to also print the field names.

Usage:
    long_read_accuracy.py <cross_space.phraya.queries|aln.sam> <reads.fq.truth.tsv> [--tolerance 100] [--header]
"""
import argparse
import csv
import io
import re
import sys

FIELDS = ["pa", "precision", "n_truth", "n_mapped", "n_correct", "n_correct_any"]
SAM_IGNORED_FLAGS = 0x4 | 0x100 | 0x800
CLIP_RE = re.compile(r"^(\d+)[SH]")


def chrom_to_label(chrom: str) -> str:
    """Mirror phraya-cli's sanitize_label_component: alnum/-/_/. pass through, else '_'."""
    return re.sub(r"[^0-9A-Za-z._-]", "_", chrom)


def load_truth(path: str) -> dict:
    """{read_id: (chrom, true_start)}"""
    with open(path, newline="") as f:
        return {
            row["read_id"]: (row["chrom"], int(row["true_start"]))
            for row in csv.DictReader(f, delimiter="\t")
        }


def load_phraya_queries(path: str) -> dict:
    import msgpack
    import zstandard

    with open(path, "rb") as f:
        # Frames written by phraya carry no content size; stream-decode.
        raw = zstandard.ZstdDecompressor().stream_reader(f).read()
    return msgpack.unpackb(raw, raw=False)


def load_sam_primary(path: str) -> dict:
    """{read_id: (chrom, clip-adjusted 0-based start)} for primary mapped records."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"malformed SAM record: {line[:80]!r}")
            flag = int(fields[1])
            if flag & SAM_IGNORED_FLAGS or fields[2] == "*":
                continue
            clip = CLIP_RE.match(fields[5])
            out[fields[0]] = (fields[2], int(fields[3]) - 1 - (int(clip.group(1)) if clip else 0))
    return out


def score_phraya(queries: dict, truth: dict, tol: int):
    mapped = correct = correct_any = 0
    for rid, (chrom, start) in truth.items():
        placements = queries.get(rid)
        if not placements:
            continue
        mapped += 1
        label = chrom_to_label(chrom)
        space, pos, _ = max(placements, key=lambda p: (p[2], -p[1]))
        if space == label and abs(pos - start) <= tol:
            correct += 1
        if any(p[0] == label and abs(p[1] - start) <= tol for p in placements):
            correct_any += 1
    return mapped, correct, correct_any


def score_sam(placements: dict, truth: dict, tol: int):
    mapped = correct = 0
    for rid, (chrom, start) in truth.items():
        hit = placements.get(rid)
        if hit is None:
            continue
        mapped += 1
        if hit[0] == chrom and abs(hit[1] - start) <= tol:
            correct += 1
    return mapped, correct, correct


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("alignments", help="*.queries (phraya sidecar) or *.sam")
    ap.add_argument("truth", help="truth TSV from gen_synthetic_long.py")
    ap.add_argument("--tolerance", type=int, default=100, help="bp (default 100)")
    ap.add_argument("--header", action="store_true", help="print field names first")
    args = ap.parse_args()

    truth = load_truth(args.truth)
    if args.alignments.endswith(".queries"):
        mapped, correct, correct_any = score_phraya(load_phraya_queries(args.alignments), truth, args.tolerance)
    elif args.alignments.endswith(".sam"):
        mapped, correct, correct_any = score_sam(load_sam_primary(args.alignments), truth, args.tolerance)
    else:
        sys.exit(f"unsupported extension (want .queries or .sam): {args.alignments}")

    n = len(truth)
    pa = correct / n if n else 0.0
    precision = correct / mapped if mapped else 0.0
    if args.header:
        print("\t".join(FIELDS))
    print(f"{pa:.4f}\t{precision:.4f}\t{n}\t{mapped}\t{correct}\t{correct_any}")


if __name__ == "__main__":
    main()
