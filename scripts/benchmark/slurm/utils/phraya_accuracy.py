#!/usr/bin/env python3
"""
Compute placement accuracy (PA) for phraya from .phraya.queries + wgsim read names.

wgsim/dwgsim encode the true FRAGMENT span in the read name:
    <chrom>_<start>_<end>_<frag_info>/1  or /2

Two simulator conventions coexist in this benchmark's data (see sam_accuracy.py):
wgsim encodes the fragment span (mate /2 aligns at max - read_len + 1), dwgsim encodes the
two mates' own leftmost positions directly. Scoring every read against min() alone counts
every right-end mate as misplaced, halving PA into a ~0.49 artifact unrelated to the aligner.
So a read is correct if its best alignment start is within <tolerance> bp of ANY candidate
true start: min, max, or (max - read_len + 1). This matches sam_accuracy.py exactly, so
phraya and the SAM aligners are scored on the same basis.

PA = (correctly placed reads) / (placed, parseable, non-random reads).

Also reports unaligned fraction if total_reads is supplied.

Usage:
    phraya_accuracy.py <file.phraya.queries> [--tolerance 10] [--read-len 150] [--total-reads N]

Outputs tab-separated to stdout:
    pa  n_mapped  n_correct  n_unaligned  unaligned_frac
"""
import argparse
import re
import sys
import zstandard
import msgpack


# Matches wgsim and dwgsim read names. dwgsim uses start>end for rev-strand reads.
WGSIM_RE = re.compile(r"^(.+?)_(\d+)_(\d+)_")


def parse_fragment(read_name: str):
    """Return (chrom, frag_lo, frag_hi) from a wgsim/dwgsim read name, or None."""
    # Strip /1 /2 suffix if present
    name = read_name.rsplit("/", 1)[0] if "/" in read_name else read_name
    if name.startswith("rand_"):  # dwgsim random/unplaceable reads
        return None
    m = WGSIM_RE.match(name)
    if m:
        pos1, pos2 = int(m.group(2)), int(m.group(3))
        return m.group(1), min(pos1, pos2), max(pos1, pos2)
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("queries_file")
    parser.add_argument("--tolerance", type=int, default=10)
    parser.add_argument("--read-len", type=int, default=150)
    parser.add_argument("--total-reads", type=int, default=0)
    args = parser.parse_args()

    dctx = zstandard.ZstdDecompressor()
    with open(args.queries_file, "rb") as f:
        with dctx.stream_reader(f) as reader:
            raw = reader.read()
    # QueryIndex = HashMap<String, Vec<(u32, f64)>>
    # msgpack decodes tuples as lists: [[pos, score], ...]
    data: dict = msgpack.unpackb(raw, raw=False)

    n_mapped = len(data)
    n_parseable = 0
    n_correct = 0

    for read_name, alignments in data.items():
        parsed = parse_fragment(read_name)
        if parsed is None:
            # rand_ reads or unparseable names — skip entirely (can't evaluate)
            continue
        if not alignments:
            continue  # no positions recorded — count as unmapped, skip
        n_parseable += 1
        _, frag_lo, frag_hi = parsed
        # Best alignment = highest score (second element)
        best_pos = int(max(alignments, key=lambda x: x[1])[0])
        # Correct if near any candidate true start (covers wgsim span + dwgsim direct):
        # frag_lo, frag_hi, or the wgsim right-end read start (frag_hi - read_len + 1).
        candidates = (frag_lo, frag_hi, frag_hi - args.read_len + 1)
        if any(abs(best_pos - c) <= args.tolerance for c in candidates):
            n_correct += 1

    pa = n_correct / n_parseable if n_parseable > 0 else 0.0

    total = args.total_reads if args.total_reads > 0 else n_mapped
    n_unaligned = max(0, total - n_mapped)
    unaligned_frac = n_unaligned / total if total > 0 else 0.0

    print(f"{pa:.4f}\t{n_mapped}\t{n_correct}\t{n_unaligned}\t{unaligned_frac:.4f}")


if __name__ == "__main__":
    main()
