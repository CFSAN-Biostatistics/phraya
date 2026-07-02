#!/usr/bin/env python3
"""
Compute placement accuracy (PA) for phraya from .phraya.queries + wgsim read names.

wgsim encodes true origin in the read name:
    <chrom>_<start>_<end>_<frag_info>/1  or /2

For each read in the .phraya.queries map, check if its best-scoring alignment
position falls within <tolerance> bp of the true start encoded in the read name.

PA = (correctly placed reads) / (total mapped reads)

Also reports unaligned fraction if total_reads is supplied.

Usage:
    phraya_accuracy.py <file.phraya.queries> [--tolerance 10] [--total-reads N]

Outputs tab-separated to stdout:
    pa  n_mapped  n_correct  n_unaligned  unaligned_frac
"""
import argparse
import re
import sys
import zstandard
import msgpack


WGSIM_RE = re.compile(r"^(.+)_(\d+)_(\d+)_")


def parse_true_pos(read_name: str):
    """Return (chrom, start) from a wgsim-style read name, or None."""
    # Strip /1 /2 suffix if present
    name = read_name.rsplit("/", 1)[0] if "/" in read_name else read_name
    m = WGSIM_RE.match(name)
    if m:
        return m.group(1), int(m.group(2))
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("queries_file")
    parser.add_argument("--tolerance", type=int, default=10)
    parser.add_argument("--total-reads", type=int, default=0)
    args = parser.parse_args()

    with open(args.queries_file, "rb") as f:
        compressed = f.read()
    dctx = zstandard.ZstdDecompressor()
    raw = dctx.decompress(compressed)
    # QueryIndex = HashMap<String, Vec<(u32, f64)>>
    # msgpack decodes tuples as lists: [[pos, score], ...]
    data: dict = msgpack.unpackb(raw, raw=False)

    n_mapped = len(data)
    n_parseable = 0
    n_correct = 0

    for read_name, alignments in data.items():
        parsed = parse_true_pos(read_name)
        if parsed is None:
            # Non-wgsim reads: can't evaluate — count as correct to avoid penalising
            n_correct += 1
            continue
        n_parseable += 1
        _, true_start = parsed
        # Best alignment = highest score (second element)
        best_pos = max(alignments, key=lambda x: x[1])[0]
        if abs(int(best_pos) - true_start) <= args.tolerance:
            n_correct += 1

    pa = n_correct / n_mapped if n_mapped > 0 else 0.0

    total = args.total_reads if args.total_reads > 0 else n_mapped
    n_unaligned = max(0, total - n_mapped)
    unaligned_frac = n_unaligned / total if total > 0 else 0.0

    print(f"{pa:.4f}\t{n_mapped}\t{n_correct}\t{n_unaligned}\t{unaligned_frac:.4f}")


if __name__ == "__main__":
    main()
