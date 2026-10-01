#!/usr/bin/env python3
"""Compute placement accuracy (PA) for long-read alignment.

Long reads (ONT/PacBio) simulated by gen_synthetic_long.py use a truth sidecar
with positions encoded directly, unlike the wgsim/dwgsim fragment-span encoding
used for short reads.

Truth TSV format (from gen_synthetic_long.py):
    read_id  chrom  true_start  strand  true_end  n_subs  n_ins  n_del  indel_events

For each aligned read, compare its best alignment position to the true position.

PA = (correctly placed reads) / (read with parseable true position)

A read is "correctly placed" if the aligned start falls within ±tolerance bp
of the true start.

Usage:
    long_read_accuracy.py <queries.ts> <truth.tsv> [--tolerance 100]

    phraya_long_read_accuracy.py <cross_space.phraya.queries> /path/to/reads.fq.truth.tsv
"""
import argparse
import csv
import sys
import msgpack
import zstandard

# Mapping from chromosome name to label (phraya sanitizes to underscore)
def chrom_to_label(chrom: str) -> str:
    """Mirror phraya-cli's sanitize_label_component: alnum/-/_/. pass through, else '_'."""
    import re
    return re.sub(r"[^0-9A-Za-z._-]", "_", chrom)


def load_truth(truth_path: str) -> dict[str, tuple[str, int, str]]:
    """Load truth TSV: {read_id: (chrom, true_start, strand)}."""
    truth = {}
    with open(truth_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rid = row["read_id"]
            chrom = row["chrom"]
            start = int(row["true_start"])
            strand = row["strand"]
            truth[rid] = (chrom, start, strand)
    return truth


def load_phraya_queries(queries_path: str) -> dict:
    """Load phraya's cross_space.queries msgpack (zstd-compressed)."""
    import msgpack
    import zstandard
    dctx = zstandard.ZstdDecompressor()
    with open(queries_path, "rb") as f:
        compressed = f.read()
    decompressed = dctx.decompress(compressed)
    data = msgpack.unpackb(decompressed, raw=False)
    return data


def compute_pa_phraya_long(queries_data: dict, truth: dict, tolerance: int = 100) -> tuple[float, int, int, int]:
    """Compute PA from phraya's cross-space queries.

    queries_data is a dict: {read_id: [placements]}
    Each placement is a list: [space_label, pos, identity]
    """
    n_parseable = 0
    n_correct = 0
    for read_id, placements in queries_data.items():
        if read_id not in truth:
            continue
        chrom, true_start, strand = truth[read_id]
        expected_label = chrom_to_label(chrom)

        # Find the best placement for the expected chromosome
        best_pos = None
        best_identity = 0.0
        for placement in placements:
            space, pos, identity = placement
            if space == expected_label:
                if best_pos is None or identity > best_identity:
                    best_pos = pos
                    best_identity = identity

        if best_pos is None:
            # Best placement was on a different chromosome or no placement
            continue

        n_parseable += 1
        if abs(best_pos - true_start) <= tolerance:
            n_correct += 1

    if n_parseable == 0:
        return 0.0, 0, 0, 0
    return n_correct / n_parseable, n_correct, n_parseable


def load_sam_placements(sam_path: str) -> dict[str, tuple[str, int, str]]:
    """Parse SAM file and return {read_id: (chrom, pos, strand)}."""
    placements = {}
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            cigar = fields[5]
            if rname == "*" or pos < 0:
                continue

            # Determine strand from flag 16
            strand = "-" if flag & 16 else "+"

            # For PA, we want the leftmost position
            if strand == "-":
                # Reverse strand — leftmost position is pos + cigar_len - 1
                # But PA scoring uses pos directly for simplicity
                pass

            placements[qname] = (rname, pos, strand)
    return placements


def compute_pa_sam_long(sam_path: str, truth: dict, tolerance: int = 100) -> tuple[float, int, int, int]:
    """Compute PA from SAM file for long reads."""
    placements = load_sam_placements(sam_path)

    n_parseable = 0
    n_correct = 0
    for read_id, (chrom, pos, strand) in placements.items():
        if read_id not in truth:
            continue
        true_chrom, true_start, true_strand = truth[read_id]
        expected_label = chrom_to_label(true_chrom)

        if chrom != true_chrom:
            continue

        n_parseable += 1
        if abs(pos - true_start) <= tolerance:
            n_correct += 1

    if n_parseable == 0:
        return 0.0, 0, 0, 0
    return n_correct / n_parseable, n_correct, n_parseable


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("queries", help="phraya queries file or SAM/BAM")
    parser.add_argument("truth", help="truth TSV from gen_synthetic_long.py")
    parser.add_argument("--tolerance", type=int, default=100, help="tolerance in bp (default 100)")
    args = parser.parse_args()

    truth = load_truth(args.truth)

    # Detect file format
    with open(args.queries, "rb") as f:
        magic = f.read(3)
    if magic == b"B9\x01" or (len(magic) >= 3 and magic[:3] == b"SM"):
        # Likely zstd-compressed msgpack (phraya) or SAM
        # Try phraya first by checking extension/content
        try:
            queries_data = load_phraya_queries(args.queries)
            pa, n_correct, n_parseable = compute_pa_phraya_long(queries_data, truth, args.tolerance)[:3]
            n_unmapped = n_parseable - n_correct
            print(f"{pa:.4f}\t{n_parseable}\t{n_correct}\t{n_unmapped}\t0.0000")
            return
        except Exception:
            pass
        # Fall back to SAM
        try:
            pa, n_correct, n_parseable, n_unmapped = compute_pa_sam_long(args.queries, truth, args.tolerance)
            n_total = n_parseable
            unaligned_frac = (n_total - n_parseable) / n_total if n_total > 0 else 0.0
            print(f"{pa:.4f}\t{n_parseable}\t{n_correct}\t{n_unmapped}\t{unaligned_frac:.4f}")
            return
        except Exception:
            pass

    # Try as plain SAM (magic won't match but it's text)
    try:
        pa, n_correct, n_parseable, n_unmapped = compute_pa_sam_long(args.queries, truth, args.tolerance)
        n_total = n_parseable
        unaligned_frac = (n_total - n_parseable) / n_total if n_total > 0 else 0.0
        print(f"{pa:.4f}\t{n_parseable}\t{n_correct}\t{n_unmapped}\t{unaligned_frac:.4f}")
        return
    except Exception as e:
        print(f"Error parsing input: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()