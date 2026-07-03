#!/usr/bin/env python3
"""
Compute placement accuracy (PA) from a SAM/BAM file with wgsim-simulated reads.

wgsim encodes true origin in the read name:
    <chrom>_<start>_<end>_<rng1>_<rng2>_<strand>/1  or /2

PA = fraction of mapped reads whose alignment start is within <tolerance> bp
of the true start encoded in the read name.

Usage:
    sam_accuracy.py <file.sam|file.bam> [--tolerance 10] [--samtools PATH]

Outputs tab-separated to stdout:
    pa  n_mapped  n_correct  n_total
"""
import argparse
import re
import subprocess
import sys

SAMTOOLS = "/nfs/software/apps/micromamba/1.5.8/envs/samtools-v1.20/bin/samtools"
# Matches both wgsim (chr_start_end_...) and dwgsim (chr_start_end_strand_...) read names.
# dwgsim uses start>end for reverse-strand reads; true left pos = min(start,end).
# "rand_..." names are dwgsim's randomly-placed reads — exclude from PA evaluation.
WGSIM_RE = re.compile(r"^(.+?)_(\d+)_(\d+)_")


def parse_true_pos(qname):
    name = qname.rsplit("/", 1)[0] if "/" in qname else qname
    # dwgsim random reads: name starts with "rand_"
    if name.startswith("rand_"):
        return None
    m = WGSIM_RE.match(name)
    if m:
        pos1, pos2 = int(m.group(2)), int(m.group(3))
        return m.group(1), min(pos1, pos2)  # leftmost genomic position
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sam_file")
    parser.add_argument("--tolerance", type=int, default=10)
    parser.add_argument("--samtools", default=SAMTOOLS)
    args = parser.parse_args()

    proc = subprocess.Popen(
        [args.samtools, "view", "-F4", args.sam_file],
        stdout=subprocess.PIPE, text=True,
    )

    n_mapped = 0
    n_parseable = 0
    n_correct = 0

    for line in proc.stdout:
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        qname = fields[0]
        flag = int(fields[1])
        # Skip secondary/supplementary
        if flag & 0x900:
            continue
        aln_pos = int(fields[3])  # 1-based leftmost mapping position
        n_mapped += 1

        parsed = parse_true_pos(qname)
        if parsed is None:
            continue  # rand_/unparseable: skip (can't evaluate)
        n_parseable += 1
        _, true_start = parsed
        # Both positions are 1-based leftmost; compare directly
        if abs(aln_pos - true_start) <= args.tolerance:
            n_correct += 1

    proc.wait()
    pa = n_correct / n_parseable if n_parseable > 0 else 0.0
    print(f"{pa:.4f}\t{n_mapped}\t{n_correct}\t{n_mapped}")


if __name__ == "__main__":
    main()
