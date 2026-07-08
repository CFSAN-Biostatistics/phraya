#!/usr/bin/env python3
"""
Compute placement accuracy (PA) from a SAM/BAM file with wgsim-simulated reads.

The read name encodes the true origin, but TWO simulator conventions coexist in this
benchmark's data and must both be handled:

  * wgsim   — `<chrom>_<fragStart>_<fragEnd>_<rng1>_<rng2>_<strand>` where the two numbers
    are the FRAGMENT span. Mate /1 aligns at min; mate /2 aligns at its own leftmost base
    ~ max - read_len + 1.
  * dwgsim  — `<chrom>_<pos1>_<pos2>_<...>` where the two numbers are the two mates'
    OWN leftmost positions directly (reverse-strand fragments have pos1 > pos2).

Scoring every read against min() alone counts every right-end mate as misplaced, collapsing
PA to a ~0.49 artifact unrelated to the aligner. To cover both conventions mate-agnostically,
a read is correct if its aligned start is within <tolerance> bp of ANY of the candidate true
starts: min, max, or (max - read_len + 1). The fragment is far longer than the tolerance, so
these candidates never overlap and no spurious credit is given.

PA = fraction of evaluated (mapped, parseable, non-random) reads correctly placed.

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
# dwgsim uses start>end for reverse-strand reads; the fragment span is [min, max].
# "rand_..." names are dwgsim's randomly-placed reads — exclude from PA evaluation.
WGSIM_RE = re.compile(r"^(.+?)_(\d+)_(\d+)_")


def parse_fragment(qname):
    """Return (chrom, frag_lo, frag_hi) from a wgsim/dwgsim read name, or None."""
    name = qname.rsplit("/", 1)[0] if "/" in qname else qname
    # dwgsim random reads: name starts with "rand_"
    if name.startswith("rand_"):
        return None
    m = WGSIM_RE.match(name)
    if m:
        pos1, pos2 = int(m.group(2)), int(m.group(3))
        return m.group(1), min(pos1, pos2), max(pos1, pos2)
    return None


def is_correct(aln_pos, frag_lo, frag_hi, read_len, tolerance):
    """A read is correct if it starts near any candidate true position.

    Covers both simulator conventions (see module docstring):
      * frag_lo, frag_hi              — dwgsim: the two mates' own leftmost positions
      * frag_hi - read_len + 1        — wgsim: the right-end mate of a fragment span
    All 1-based leftmost coordinates; candidates are fragment-length apart so they
    never collide within the tolerance.
    """
    candidates = (frag_lo, frag_hi, frag_hi - read_len + 1)
    return any(abs(aln_pos - c) <= tolerance for c in candidates)


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
        if len(fields) < 10:
            continue
        qname = fields[0]
        flag = int(fields[1])
        # Skip secondary/supplementary
        if flag & 0x900:
            continue
        aln_pos = int(fields[3])  # 1-based leftmost mapping position
        # SEQ length = read length; '*' (absent) falls back to 0 -> right_start check no-ops.
        read_len = len(fields[9]) if fields[9] != "*" else 0
        n_mapped += 1

        parsed = parse_fragment(qname)
        if parsed is None:
            continue  # rand_/unparseable: skip (can't evaluate)
        n_parseable += 1
        _, frag_lo, frag_hi = parsed
        if is_correct(aln_pos, frag_lo, frag_hi, read_len, args.tolerance):
            n_correct += 1

    proc.wait()
    pa = n_correct / n_parseable if n_parseable > 0 else 0.0
    print(f"{pa:.4f}\t{n_mapped}\t{n_correct}\t{n_mapped}")


if __name__ == "__main__":
    main()
