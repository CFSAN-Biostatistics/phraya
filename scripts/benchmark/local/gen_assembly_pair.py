#!/usr/bin/env python3
"""Synthetic assembly pairs with exact SNP truth (CSP2 handoff section 6.3, cases S1-S7).

A random genome is cut into reference contigs; the query is the same genome with planted
SNPs, cut independently. Writes to --out:

  ref.fa, qry.fa  the two assemblies
  truth.tsv       ref_contig, ref_pos (1-based), ref_base, alt_base, type, callable
  case.json       parameters and contig layout

`callable` is 1 when the SNP is at least --edge bp from every end of its reference contig
and of its query contig (a local aligner cannot place a mismatch at an alignment terminus).
Every planted SNP is listed, callable or not, so callers can score precision against all
of them and recall against the callable subset.

Cases:
  S1  identical breakpoints
  S2  shifted breakpoints: query contigs span 2-3 reference contigs
  S3  S2 with every other query contig reverse-complemented
  S4  S2 plus unrelated decoy contigs in both assemblies
  S5  identical breakpoints at contig lengths 10 kb, 50 kb, 200 kb, 1 Mb
  S6  S2 layout with --snp-count SNPs (run once per divergence level)
  S7  S1 plus a reference contig duplicating a region of another (one-to-one choice);
      SNPs inside the duplicated region are also listed on the duplicate copy with
      callable=0, so a call on either copy is not a false positive

Usage: gen_assembly_pair.py S2 --out DIR [--seed N] [--genome-size N] [--snp-period N]
"""
import argparse
import json
import os
import random

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def alt_base(base: str) -> str:
    return BASES[(BASES.index(base) + 1) % 4]


def random_cuts(rng: random.Random, length: int, pieces: int, min_len: int) -> list:
    """`pieces` contigs covering [0, length), each at least `min_len` long."""
    while True:
        inner = sorted(rng.sample(range(min_len, length - min_len), pieces - 1))
        points = [0] + inner + [length]
        if all(b - a >= min_len for a, b in zip(points, points[1:])):
            return points


def snp_positions(genome_len: int, period: int, count: int, rng: random.Random) -> list:
    if count:
        return sorted(rng.sample(range(genome_len), count))
    return list(range(period // 2, genome_len, period))


def avoid_cut_points(points: list, snps: set) -> list:
    """Nudge interior cut points so no planted SNP sits exactly at a contig boundary."""
    out = [points[0]]
    for p in points[1:-1]:
        while p in snps or p - 1 in snps:
            p += 1
        out.append(p)
    return out + [points[-1]]


def contig_index(points: list, g: int) -> int:
    lo, hi = 0, len(points) - 2
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if points[mid] <= g:
            lo = mid
        else:
            hi = mid - 1
    return lo


def build_case(case: str, seed: int, genome_size: int, period: int, snp_count: int, edge: int) -> dict:
    rng = random.Random(seed)
    if case == "S5":
        lengths = [10_000, 50_000, 200_000, 1_000_000]
        genome_size = max(genome_size, sum(lengths) + 20_000)
    genome = "".join(rng.choice(BASES) for _ in range(genome_size))
    snps = snp_positions(genome_size, period, snp_count if case == "S6" else 0, rng)
    snp_set = set(snps)
    mutated = list(genome)
    for g in snps:
        mutated[g] = alt_base(genome[g])
    mutated = "".join(mutated)

    if case == "S5":
        ref_points = [0]
        for n in lengths:
            ref_points.append(ref_points[-1] + n)
        ref_points.append(genome_size)
        qry_points = list(ref_points)
    else:
        ref_points = random_cuts(rng, genome_size, 8, 5_000)
        if case in ("S1", "S7"):
            qry_points = list(ref_points)
        else:
            qry_points = random_cuts(rng, genome_size, 5, 20_000)
    ref_points = avoid_cut_points(ref_points, snp_set)
    qry_points = avoid_cut_points(qry_points, snp_set)

    reversed_qry = set(range(1, len(qry_points) - 1, 2)) if case == "S3" else set()
    ref_records = [(f"r{i}", genome[a:b]) for i, (a, b) in enumerate(zip(ref_points, ref_points[1:]))]
    qry_records = []
    for i, (a, b) in enumerate(zip(qry_points, qry_points[1:])):
        seq = mutated[a:b]
        qry_records.append((f"q{i}", revcomp(seq) if i in reversed_qry else seq))

    truth = []
    for g in snps:
        ri, qi = contig_index(ref_points, g), contig_index(qry_points, g)
        ra, rb = ref_points[ri], ref_points[ri + 1]
        qa, qb = qry_points[qi], qry_points[qi + 1]
        callable_ = min(g - ra, rb - 1 - g) >= edge and min(g - qa, qb - 1 - g) >= edge
        truth.append((f"r{ri}", g - ra + 1, genome[g], mutated[g], "SNP", int(callable_)))

    if case == "S4":
        for i in range(2):
            ref_records.append((f"decoy_ref{i}", "".join(rng.choice(BASES) for _ in range(20_000))))
            qry_records.append((f"decoy_qry{i}", "".join(rng.choice(BASES) for _ in range(20_000))))

    if case == "S7":
        # Duplicate 5 kb from the middle of the largest reference contig as its own contig.
        ri = max(range(len(ref_records)), key=lambda i: len(ref_records[i][1]))
        ra = ref_points[ri]
        mid = ra + len(ref_records[ri][1]) // 2
        dup_start, dup_end = mid - 2_500, mid + 2_500
        ref_records.append(("dup", genome[dup_start:dup_end]))
        for g in snps:
            if dup_start <= g < dup_end:
                truth.append(("dup", g - dup_start + 1, genome[g], mutated[g], "SNP", 0))

    return {
        "ref_records": ref_records,
        "qry_records": qry_records,
        "truth": truth,
        "meta": {
            "case": case,
            "seed": seed,
            "genome_size": genome_size,
            "snp_count": len(snps),
            "snp_period": None if case == "S6" else period,
            "edge": edge,
            "ref_contigs": [{"id": n, "length": len(s)} for n, s in ref_records],
            "qry_contigs": [
                {"id": n, "length": len(s), "reverse_complemented": i in reversed_qry}
                for i, (n, s) in enumerate(qry_records)
            ],
            "callable_snps": sum(t[5] for t in truth if t[0] != "dup"),
        },
    }


def write_fasta(path: str, records: list) -> None:
    with open(path, "w", newline="\n") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=["S1", "S2", "S3", "S4", "S5", "S6", "S7"])
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--genome-size", type=int, default=400_000)
    ap.add_argument("--snp-period", type=int, default=200, help="one SNP every N bp (all cases but S6)")
    ap.add_argument("--snp-count", type=int, default=1_000, help="planted SNPs for S6")
    ap.add_argument("--edge", type=int, default=50, help="callable distance from contig ends")
    args = ap.parse_args()

    built = build_case(args.case, args.seed, args.genome_size, args.snp_period, args.snp_count, args.edge)
    os.makedirs(args.out, exist_ok=True)
    write_fasta(os.path.join(args.out, "ref.fa"), built["ref_records"])
    write_fasta(os.path.join(args.out, "qry.fa"), built["qry_records"])
    with open(os.path.join(args.out, "truth.tsv"), "w", newline="\n") as f:
        f.write("ref_contig\tref_pos\tref_base\talt_base\ttype\tcallable\n")
        for row in built["truth"]:
            f.write("\t".join(str(x) for x in row) + "\n")
    with open(os.path.join(args.out, "case.json"), "w", newline="\n") as f:
        json.dump(built["meta"], f, indent=2)
    m = built["meta"]
    print(f"{m['case']}: {len(m['ref_contigs'])} ref contigs, {len(m['qry_contigs'])} query contigs, "
          f"{m['snp_count']} SNPs ({m['callable_snps']} callable) -> {args.out}")


if __name__ == "__main__":
    main()
