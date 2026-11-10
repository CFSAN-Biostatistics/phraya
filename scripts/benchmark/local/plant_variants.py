#!/usr/bin/env python3
"""Plant SNPs, SNP clusters and small indels in a real assembly to make a mid-band pair.

The CSP2_Test assemblies are bimodal: within a clade pairs differ by 0-7 variants, between
clades by more than 20,000, so no natural pair sits in the 100-10,000 range CSP2 actually
screens (SNP_Cutoff halts a pair above 10,000 SNPs + indels + invalid sites). This makes one:
the query is the reference with N events planted at random, so a Phraya-vs-dnadiff comparison
exercises CSP2's identity, length and density filters on a real genome.

Background events: 83% single SNPs, 15% clusters of 4 SNPs inside 60 bp (so the density
windows have something to purge), 2% indels of 1-6 bases. Events keep 30 bp apart and 200 bp
from contig ends.

Optional events probe the two known Phraya-vs-dnadiff differences. They draw from their own
RNG and sit at least 2,000 bp from contig ends and 1,100 bp from every other event (clear of
CSP2's 1,000 bp density window), so with both counts at 0 the output is byte-identical to
earlier versions of this script.

  --breaks N   a 400 bp random insertion, which splits an alignment block, with one SNP d bp
               from it (d cycles 1..50, side alternates): Phraya's blocks start at the first
               exact seed, so SNPs close to an internal block edge can be lost.
  --islands N  four indels within 35 bp (alternately all insertions of 4-17 bp, or
               insertion/deletion/insertion/deletion) with flank SNPs 90 and 250 bp after the
               island start: aligners that represent part of the island as mismatches can push
               the flanks over CSP2's density limit.
  --truth TSV  every planted event: ref_contig, ref_pos (1-based), ref_base, alt_base, type
               (SNP/INS/DEL), callable, tag, detail. The first six columns are the
               snpdiffs_concordance.py --truth format.

Deterministic for a given seed. Stdlib only.

Usage: plant_variants.py REF.fasta OUT.fasta N_EVENTS SEED [--breaks N] [--islands N] [--truth TSV]
"""
import argparse
import random
from bisect import bisect_left, insort

BREAK_INSERT = 400
EDGE_DISTANCES = (1, 2, 3, 4, 6, 8, 11, 15, 20, 26, 33, 41, 50)
ISLAND_OFFSETS = (0, 13, 20, 29)
ISLAND_FLANKS = (90, 250)
END_MARGIN = 2000
EVENT_MARGIN = 1100


def main(src: str, out: str, n: int, seed: int, breaks: int = 0, islands: int = 0, truth_path: str = None) -> None:
    rng = random.Random(seed)
    names, seqs = [], {}
    for line in open(src):
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            seqs[names[-1]] = []
        else:
            seqs[names[-1]].append(line.strip().upper())
    seqs = {name: list("".join(parts)) for name, parts in seqs.items()}
    weights = [max(0, len(seqs[name]) - 400) for name in names]

    used = set()
    events = []

    def pick():
        while True:
            contig = rng.choices(names, weights)[0]
            pos = rng.randrange(200, len(seqs[contig]) - 200)
            if all((contig, pos + d) not in used for d in range(-30, 31)):
                return contig, pos

    while len(events) < n:
        r = rng.random()
        contig, pos = pick()
        if r < 0.02:
            events.append(("indel", contig, pos, rng.randint(1, 6), rng.choice("ID")))
            used.add((contig, pos))
        elif r < 0.17:
            for d in sorted(rng.sample(range(60), 4)):
                events.append(("snp", contig, pos + d, "cluster"))
                used.add((contig, pos + d))
        else:
            events.append(("snp", contig, pos, "single"))
            used.add((contig, pos))

    truth = []  # (contig, 0-based pos, ref_base, alt_base, type, tag, detail)

    def plant_snp(r, contig, pos, tag, detail=""):
        base = seqs[contig][pos]
        if base not in "ACGT":
            return False
        alt = r.choice([b for b in "ACGT" if b != base])
        seqs[contig][pos] = alt
        truth.append((contig, pos, base, alt, "SNP", tag, detail))
        return True

    planted = sum(plant_snp(rng, e[1], e[2], e[3]) for e in events if e[0] == "snp")
    indels = [e for e in events if e[0] == "indel"]
    # Edits are applied right to left once every event is placed; background insertion bases are
    # drawn here in the same order as before, so the background RNG stream is unchanged.
    edits = []
    for _, contig, pos, length, kind in sorted(indels, key=lambda e: (e[1], -e[2])):
        if kind == "I":
            edits.append((contig, pos, "I", [rng.choice("ACGT") for _ in range(length)], "indel", ""))
        else:
            edits.append((contig, pos, "D", length, "indel", ""))

    extra = {"breaks": 0, "edge SNPs": 0, "islands": 0, "flank SNPs": 0}
    if breaks or islands:
        rng2 = random.Random(seed * 1_000_003 + 7)
        reserved = {name: [] for name in names}
        for contig, pos in used:
            insort(reserved[contig], pos)

        def place(lo_off, hi_off):
            span = hi_off - lo_off
            w = [max(0, len(seqs[name]) - 2 * END_MARGIN - span) for name in names]
            for _ in range(200_000):
                contig = rng2.choices(names, w)[0]
                anchor = rng2.randrange(END_MARGIN - lo_off, len(seqs[contig]) - END_MARGIN - hi_off)
                lo, hi = anchor + lo_off, anchor + hi_off
                arr = reserved[contig]
                i = bisect_left(arr, lo - EVENT_MARGIN)
                if i == len(arr) or arr[i] > hi + EVENT_MARGIN:
                    insort(arr, lo)
                    insort(arr, hi)
                    return contig, anchor
            raise SystemExit("plant_variants: no room left for more breaks/islands; lower --breaks/--islands")

        for i in range(breaks):
            d = EDGE_DISTANCES[(i // 2) % len(EDGE_DISTANCES)]
            side = "L" if i % 2 == 0 else "R"
            contig, b = place(-EDGE_DISTANCES[-1], EDGE_DISTANCES[-1])
            edits.append((contig, b, "I", [rng2.choice("ACGT") for _ in range(BREAK_INSERT)], "break", ""))
            extra["breaks"] += 1
            if plant_snp(rng2, contig, b - d if side == "L" else b + d - 1, "edge", f"{side}{d}"):
                extra["edge SNPs"] += 1

        for i in range(islands):
            kind = "ins" if i % 2 == 0 else "mix"
            contig, s = place(0, ISLAND_FLANKS[-1])
            for j, off in enumerate(ISLAND_OFFSETS):
                if kind == "mix" and j % 2 == 1:
                    edits.append((contig, s + off, "D", rng2.randint(3, 6), "island", kind))
                else:
                    edits.append((contig, s + off, "I", [rng2.choice("ACGT") for _ in range(rng2.randint(4, 17))], "island", kind))
            extra["islands"] += 1
            for off in ISLAND_FLANKS:
                if plant_snp(rng2, contig, s + off, "island_flank", f"{kind}:{s + 1}"):
                    extra["flank SNPs"] += 1

    for contig, pos, kind, what, tag, detail in sorted(edits, key=lambda e: (e[0], -e[1])):  # right to left
        if kind == "I":
            seqs[contig][pos:pos] = what
            truth.append((contig, pos - 1, ".", "".join(what), "INS", tag, detail))
        else:
            truth.append((contig, pos, "".join(seqs[contig][pos:pos + what]), ".", "DEL", tag, detail))
            del seqs[contig][pos:pos + what]

    with open(out, "w") as f:
        for name in names:
            f.write(f">{name}\n")
            s = "".join(seqs[name])
            for i in range(0, len(s), 80):
                f.write(s[i:i + 80] + "\n")
    print(f"{out}: planted SNPs {planted}, indels {len(indels)}")
    if breaks or islands:
        print("  " + ", ".join(f"{k} {v}" for k, v in extra.items()))
    if truth_path:
        with open(truth_path, "w") as f:
            f.write("ref_contig\tref_pos\tref_base\talt_base\ttype\tcallable\ttag\tdetail\n")
            for contig, pos, ref_base, alt, kind, tag, detail in sorted(truth, key=lambda t: (t[0], t[1])):
                f.write(f"{contig}\t{pos + 1}\t{ref_base}\t{alt}\t{kind}\t1\t{tag}\t{detail}\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref")
    ap.add_argument("out")
    ap.add_argument("n_events", type=int)
    ap.add_argument("seed", type=int)
    ap.add_argument("--breaks", type=int, default=0)
    ap.add_argument("--islands", type=int, default=0)
    ap.add_argument("--truth")
    a = ap.parse_args()
    main(a.ref, a.out, a.n_events, a.seed, a.breaks, a.islands, a.truth)
