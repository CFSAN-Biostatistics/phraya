#!/usr/bin/env python3
"""Generate a synthetic reference + simulated reads for local perf benchmarking.

Dependency-free and fully seeded, so a given set of parameters always produces
byte-identical inputs. This lets us compare `phraya align` speed/RSS across code
revisions on the exact same data, and diff the resulting `.phraya` output to prove
an optimization did not change results.

Reads are single-end: fragments of the reference at random positions, mutated with
independent per-base substitutions at the requested divergence. Quality is a flat
Q40 ('I'). That is enough to exercise the Case 2 (reads vs reference) hot path.
"""
import argparse
import random

BASES = "ACGT"
# Complement for the (optional) reverse-strand reads.
COMP = str.maketrans("ACGT", "TGCA")


def gen_reference(size: int, rng: random.Random) -> str:
    # random.choices over the whole genome in one call is far faster than a loop.
    return "".join(rng.choices(BASES, k=size))


def write_fasta(path: str, name: str, seq: str, width: int = 70) -> None:
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i : i + width])
            fh.write("\n")


def mutate(frag: str, divergence: float, rng: random.Random) -> str:
    if divergence <= 0.0:
        return frag
    out = list(frag)
    for i, b in enumerate(out):
        if rng.random() < divergence:
            # Substitute to a different base.
            alt = rng.choice(BASES)
            while alt == b:
                alt = rng.choice(BASES)
            out[i] = alt
    return "".join(out)


def write_reads(
    path: str,
    reference: str,
    num_reads: int,
    read_len: int,
    divergence: float,
    rng: random.Random,
) -> None:
    genome = len(reference)
    max_start = genome - read_len
    if max_start < 0:
        raise SystemExit(f"read_len {read_len} exceeds genome size {genome}")
    qual = "I" * read_len
    with open(path, "w") as fh:
        for i in range(num_reads):
            start = rng.randint(0, max_start)
            frag = reference[start : start + read_len]
            if rng.random() < 0.5:
                frag = frag.translate(COMP)[::-1]  # reverse-complement strand
            frag = mutate(frag, divergence, rng)
            fh.write(f"@read_{i} pos={start}\n{frag}\n+\n{qual}\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome-size", type=int, default=2_000_000)
    ap.add_argument("--num-reads", type=int, default=20_000)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--divergence", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ref-out", required=True)
    ap.add_argument("--reads-out", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    reference = gen_reference(args.genome_size, rng)
    write_fasta(args.ref_out, "synthetic_ref", reference)
    write_reads(
        args.reads_out,
        reference,
        args.num_reads,
        args.read_len,
        args.divergence,
        rng,
    )
    print(
        f"wrote {args.ref_out} ({args.genome_size} bp) and "
        f"{args.reads_out} ({args.num_reads} x {args.read_len} bp, "
        f"div={args.divergence}, seed={args.seed})"
    )


if __name__ == "__main__":
    main()
