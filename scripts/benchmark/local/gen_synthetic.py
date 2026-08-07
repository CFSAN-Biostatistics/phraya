#!/usr/bin/env python3
"""Generate a synthetic reference + simulated reads for local perf benchmarking.

Dependency-free and fully seeded, so a given set of parameters always produces
byte-identical inputs. This lets us compare `phraya align` speed/RSS across code
revisions on the exact same data, and diff the resulting `.phraya` output to prove
an optimization did not change results.

Reads are single-end: fragments of the reference at random positions, mutated with
independent per-base substitutions at the requested divergence. Quality is a flat
Q40 ('I'). That is enough to exercise the Case 2 (reads vs reference) hot path.

Optionally, sites can also mutate as indel events (`--indel-rate`) instead of
substitutions — mutually exclusive per site, so `--divergence` keeps meaning "rate
of point substitutions" regardless of whether indels are enabled. A `<reads-out>.truth.tsv`
sidecar records the true mutation events per read (substitution count + indel
events), consumed by `compute_indel_recovery.py` to score whether an aligner's
reported CIGAR recovers indels as single ops (see BENCHMARK_EXPANSION.md, "Indel
Event Concordance"). `--indel-rate 0` (the default) never calls the extra RNG draw,
so existing substitution-only invocations reproduce byte-identical output to before
indel support existed.
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


def mutate_with_truth(
    frag: str,
    divergence: float,
    indel_rate: float,
    indel_max_len: int,
    rng: random.Random,
    alphabet: str = BASES,
) -> tuple[str, int, list[tuple[int, str, int]]]:
    """Mutate `frag`, returning (mutated_seq, n_subs, indel_events).

    Walks `frag` left to right. At each original-fragment offset, at most one
    event fires: an indel (checked first, via `indel_rate`) or a substitution
    (via `divergence`) — never both at the same site, so `divergence` keeps its
    existing meaning (point-substitution rate) independent of whether indels are
    enabled at all. `indel_events` records `(offset, op, length)` with `offset`
    in *original fragment* coordinates (pre-mutation) and `op` in `{"I", "D"}`.

    `alphabet` is the symbol set substitutions/insertions draw from — `BASES`
    (ACGT) for DNA, the 20 canonical amino acids for `gen_synthetic_protein.py`
    (which imports this function rather than reimplementing it).

    When `indel_rate <= 0.0`, the indel check's `rng.random()` call is skipped
    entirely (short-circuited), so the substitution RNG stream — and therefore
    the output — is byte-identical to the pre-indel-support `mutate()`.
    """
    out: list[str] = []
    n_subs = 0
    indel_events: list[tuple[int, str, int]] = []
    n = len(frag)
    i = 0
    while i < n:
        b = frag[i]
        if indel_rate > 0.0 and rng.random() < indel_rate:
            length = rng.randint(1, indel_max_len)
            if rng.random() < 0.5:
                # Insertion: splice `length` random symbols before this site; the
                # original symbol at `i` is still processed next iteration.
                out.append("".join(rng.choices(alphabet, k=length)))
                indel_events.append((i, "I", length))
                continue
            # Deletion: drop up to `length` original symbols (bounded by what's left).
            length = min(length, n - i)
            indel_events.append((i, "D", length))
            i += length
            continue
        if divergence > 0.0 and rng.random() < divergence:
            alt = rng.choice(alphabet)
            while alt == b:
                alt = rng.choice(alphabet)
            out.append(alt)
            n_subs += 1
        else:
            out.append(b)
        i += 1
    return "".join(out), n_subs, indel_events


def format_indel_events(events: list[tuple[int, str, int]]) -> str:
    return ";".join(f"{offset}:{op}:{length}" for offset, op, length in events)


def write_reads(
    path: str,
    truth_path: str,
    reference: str,
    num_reads: int,
    read_len: int,
    divergence: float,
    indel_rate: float,
    indel_max_len: int,
    rng: random.Random,
) -> None:
    genome = len(reference)
    max_start = genome - read_len
    if max_start < 0:
        raise SystemExit(f"read_len {read_len} exceeds genome size {genome}")
    qual_char = "I"
    with open(path, "w") as fh, open(truth_path, "w") as truth_fh:
        truth_fh.write("read_id\tstart\tstrand\tn_subs\tindel_events\n")
        for i in range(num_reads):
            start = rng.randint(0, max_start)
            frag = reference[start : start + read_len]
            strand = "+"
            if rng.random() < 0.5:
                frag = frag.translate(COMP)[::-1]  # reverse-complement strand
                strand = "-"
            frag, n_subs, indel_events = mutate_with_truth(
                frag, divergence, indel_rate, indel_max_len, rng
            )
            read_id = f"read_{i}"
            fh.write(f"@{read_id} pos={start}\n{frag}\n+\n{qual_char * len(frag)}\n")
            truth_fh.write(
                f"{read_id}\t{start}\t{strand}\t{n_subs}\t{format_indel_events(indel_events)}\n"
            )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome-size", type=int, default=2_000_000)
    ap.add_argument("--num-reads", type=int, default=20_000)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--divergence", type=float, default=0.01)
    ap.add_argument(
        "--indel-rate",
        type=float,
        default=0.0,
        help="Per-site probability of an indel event instead of a substitution "
        "(default 0.0 — no indels, byte-identical to the substitution-only generator)",
    )
    ap.add_argument(
        "--indel-max-len",
        type=int,
        default=5,
        help="Max indel length in bp; actual length is uniform 1..=max (default 5)",
    )
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ref-out", required=True)
    ap.add_argument("--reads-out", required=True)
    args = ap.parse_args()

    if args.indel_rate < 0.0 or args.indel_rate > 1.0:
        raise SystemExit(f"--indel-rate must be in [0, 1], got {args.indel_rate}")
    if args.indel_max_len < 1:
        raise SystemExit(f"--indel-max-len must be >= 1, got {args.indel_max_len}")

    rng = random.Random(args.seed)
    reference = gen_reference(args.genome_size, rng)
    write_fasta(args.ref_out, "synthetic_ref", reference)
    truth_out = f"{args.reads_out}.truth.tsv"
    write_reads(
        args.reads_out,
        truth_out,
        reference,
        args.num_reads,
        args.read_len,
        args.divergence,
        args.indel_rate,
        args.indel_max_len,
        rng,
    )
    print(
        f"wrote {args.ref_out} ({args.genome_size} bp), "
        f"{args.reads_out} ({args.num_reads} x {args.read_len} bp, "
        f"div={args.divergence}, indel_rate={args.indel_rate}, "
        f"indel_max_len={args.indel_max_len}, seed={args.seed}), and {truth_out}"
    )


if __name__ == "__main__":
    main()
