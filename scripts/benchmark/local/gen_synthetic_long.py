#!/usr/bin/env python3
"""Generate synthetic long reads (ONT/PacBio style) for long-read benchmarking.

Produces a random reference sequence and a set of long reads sampled from it,
mutated with a realistic long-read error profile: indels dominate (5–15% total
error rate, ~70% of errors are indels with a geometric length distribution),
substitutions are a smaller fraction (~1–3%), and reads range 10–50 kb.

Unlike gen_synthetic.py (short-read, wgsim-style), the truth sidecar records the
*original fragment* coordinates in reference space, not a fragment-span encoding
— long reads are single molecules, so each read has a single true position.

Truth sidecar format (TSV, one row per read):
    read_id  chrom  true_start  strand  true_end  n_subs  n_ins  n_del  indel_events

`indel_events` is a `;-` separated list of `offset:op:len` (1-based offset into
the original fragment, before mutation). `op` ∈ {I, D}. `len` is the event length.

Reads are single-end (one molecule per read).  ~50% are reverse-complemented
to simulate strand sampling.  Quality scores are a flat Q12 ('?') matching a
typical ONT R10.4 / PacBio HiFi-polished raw read.

Usage:
    gen_synthetic_long.py --genome-size 5000000 --num-reads 200 \
        --min-len 10000 --max-len 30000 \
        --read-error 0.12 --indel-frac 0.70 \
        --ref-out ref.fa --reads-out reads.fq
"""
import argparse
import random
import sys

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def gen_reference(size: int, rng: random.Random) -> str:
    """Generate a random DNA sequence of the given size."""
    # Chunk the choices call to avoid huge memory spikes for very large genomes.
    chunk = 1_000_000
    parts = []
    remaining = size
    while remaining > 0:
        n = min(chunk, remaining)
        parts.append("".join(rng.choices(BASES, k=n)))
        remaining -= n
    return "".join(parts)


def write_fasta(path: str, name: str, seq: str, width: int = 70) -> None:
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i : i + width])
            fh.write("\n")


def mutate_with_truth(
    frag: str,
    read_error: float,
    indel_frac: float,
    max_indel_len: int,
    rng: random.Random,
) -> tuple[str, int, int, int, str]:
    """Mutate a fragment with a long-read error model.

    Walks ``frag`` left to right. At each position, with probability
    ``read_error`` an error fires; ``indel_frac`` of errors are indels
    (insertion or deletion, 50/50), the rest are substitutions.

    Indel length is drawn from a geometric distribution (p=0.5, capped at
    ``max_indel_len``) — a simple model that produces the long-tailed length
    spectrum seen in ONT/PacBio data (most indels are 1–5 bp, a few are tens).

    Returns ``(mutated_seq, n_subs, n_ins, n_del, indel_events_str)``.
    """
    out: list[str] = []
    n_subs = 0
    n_ins = 0
    n_del = 0
    indel_events: list[str] = []
    n = len(frag)
    i = 0
    while i < n:
        if rng.random() < read_error:
            if rng.random() < indel_frac:
                # Indel
                is_insertion = rng.random() < 0.5
                # Geometric length: p=0.5, capped
                length = 1
                while length < max_indel_len and rng.random() < 0.5:
                    length += 1
                if is_insertion:
                    out.append("".join(rng.choices(BASES, k=length)))
                    n_ins += length
                    indel_events.append(f"{i + 1}:I:{length}")
                    # do NOT advance i — insertion is before the current base
                    continue
                else:
                    actual_len = min(length, n - i)
                    indel_events.append(f"{i + 1}:D:{actual_len}")
                    n_del += actual_len
                    i += actual_len
                    continue
            else:
                # Substitution
                b = frag[i]
                alt = rng.choice(BASES)
                while alt == b:
                    alt = rng.choice(BASES)
                out.append(alt)
                n_subs += 1
        else:
            out.append(frag[i])
        i += 1

    return "".join(out), n_subs, n_ins, n_del, ";".join(indel_events)


def write_reads(
    path: str,
    truth_path: str,
    reference: str,
    num_reads: int,
    min_len: int,
    max_len: int,
    read_error: float,
    indel_frac: float,
    max_indel_len: int,
    rng: random.Random,
    chrom_name: str,
) -> None:
    genome = len(reference)
    max_start = genome - 1  # just need to be within the genome

    with open(path, "w") as fh, open(truth_path, "w") as truth_fh:
        truth_fh.write(
            "read_id\tchrom\ttrue_start\tstrand\ttrue_end\tn_subs\tn_ins\tn_del\tindel_events\n"
        )
        for i in range(num_reads):
            read_len = rng.randint(min_len, max_len)
            # Ensure the fragment fits within the genome
            max_pos = genome - read_len
            if max_pos < 0:
                # Read longer than genome — skip (shouldn't happen with sane params)
                continue
            start = rng.randint(0, max_pos)
            frag = reference[start : start + read_len]
            strand = "+"
            if rng.random() < 0.5:
                frag = frag.translate(COMP)[::-1]
                strand = "-"

            mutated, n_subs, n_ins, n_del, events = mutate_with_truth(
                frag, read_error, indel_frac, max_indel_len, rng
            )

            read_id = f"longread_{i}"
            qual_char = "?"  # Q12 (~6% base error floor)
            fh.write(f"@{read_id} pos={start} strand={strand}\n{mutated}\n+\n{qual_char * len(mutated)}\n")
            truth_fh.write(
                f"{read_id}\t{chrom_name}\t{start}\t{strand}\t{start + read_len}\t"
                f"{n_subs}\t{n_ins}\t{n_del}\t{events}\n"
            )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome-size", type=int, default=5_000_000)
    ap.add_argument("--num-reads", type=int, default=200)
    ap.add_argument("--min-len", type=int, default=10_000)
    ap.add_argument("--max-len", type=int, default=30_000)
    ap.add_argument(
        "--read-error",
        type=float,
        default=0.12,
        help="Overall per-base error rate (default 0.12 = 12%%, typical ONT R10.4 raw)",
    )
    ap.add_argument(
        "--indel-frac",
        type=float,
        default=0.70,
        help="Fraction of errors that are indels (default 0.70)",
    )
    ap.add_argument(
        "--max-indel-len",
        type=int,
        default=50,
        help="Maximum indel length in bp (default 50; geometric length distribution)",
    )
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--chrom-name", default="synthetic_chr")
    ap.add_argument("--ref-out", required=True)
    ap.add_argument("--reads-out", required=True)
    args = ap.parse_args()

    if args.read_error < 0.0 or args.read_error > 1.0:
        raise SystemExit(f"--read-error must be in [0, 1], got {args.read_error}")
    if args.indel_frac < 0.0 or args.indel_frac > 1.0:
        raise SystemExit(f"--indel-frac must be in [0, 1], got {args.indel_frac}")
    if args.min_len > args.max_len:
        raise SystemExit(f"--min-len ({args.min_len}) > --max-len ({args.max_len})")
    if args.genome_size < args.max_len:
        raise SystemExit(f"genome-size ({args.genome_size}) < max-len ({args.max_len})")

    rng = random.Random(args.seed)
    reference = gen_reference(args.genome_size, rng)
    write_fasta(args.ref_out, args.chrom_name, reference)
    truth_out = f"{args.reads_out}.truth.tsv"
    write_reads(
        args.reads_out,
        truth_out,
        reference,
        args.num_reads,
        args.min_len,
        args.max_len,
        args.read_error,
        args.indel_frac,
        args.max_indel_len,
        rng,
        args.chrom_name,
    )
    print(
        f"wrote {args.ref_out} ({args.genome_size} bp), "
        f"{args.reads_out} ({args.num_reads} reads, {args.min_len}-{args.max_len} bp, "
        f"err={args.read_error}, indel_frac={args.indel_frac}, seed={args.seed}), "
        f"and {truth_out}"
    )


if __name__ == "__main__":
    main()
