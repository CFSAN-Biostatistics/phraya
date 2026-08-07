#!/usr/bin/env python3
"""Generate a synthetic reference proteome + query proteins for protein-space
benchmarking (ADR-0013).

Mirrors `gen_synthetic.py`'s design (dependency-free, fully seeded, byte-identical
across runs of the same parameters) but over the 20-letter canonical amino-acid
alphabet, with **no reverse-complement step** — protein has no reverse strand
(ADR-0013's decision: `align_read` skips dual-strand search entirely in protein
mode), so every query is generated and expected to align in a single orientation.

Reuses `mutate_with_truth`/`format_indel_events` from `gen_synthetic.py` rather
than reimplementing mutation, so DNA and protein data share one (tested)
mutation/truth-sidecar implementation and can't drift apart.

Shape: a reference "proteome" of `--num-proteins` independent random proteins
(mean length `--protein-len`), and `--num-queries` query proteins, each a
mutated (substitution + optional indel) full-length or fragment homolog of a
randomly chosen proteome protein — closer to how DIAMOND/BLASTP/Phraya-protein
would actually be exercised (whole-sequence homolog search) than a short-read
fragment-sampling model. `--query-len` below the source protein's length
produces a fragment (Case-2-style); set it >= `--protein-len` (the default) to
usually get whole-protein queries.
"""
import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gen_synthetic import format_indel_events, mutate_with_truth  # noqa: E402

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"  # 20 canonical residues, no ambiguity codes


def gen_protein(length: int, rng: random.Random) -> str:
    return "".join(rng.choices(AMINO_ACIDS, k=length))


def gen_proteome(num_proteins: int, mean_len: int, rng: random.Random) -> list[tuple[str, str]]:
    """Returns [(protein_id, sequence), ...]. Length varies +/-30% around mean_len
    (uniform), like real proteomes, so proteins aren't all suspiciously identical
    in size — avoids a length-based seeding artifact that a fixed-length proteome
    would introduce."""
    proteins = []
    lo, hi = max(10, int(mean_len * 0.7)), int(mean_len * 1.3)
    for i in range(num_proteins):
        length = rng.randint(lo, hi)
        proteins.append((f"protein_{i}", gen_protein(length, rng)))
    return proteins


def write_fasta_multi(path: str, records: list[tuple[str, str]], width: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width])
                fh.write("\n")


def write_queries(
    path: str,
    truth_path: str,
    proteome: list[tuple[str, str]],
    num_queries: int,
    query_len: int,
    divergence: float,
    indel_rate: float,
    indel_max_len: int,
    rng: random.Random,
) -> None:
    records = []
    with open(truth_path, "w") as truth_fh:
        truth_fh.write("query_id\tsource_protein\tstart\tn_subs\tindel_events\n")
        for i in range(num_queries):
            source_id, source_seq = rng.choice(proteome)
            frag_len = min(query_len, len(source_seq))
            max_start = len(source_seq) - frag_len
            start = rng.randint(0, max_start) if max_start > 0 else 0
            frag = source_seq[start : start + frag_len]
            mutated, n_subs, indel_events = mutate_with_truth(
                frag, divergence, indel_rate, indel_max_len, rng, alphabet=AMINO_ACIDS
            )
            query_id = f"query_{i}"
            records.append((f"{query_id} source={source_id} pos={start}", mutated))
            truth_fh.write(
                f"{query_id}\t{source_id}\t{start}\t{n_subs}\t{format_indel_events(indel_events)}\n"
            )
    write_fasta_multi(path, records)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--num-proteins", type=int, default=2_000)
    ap.add_argument("--protein-len", type=int, default=300, help="Mean proteome protein length (aa)")
    ap.add_argument("--num-queries", type=int, default=500)
    ap.add_argument(
        "--query-len",
        type=int,
        default=300,
        help="Query length (aa); clamped to the source protein's length. "
        "Set below --protein-len for fragment queries, >= for whole-protein queries.",
    )
    ap.add_argument("--divergence", type=float, default=0.02, help="Point-substitution rate")
    ap.add_argument(
        "--indel-rate",
        type=float,
        default=0.0,
        help="Per-site probability of an indel event instead of a substitution "
        "(default 0.0 — no indels)",
    )
    ap.add_argument(
        "--indel-max-len",
        type=int,
        default=3,
        help="Max indel length in residues; actual length is uniform 1..=max (default 3)",
    )
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--proteome-out", required=True)
    ap.add_argument("--queries-out", required=True)
    args = ap.parse_args()

    if args.indel_rate < 0.0 or args.indel_rate > 1.0:
        raise SystemExit(f"--indel-rate must be in [0, 1], got {args.indel_rate}")
    if args.indel_max_len < 1:
        raise SystemExit(f"--indel-max-len must be >= 1, got {args.indel_max_len}")

    rng = random.Random(args.seed)
    proteome = gen_proteome(args.num_proteins, args.protein_len, rng)
    write_fasta_multi(args.proteome_out, proteome)

    truth_out = f"{args.queries_out}.truth.tsv"
    write_queries(
        args.queries_out,
        truth_out,
        proteome,
        args.num_queries,
        args.query_len,
        args.divergence,
        args.indel_rate,
        args.indel_max_len,
        rng,
    )
    print(
        f"wrote {args.proteome_out} ({args.num_proteins} proteins, mean {args.protein_len} aa), "
        f"{args.queries_out} ({args.num_queries} x ~{args.query_len} aa, "
        f"div={args.divergence}, indel_rate={args.indel_rate}, "
        f"indel_max_len={args.indel_max_len}, seed={args.seed}), and {truth_out}"
    )


if __name__ == "__main__":
    main()
