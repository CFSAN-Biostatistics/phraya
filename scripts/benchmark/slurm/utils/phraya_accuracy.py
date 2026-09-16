#!/usr/bin/env python3
"""
Compute placement accuracy (PA) for phraya from .phraya.queries + wgsim read names.

wgsim/dwgsim encode the true FRAGMENT span in the read name:
    <chrom>_<start>_<end>_<frag_info>/1  or /2

Two simulator conventions coexist in this benchmark's data (see sam_accuracy.py):
wgsim encodes the fragment span (mate /2 aligns at max - read_len + 1), dwgsim encodes the
two mates' own leftmost positions directly. Scoring every read against min() alone counts
every right-end mate as misplaced, halving PA into a ~0.49 artifact unrelated to the aligner.
So a read is correct if its best alignment start is within <tolerance> bp of ANY candidate
true start: min, max, or (max - read_len + 1). This matches sam_accuracy.py exactly, so
phraya and the SAM aligners are scored on the same basis.

Two sidecar schemas are supported, detected automatically from the file's content:
  * legacy, single-space: QueryIndex = HashMap<String, Vec<(u32, f64)>>
    msgpack decodes each placement as [pos, score].
  * cross-space (ADR-0011, issue #198): CrossSpaceQueryIndex =
    HashMap<String, Vec<CrossSpacePlacement>>, msgpack decodes each placement as
    [space, pos, identity] (phraya-io/src/queries.rs). `space` is the reference space's
    label — its name (sanitized the same way phraya-cli sanitizes filenames), or, for a
    sidecar written before reference spaces carried names, a bare content-hash prefix.
    Scoring requires both the position AND the space to match the read's true chromosome;
    with hundreds of reference spaces (e.g. a per-chromosome palette), position-only
    scoring would credit a read placed on the wrong chromosome at a coincidentally
    matching offset. A hash-labelled sidecar can never match by chromosome name, so this
    script falls back to position-only scoring for such files (with a warning) rather than
    reporting a spurious PA of 0.

PA = (correctly placed reads) / (placed, parseable, non-random reads).

Also reports unaligned fraction if total_reads is supplied.

Usage:
    phraya_accuracy.py <file.phraya.queries> [--tolerance 10] [--read-len 150] [--total-reads N]

Outputs tab-separated to stdout:
    pa  n_mapped  n_correct  n_unaligned  unaligned_frac
"""
import argparse
import re
import sys
import zstandard
import msgpack


# Matches wgsim and dwgsim read names. dwgsim uses start>end for rev-strand reads.
WGSIM_RE = re.compile(r"^(.+?)_(\d+)_(\d+)_")


def parse_fragment(read_name: str):
    """Return (chrom, frag_lo, frag_hi) from a wgsim/dwgsim read name, or None."""
    # Strip /1 /2 suffix if present
    name = read_name.rsplit("/", 1)[0] if "/" in read_name else read_name
    if name.startswith("rand_"):  # dwgsim random/unplaceable reads
        return None
    m = WGSIM_RE.match(name)
    if m:
        pos1, pos2 = int(m.group(2)), int(m.group(3))
        return m.group(1), min(pos1, pos2), max(pos1, pos2)
    return None


def label_for_chrom(chrom: str) -> str:
    """Mirror phraya-cli's sanitize_label_component: alnum/-/_/. pass through, else '_'."""
    return re.sub(r"[^0-9A-Za-z._-]", "_", chrom)


def strip_collision_suffix(space: str) -> str:
    """Undo phraya-cli's `<label>__<16-hex>` disambiguation suffix (added only when two
    reference space names sanitize to the same label) so a chromosome comparison isn't
    defeated by a collision-safe label."""
    return re.sub(r"__[0-9a-f]{16}$", "", space)


def score_pass(data: dict, read_len: int, tolerance: int, cross_space: bool, chrom_aware: bool):
    """One scoring pass over every read. `cross_space` selects the 3-element
    (space, pos, identity) decoding; `chrom_aware` additionally requires the placement's
    space to match the read's true chromosome (ignored when `cross_space` is False, since
    a legacy single-space sidecar carries no space label to compare)."""
    n_parseable = 0
    n_correct = 0
    for read_name, alignments in data.items():
        parsed = parse_fragment(read_name)
        if parsed is None:
            # rand_ reads or unparseable names — skip entirely (can't evaluate)
            continue
        if not alignments:
            continue  # no positions recorded — count as unmapped, skip
        n_parseable += 1
        chrom, frag_lo, frag_hi = parsed
        if cross_space:
            # Best alignment = highest identity (third element)
            space, pos, _ = max(alignments, key=lambda p: p[2])
            if chrom_aware and strip_collision_suffix(space) != label_for_chrom(chrom):
                continue
            pos = int(pos)
        else:
            # Best alignment = highest score (second element)
            pos = int(max(alignments, key=lambda p: p[1])[0])
        # Correct if near any candidate true start (covers wgsim span + dwgsim direct):
        # frag_lo, frag_hi, or the wgsim right-end read start (frag_hi - read_len + 1).
        candidates = (frag_lo, frag_hi, frag_hi - read_len + 1)
        if any(abs(pos - c) <= tolerance for c in candidates):
            n_correct += 1
    return n_parseable, n_correct


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("queries_file")
    parser.add_argument("--tolerance", type=int, default=10)
    parser.add_argument("--read-len", type=int, default=150)
    parser.add_argument("--total-reads", type=int, default=0)
    args = parser.parse_args()

    dctx = zstandard.ZstdDecompressor()
    with open(args.queries_file, "rb") as f:
        with dctx.stream_reader(f) as reader:
            raw = reader.read()
    data: dict = msgpack.unpackb(raw, raw=False)

    n_mapped = len(data)

    # Schema is uniform across one file: either every non-empty entry decodes to 3-element
    # (space, pos, identity) triples (cross-space sidecar) or 2-element (pos, score) pairs
    # (legacy single-space sidecar). Detect once from the first non-empty entry.
    cross_space_schema = any(len(v[0]) == 3 for v in data.values() if v)

    if cross_space_schema:
        n_parseable, n_correct = score_pass(
            data, args.read_len, args.tolerance, cross_space=True, chrom_aware=True
        )
        if n_correct == 0 and n_parseable > 0:
            _, n_correct_position_only = score_pass(
                data, args.read_len, args.tolerance, cross_space=True, chrom_aware=False
            )
            if n_correct_position_only > 0:
                print(
                    "warning: 0 reads matched by chromosome label — sidecar likely "
                    "predates named reference spaces (bare content-hash labels); "
                    "falling back to position-only scoring",
                    file=sys.stderr,
                )
                n_correct = n_correct_position_only
    else:
        n_parseable, n_correct = score_pass(
            data, args.read_len, args.tolerance, cross_space=False, chrom_aware=False
        )

    pa = n_correct / n_parseable if n_parseable > 0 else 0.0

    total = args.total_reads if args.total_reads > 0 else n_mapped
    n_unaligned = max(0, total - n_mapped)
    unaligned_frac = n_unaligned / total if total > 0 else 0.0

    print(f"{pa:.4f}\t{n_mapped}\t{n_correct}\t{n_unaligned}\t{unaligned_frac:.4f}")


if __name__ == "__main__":
    main()
