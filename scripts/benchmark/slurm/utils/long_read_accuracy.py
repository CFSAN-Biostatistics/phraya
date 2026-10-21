#!/usr/bin/env python3
"""Placement accuracy (PA) for long-read alignment, scored against ground truth.

Truth TSV (from gen_synthetic_long.py / data-commons gen_long_reads.py) columns:
    read_id  chrom  true_start  strand  true_end  n_subs  n_ins  n_del  indel_events
`true_start` is the leftmost reference coordinate (0-based) of the simulated read; `true_end`
is half-open.

Two input kinds, chosen by file extension (no content sniffing, no silent fallback):
    *.queries        phraya cross-space queries sidecar: zstd(msgpack({read_id: [[space, pos, identity], ...]}))
    *.sam            SAM text (e.g. minimap2 -a). Primary records give the best placement;
                     primary + secondary records give the retained placements (flag 0x4
                     unmapped and 0x800 supplementary are never used).

The same rule is applied to both tools. Per truth read:
    mapped        >= 1 placement/record (on any reference)
    best          the read's single best placement.
                  phraya: highest-identity placement (earliest position on ties).
                  SAM:    the primary record; the leading soft/hard clip is subtracted so the
                          position is where the read's first base would sit.
    correct       best is on the true chrom and its START is within --tolerance bp of true_start.
    correct_any   at least one retained placement (phraya: every placement in the sidecar;
                  SAM: primary + secondary records) is on the true chrom within tolerance.

    PA        = correct / n_truth       (unmapped reads count against the aligner)
    precision = correct / mapped

Choosing --tolerance: a start test measures two things at once, locus correctness and how
precisely the aligner reproduces the read's first base. Noisy read ends are soft-clipped, so a
tight tolerance (100 bp) penalises clipping, and more so at higher error rates (minimap2 on the
data-commons T9 sets: within 100 bp for 0.48-0.96 of reads, within 2 kb for ~1.0). The harness
therefore reports PA at both: ~2000 bp ("locus correct", used for CAS) and 100 bp ("start
precision"). Run this script once per tolerance.

Optional modes:
    --segdup-bed BED   adds four fields for a planted-repeat target (data-commons
                       segdup_truth.bed: chrom start end name identity). Reads wholly inside a copy
                       are an intentional near-tie; reads overhanging into unique flank are not.
                         n_repeat_inside      truth reads wholly inside a repeat copy
                         n_repeat_both_copies ... of those, with a retained placement at the true
                                              position AND at the mirrored position in every
                                              other copy (within tolerance)
                         n_overhang           truth reads overlapping a copy but not wholly inside
                         n_overhang_correct   ... of those, with `correct` (unambiguous subset)
    --agree-with SAM   no-truth mode (real reads): compare this aligner's best placement with the
                       primary placements in SAM (use minimap2). Prints
                       n_both n_agree n_only_this n_only_ref. The truth argument is not used.

Output: one TSV line, or `--header` to also print the field names.

Usage:
    long_read_accuracy.py <cross_space.phraya.queries|aln.sam> <reads.fq.truth.tsv> [--tolerance 100] [--segdup-bed BED] [--header]
    long_read_accuracy.py <cross_space.phraya.queries|aln.sam> --agree-with ref.sam [--tolerance 2000] [--header]
"""
import argparse
import csv
import re
import sys

FIELDS = ["pa", "precision", "n_truth", "n_mapped", "n_correct", "n_correct_any"]
SEGDUP_FIELDS = ["n_repeat_inside", "n_repeat_both_copies", "n_overhang", "n_overhang_correct"]
AGREE_FIELDS = ["n_both", "n_agree", "n_only_this", "n_only_ref"]
SAM_UNUSABLE_FLAGS = 0x4 | 0x800
CLIP_RE = re.compile(r"^(\d+)[SH]")


def chrom_to_label(chrom: str) -> str:
    """Mirror phraya-cli's sanitize_label_component: alnum/-/_/. pass through, else '_'."""
    return re.sub(r"[^0-9A-Za-z._-]", "_", chrom)


def load_truth(path: str) -> dict:
    """{read_id: (chrom_label, true_start, true_end)}"""
    with open(path, newline="") as f:
        return {
            row["read_id"]: (chrom_to_label(row["chrom"]), int(row["true_start"]), int(row["true_end"]))
            for row in csv.DictReader(f, delimiter="\t")
        }


def load_segdup_bed(path: str) -> list:
    """[(chrom_label, start, end)] for each repeat copy."""
    copies = []
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                p = line.split("\t")
                copies.append((chrom_to_label(p[0]), int(p[1]), int(p[2])))
    return copies


def load_phraya_queries(path: str) -> dict:
    """{read_id: (best (label, pos), [(label, pos), ...])} from the phraya sidecar."""
    import msgpack
    import zstandard

    with open(path, "rb") as f:
        # Frames written by phraya carry no content size; stream-decode.
        raw = zstandard.ZstdDecompressor().stream_reader(f).read()
    out = {}
    for rid, placements in msgpack.unpackb(raw, raw=False).items():
        if placements:
            space, pos, _ = max(placements, key=lambda p: (p[2], -p[1]))
            out[rid] = ((space, pos), [(p[0], p[1]) for p in placements])
    return out


def load_sam(path: str) -> dict:
    """{read_id: (best (label, clip-adjusted 0-based start), [(label, start), ...])}.
    best = the primary record; the list holds primary + secondary records."""
    best, allp = {}, {}
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"malformed SAM record: {line[:80]!r}")
            flag = int(fields[1])
            if flag & SAM_UNUSABLE_FLAGS or fields[2] == "*":
                continue
            clip = CLIP_RE.match(fields[5])
            place = (chrom_to_label(fields[2]), int(fields[3]) - 1 - (int(clip.group(1)) if clip else 0))
            allp.setdefault(fields[0], []).append(place)
            if not flag & 0x100:
                best[fields[0]] = place
    return {rid: (best[rid], allp[rid]) for rid in best}


def near(place, label, pos, tol):
    return place[0] == label and abs(place[1] - pos) <= tol


def score_truth(placements: dict, truth: dict, tol: int):
    mapped = correct = correct_any = 0
    for rid, (chrom, start, _end) in truth.items():
        hit = placements.get(rid)
        if hit is None:
            continue
        mapped += 1
        best, allp = hit
        if near(best, chrom, start, tol):
            correct += 1
        if any(near(p, chrom, start, tol) for p in allp):
            correct_any += 1
    return mapped, correct, correct_any


def score_segdup(placements: dict, truth: dict, copies: list, tol: int):
    inside = both = overhang = overhang_correct = 0
    for rid, (chrom, start, end) in truth.items():
        here = next((c for c in copies if c[0] == chrom and c[1] <= start and end <= c[2]), None)
        overlaps = here is not None or any(c[0] == chrom and start < c[2] and end > c[1] for c in copies)
        if not overlaps:
            continue
        hit = placements.get(rid)
        if here is not None:
            inside += 1
            if hit is not None:
                mirrors = [start + (c[1] - here[1]) for c in copies if c is not here]
                targets = [start] + mirrors
                if all(any(near(p, chrom, t, tol) for p in hit[1]) for t in targets):
                    both += 1
        else:
            overhang += 1
            if hit is not None and near(hit[0], chrom, start, tol):
                overhang_correct += 1
    return inside, both, overhang, overhang_correct


def agreement(this: dict, ref: dict, tol: int):
    both = [rid for rid in this if rid in ref]
    agree = sum(1 for rid in both if near(this[rid][0], ref[rid][0][0], ref[rid][0][1], tol))
    return len(both), agree, len(this) - len(both), len(ref) - len(both)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("alignments", help="*.queries (phraya sidecar) or *.sam")
    ap.add_argument("truth", nargs="?", help="truth TSV (omit with --agree-with)")
    ap.add_argument("--tolerance", type=int, default=100, help="bp (default 100)")
    ap.add_argument("--segdup-bed", help="planted-repeat BED; appends repeat fields")
    ap.add_argument("--agree-with", metavar="SAM", help="no-truth mode: agreement with this SAM's primary placements")
    ap.add_argument("--header", action="store_true", help="print field names first")
    args = ap.parse_args()

    if args.alignments.endswith(".queries"):
        placements = load_phraya_queries(args.alignments)
    elif args.alignments.endswith(".sam"):
        placements = load_sam(args.alignments)
    else:
        sys.exit(f"unsupported extension (want .queries or .sam): {args.alignments}")

    if args.agree_with:
        if not args.agree_with.endswith(".sam"):
            sys.exit(f"--agree-with wants a .sam: {args.agree_with}")
        if args.segdup_bed:
            sys.exit("--segdup-bed needs truth; it cannot be combined with --agree-with")
        n_both, n_agree, only_this, only_ref = agreement(placements, load_sam(args.agree_with), args.tolerance)
        if args.header:
            print("\t".join(AGREE_FIELDS))
        print(f"{n_both}\t{n_agree}\t{only_this}\t{only_ref}")
        return

    if not args.truth:
        sys.exit("truth TSV is required unless --agree-with is given")
    truth = load_truth(args.truth)
    mapped, correct, correct_any = score_truth(placements, truth, args.tolerance)
    n = len(truth)
    pa = correct / n if n else 0.0
    precision = correct / mapped if mapped else 0.0
    names, row = list(FIELDS), [f"{pa:.4f}", f"{precision:.4f}", n, mapped, correct, correct_any]
    if args.segdup_bed:
        names += SEGDUP_FIELDS
        row += score_segdup(placements, truth, load_segdup_bed(args.segdup_bed), args.tolerance)
    if args.header:
        print("\t".join(names))
    print("\t".join(str(v) for v in row))


if __name__ == "__main__":
    main()
