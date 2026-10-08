#!/usr/bin/env python3
"""Score a `.snpdiffs` file against exact synthetic truth or a dnadiff (1-to-1) baseline.

Implements the metrics of the CSP2 handoff (`phraya-csp2-snp-discrepancy.md` §6.2) and the
pass criteria of §6.5:

  1. SNP recall and precision by position and by position + allele
  2. SNP rows, truth-matching rows and indel rows per alignment-identity bucket
     (<60, 60-90, 90-99, >=99 Perc_Iden), which separates under-alignment from noise
  3. reference and query breadth: union of covered BED intervals over assembly bases,
     compared with dnadiff AlignedBases
  4. pair-level recall: (reference contig, query contig) pairs in dnadiff .1coords vs
     covered BED pairs at >= 90% identity
  5. reference positions reported in more than one SNP row
  6. header invariants: *_Percent_Aligned in [0, 100], header counts equal row counts,
     header breadth equal to the BED union

Truth sources (exactly one):
  --truth truth.tsv   from gen_assembly_pair.py; recall uses the callable subset when
                      --callable-only is set, precision always uses every planted SNP
  --dnadiff PREFIX    PREFIX.snps, PREFIX.1coords and PREFIX.report from dnadiff

Prints a JSON object on stdout and a short summary on stderr; exits 1 if --gate is set
and a pass criterion fails. Stdlib only.
"""
import argparse
import json
import sys
from collections import defaultdict

BUCKETS = ("<60", "60-90", "90-99", ">=99")
BED_COLUMNS = ("Ref_Contig", "Ref_Start", "Ref_End", "Ref_Length", "Ref_Aligned", "Query_Contig",
               "Query_Start", "Query_End", "Query_Length", "Query_Aligned", "Perc_Iden")
VARIANT_COLUMNS = ("Ref_Contig", "Start_Ref", "Ref_Pos", "Query_Contig", "Start_Query", "Query_Pos",
                   "Ref_Loc", "Query_Loc", "Ref_Start", "Ref_End", "Query_Start", "Query_End",
                   "Ref_Base", "Query_Base", "Dist_to_Ref_End", "Dist_to_Query_End", "Ref_Aligned",
                   "Query_Aligned", "Query_Direction", "Perc_Iden", "Cat")


def bucket(perc_iden: float) -> str:
    if perc_iden < 60:
        return "<60"
    if perc_iden < 90:
        return "60-90"
    if perc_iden < 99:
        return "90-99"
    return ">=99"


def parse_snpdiffs(path: str):
    header, bed, variants = {}, [], []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("##\t"):
                fields = line[3:].split("\t")
                if len(fields) != len(BED_COLUMNS):
                    raise ValueError(f"BED row with {len(fields)} columns: {line}")
                bed.append(dict(zip(BED_COLUMNS, fields)))
            elif line.startswith("#"):
                for item in line[1:].strip("\t").split("\t"):
                    if ":" in item:
                        key, value = item.split(":", 1)
                        header[key] = value
            else:
                fields = line.split("\t")
                if len(fields) != len(VARIANT_COLUMNS):
                    raise ValueError(f"variant row with {len(fields)} columns: {line}")
                variants.append(dict(zip(VARIANT_COLUMNS, fields)))
    return header, bed, variants


def covered(bed_rows):
    return [r for r in bed_rows if r["Ref_Start"] != "." and r["Query_Start"] != "."]


def union_length(intervals) -> int:
    """Total length of the union of (name, start, end) half-open intervals."""
    by_name = defaultdict(list)
    for name, start, end in intervals:
        by_name[name].append((start, end))
    total = 0
    for spans in by_name.values():
        spans.sort()
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start > cur_end:
                total += cur_end - cur_start
                cur_start, cur_end = start, end
            else:
                cur_end = max(cur_end, end)
        total += cur_end - cur_start
    return total


def load_truth_tsv(path: str):
    """Returns {(contig, pos): (ref_base, alt_base, callable)} for planted SNPs."""
    truth = {}
    with open(path) as f:
        columns = f.readline().rstrip("\r\n").split("\t")
        for line in f:
            row = dict(zip(columns, line.rstrip("\r\n").split("\t")))
            if row["type"] == "SNP":
                truth[(row["ref_contig"], int(row["ref_pos"]))] = (row["ref_base"], row["alt_base"], row["callable"] == "1")
    return truth


def load_dnadiff(prefix: str):
    """dnadiff 1-to-1 outputs: SNPs, indel count, coords blocks, report AlignedBases."""
    snps, indels = {}, 0
    with open(prefix + ".snps") as f:
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 12 or not fields[0].strip().isdigit():
                continue
            ref_pos, ref_base, qry_base, ref_contig = int(fields[0]), fields[1], fields[2], fields[10]
            if ref_base == "." or qry_base == ".":
                indels += 1
            else:
                snps[(ref_contig, ref_pos)] = (ref_base, qry_base, True)
    coords = []
    with open(prefix + ".1coords") as f:
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 13 and fields[0].strip().isdigit():
                coords.append({"ref": fields[11], "qry": fields[12], "iden": float(fields[6])})
    report = {}
    with open(prefix + ".report") as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3 and parts[0] not in report:
                report[parts[0]] = (parts[1].split("(")[0], parts[2].split("(")[0])
    return snps, indels, coords, report


def pct(num: float, den: float):
    return round(100.0 * num / den, 4) if den else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snpdiffs", required=True)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--truth", help="truth.tsv from gen_assembly_pair.py")
    src.add_argument("--dnadiff", help="dnadiff output prefix")
    ap.add_argument("--callable-only", action="store_true", help="synthetic recall over callable SNPs only")
    ap.add_argument("--gate", action="store_true", help="exit 1 when a pass criterion fails")
    ap.add_argument("--min-recall", type=float, default=0.99)
    ap.add_argument("--min-allele-concordance", type=float, default=0.99)
    ap.add_argument("--max-extra-fraction", type=float, default=0.01)
    ap.add_argument("--max-breadth-delta", type=float, default=1.0, help="percentage points vs dnadiff")
    args = ap.parse_args()

    header, bed, variants = parse_snpdiffs(args.snpdiffs)
    snp_rows = [v for v in variants if v["Cat"] == "SNP"]
    indel_rows = [v for v in variants if v["Cat"] == "Indel"]
    found = defaultdict(list)
    for r in snp_rows:
        found[(r["Ref_Contig"], int(r["Ref_Pos"]))].append(r)

    coords, report, truth_indels = [], {}, None
    if args.truth:
        truth_all = load_truth_tsv(args.truth)
    else:
        truth_all, truth_indels, coords, report = load_dnadiff(args.dnadiff)
    truth_recall = {k: v for k, v in truth_all.items() if v[2] or not args.callable_only}

    hits = [k for k in truth_recall if k in found]
    allele_hits = [k for k in hits if any(
        (r["Ref_Base"], r["Query_Base"]) == truth_recall[k][:2] for r in found[k])]
    extra = [k for k in found if k not in truth_all]
    out = {
        "snp_rows": len(snp_rows),
        "snp_positions": len(found),
        "indel_rows": len(indel_rows),
        "truth_snps": len(truth_recall),
        "truth_indels": truth_indels,
        "recall": round(len(hits) / len(truth_recall), 6) if truth_recall else None,
        "allele_concordance": round(len(allele_hits) / len(hits), 6) if hits else None,
        "recall_with_allele": round(len(allele_hits) / len(truth_recall), 6) if truth_recall else None,
        "precision": round(1 - len(extra) / len(found), 6) if found else None,
        "extra_positions": len(extra),
        "extra_fraction": round(len(extra) / len(found), 6) if found else None,
        "first_extra": [f"{c}:{p}" for c, p in sorted(extra)[:10]],
        "first_missed": [f"{c}:{p}" for c, p in sorted(k for k in truth_recall if k not in found)[:10]],
        "duplicate_ref_positions": sum(1 for rows in found.values() if len(rows) > 1),
    }

    by_bucket = {b: {"snp_rows": 0, "matching_rows": 0, "indel_rows": 0} for b in BUCKETS}
    for r in snp_rows:
        b = by_bucket[bucket(float(r["Perc_Iden"]))]
        b["snp_rows"] += 1
        b["matching_rows"] += (r["Ref_Contig"], int(r["Ref_Pos"])) in truth_all
    for r in indel_rows:
        by_bucket[bucket(float(r["Perc_Iden"]))]["indel_rows"] += 1
    out["by_identity_bucket"] = by_bucket

    blocks = covered(bed)
    ref_bases = float(header.get("Reference_Assembly_Bases", "nan").replace("NA", "nan"))
    qry_bases = float(header.get("Query_Assembly_Bases", "nan").replace("NA", "nan"))
    ref_union = union_length((r["Ref_Contig"], int(r["Ref_Start"]), int(r["Ref_End"])) for r in blocks) if blocks else 0
    qry_union = union_length((r["Query_Contig"], int(r["Query_Start"]), int(r["Query_End"])) for r in blocks) if blocks else 0
    out["bed_blocks"] = len(blocks)
    out["reference_breadth_pct"] = pct(ref_union, ref_bases) if ref_bases == ref_bases else None
    out["query_breadth_pct"] = pct(qry_union, qry_bases) if qry_bases == qry_bases else None

    if report:
        total_r, total_q = (float(x) for x in report.get("TotalBases", ("nan", "nan")))
        aligned_r, aligned_q = (float(x) for x in report.get("AlignedBases", ("nan", "nan")))
        out["dnadiff_reference_breadth_pct"] = pct(aligned_r, total_r)
        out["dnadiff_query_breadth_pct"] = pct(aligned_q, total_q)
        for side in ("reference", "query"):
            mine, theirs = out[f"{side}_breadth_pct"], out[f"dnadiff_{side}_breadth_pct"]
            out[f"{side}_breadth_delta_pp"] = round(mine - theirs, 4) if None not in (mine, theirs) else None
    if coords:
        dnadiff_pairs = {(c["ref"], c["qry"]) for c in coords}
        phraya_pairs = {(r["Ref_Contig"], r["Query_Contig"]) for r in blocks if float(r["Perc_Iden"]) >= 90}
        out["pair_recall"] = round(len(dnadiff_pairs & phraya_pairs) / len(dnadiff_pairs), 6) if dnadiff_pairs else None

    invariants = {}
    for key in ("Reference_Percent_Aligned", "Query_Percent_Aligned"):
        value = header.get(key)
        invariants[f"{key}_in_range"] = value not in (None, "NA") and 0.0 <= float(value) <= 100.0
    invariants["header_snps_match"] = header.get("SNPs") == str(len(snp_rows))
    invariants["header_indels_match"] = header.get("Indels") == str(len(indel_rows))
    for key, mine in (("Reference_Percent_Aligned", out["reference_breadth_pct"]), ("Query_Percent_Aligned", out["query_breadth_pct"])):
        value = header.get(key)
        invariants[f"{key}_matches_bed"] = (
            value not in (None, "NA") and mine is not None and abs(float(value) - mine) <= 0.01)
    out["header_invariants"] = invariants

    gate = {
        "recall": out["recall"] is not None and out["recall"] >= args.min_recall,
        "allele_concordance": out["allele_concordance"] is not None and out["allele_concordance"] >= args.min_allele_concordance,
        "extra_positions": out["extra_fraction"] is not None and out["extra_fraction"] <= args.max_extra_fraction,
        "no_duplicates": out["duplicate_ref_positions"] == 0,
        "header_invariants": all(invariants.values()),
    }
    for side in ("reference", "query"):
        delta = out.get(f"{side}_breadth_delta_pp")
        if delta is not None:
            gate[f"{side}_breadth"] = abs(delta) <= args.max_breadth_delta
    out["gate"] = gate
    out["gate_passed"] = all(gate.values())

    json.dump(out, sys.stdout, indent=2)
    sys.stdout.write("\n")
    print(f"recall {out['recall']}  allele {out['allele_concordance']}  extra {out['extra_positions']}  "
          f"dups {out['duplicate_ref_positions']}  breadth ref {out['reference_breadth_pct']} "
          f"qry {out['query_breadth_pct']}  gate {'PASS' if out['gate_passed'] else 'FAIL'}", file=sys.stderr)
    return 1 if args.gate and not out["gate_passed"] else 0


if __name__ == "__main__":
    sys.exit(main())
