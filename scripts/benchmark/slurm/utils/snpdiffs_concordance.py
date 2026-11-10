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
  6. header invariants: *_Percent_Aligned in [0, 100], header SNPs equal to SNP rows, header
     Indels (indel bases, as CSP2's SNP_Cutoff counts them) at least the indel rows,
     header breadth equal to the BED union

Truth sources (exactly one):
  --truth truth.tsv   from gen_assembly_pair.py; recall uses the callable subset when
                      --callable-only is set, precision always uses every planted SNP
  --dnadiff PREFIX    PREFIX.snps, PREFIX.1coords and PREFIX.report from dnadiff

--csp2-filter (needs --dnadiff) is the acceptance gate: both call sets go through CSP2's
`filterSNPs` (csp2_common.py) in its order — min_len/min_iden on the row's alignment, Indel and
Invalid rows dropped, heterozygous loci purged, duplicates resolved by (Ref_Aligned, Perc_Iden),
density windows, then reference/query edge trimming. dnadiff rows get their alignment length
and identity from .1coords as compileMUMmer.makeSNPCoords assigns them. Gate: both sides agree
on CSP2's SNP_Cutoff (SNPs + indels + invalid > 10000 halts a pair); when neither halts, the
surviving SNP sets must be equal. Raw recall and extras are still reported but not gated.

Prints a JSON object on stdout and a short summary on stderr; exits 1 if --gate is set
and a pass criterion fails. Stdlib only.
"""
import argparse
import json
import sys
from collections import defaultdict
from bisect import bisect_left, bisect_right

# CSP2's screening cutoff (screenSNPDiffs.py / runSNPPipeline.py): a pair with more than this
# many raw SNPs + indels + invalid sites is halted as SNP_Cutoff.
CSP2_SNP_CUTOFF = 10000
CSP2_BASES = set("ACGTacgt")

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


def csp2_rows_from_snpdiffs(variants):
    """CSP2 filter inputs from snpdiffs variant rows."""
    return [{
        "contig": v["Ref_Contig"], "pos": int(v["Ref_Pos"]), "query_base": v["Query_Base"],
        "ref_aligned": int(v["Ref_Aligned"]), "iden": float(v["Perc_Iden"]),
        "dist_ref": int(v["Dist_to_Ref_End"]), "dist_query": int(v["Dist_to_Query_End"]), "cat": v["Cat"],
    } for v in variants]


def csp2_rows_from_dnadiff(prefix: str):
    """CSP2 filter inputs from dnadiff, as compileMUMmer.parseMUmmerSNPs/makeSNPCoords build them.

    Returns (rows, raw_snps, raw_indels, raw_invalid)."""
    aligns = defaultdict(list)
    with open(prefix + ".1coords") as f:
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 13 or not fields[0].strip().isdigit():
                continue
            r1, r2, q1, q2 = (int(x) for x in fields[:4])
            aligns[(fields[11], fields[12])].append(
                (min(r1, r2), max(r1, r2), min(q1, q2), max(q1, q2), int(fields[4]), float(fields[6])))
    rows, counts = [], {"SNP": 0, "Indel": 0, "Invalid": 0}
    with open(prefix + ".snps") as f:
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 12 or not fields[0].strip().isdigit():
                continue
            ref_pos, ref_base, qry_base, qry_pos = int(fields[0]), fields[1], fields[2], int(fields[3])
            ref_len, qry_len, ref_contig, qry_contig = int(fields[6]), int(fields[7]), fields[10], fields[11]
            if not {ref_base, qry_base} <= CSP2_BASES | {"."}:
                cat = "Invalid"
            elif "." in (ref_base, qry_base):
                cat = "Indel"
            else:
                cat = "SNP"
            counts[cat] += 1
            best = max(((length, iden) for r1, r2, q1, q2, length, iden in aligns[(ref_contig, qry_contig)]
                        if r1 <= ref_pos <= r2 and q1 <= qry_pos <= q2), default=None)
            if best is None:
                continue
            rows.append({
                "contig": ref_contig, "pos": ref_pos, "query_base": qry_base,
                "ref_aligned": best[0], "iden": best[1],
                "dist_ref": min(ref_pos, ref_len - ref_pos), "dist_query": min(qry_pos, qry_len - qry_pos), "cat": cat,
            })
    return rows, counts["SNP"], counts["Indel"], counts["Invalid"]


def csp2_filter(rows, min_len, min_iden, edge, windows, max_snps):
    """CSP2 filterSNPs (csp2_common.py): returns {(contig, pos): query_base} of surviving SNPs."""
    rows = [r for r in rows if r["ref_aligned"] >= min_len and r["iden"] >= min_iden and r["cat"] == "SNP"]
    by_loc = defaultdict(list)
    for r in rows:
        by_loc[(r["contig"], r["pos"])].append(r)
    kept = {}
    for loc, group in by_loc.items():
        if len({r["query_base"] for r in group}) > 1:
            continue  # heterozygous
        kept[loc] = max(group, key=lambda r: (r["ref_aligned"], r["iden"]))  # duplicates
    positions = defaultdict(list)
    for contig, pos in kept:
        positions[contig].append(pos)
    purged = set()
    for contig, sites in positions.items():
        sites.sort()
        for window, limit in zip(windows, max_snps):
            for pos in sites:
                lo, hi = bisect_left(sites, pos - window), bisect_right(sites, pos + window)
                if hi - lo > limit:  # bedtools window -c counts the site itself
                    purged.update((contig, p) for p in sites[lo:hi])
    return {loc: r["query_base"] for loc, r in kept.items()
            if loc not in purged and r["dist_ref"] >= edge and r["dist_query"] >= edge}


def csp2_coverage_dnadiff(prefix: str, min_len: int, min_iden: float):
    """dnadiff 1-to-1 intervals CSP2 would count as covered: Ref_Aligned >= min_len and identity >= min_iden."""
    spans = []
    with open(prefix + ".1coords") as f:
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 13 and fields[0].strip().isdigit() and int(fields[4]) >= min_len and float(fields[6]) >= min_iden:
                start, end = sorted((int(fields[0]), int(fields[1])))
                spans.append((fields[11], start - 1, end))
    return spans


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
    ap.add_argument("--csp2-filter", action="store_true",
                    help="gate on CSP2 filterSNPs output equality (needs --dnadiff)")
    ap.add_argument("--csp2-min-len", type=int, default=500)
    ap.add_argument("--csp2-min-iden", type=float, default=99.0)
    ap.add_argument("--csp2-edge", type=int, default=150, help="ref_edge and query_edge")
    ap.add_argument("--csp2-dwin", default="1000,125,15")
    ap.add_argument("--csp2-wsnps", default="3,2,1")
    args = ap.parse_args()
    if args.csp2_filter and not args.dnadiff:
        ap.error("--csp2-filter needs --dnadiff")

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
    invariants["header_indels_cover_rows"] = header.get("Indels", "").isdigit() and int(header["Indels"]) >= len(indel_rows)
    for key, mine in (("Reference_Percent_Aligned", out["reference_breadth_pct"]), ("Query_Percent_Aligned", out["query_breadth_pct"])):
        value = header.get(key)
        invariants[f"{key}_matches_bed"] = (
            value not in (None, "NA") and mine is not None and abs(float(value) - mine) <= 0.01)
    out["header_invariants"] = invariants

    if args.csp2_filter:
        windows = [int(x) for x in args.csp2_dwin.split(",")]
        limits = [int(x) for x in args.csp2_wsnps.split(",")]
        filt = lambda rows: csp2_filter(rows, args.csp2_min_len, args.csp2_min_iden, args.csp2_edge, windows, limits)
        phraya = filt(csp2_rows_from_snpdiffs(variants))
        dn_rows, dn_snps, dn_indels, dn_invalid = csp2_rows_from_dnadiff(args.dnadiff)
        dnadiff = filt(dn_rows)
        shared = [k for k in phraya if k in dnadiff]
        raw_phraya = sum(int(header.get(k, "0").replace("NA", "0")) for k in ("SNPs", "Indels", "Invalid"))
        csp2 = {
            "phraya_kept": len(phraya), "dnadiff_kept": len(dnadiff), "shared": len(shared),
            "allele_match": sum(1 for k in shared if phraya[k] == dnadiff[k]),
            "phraya_only": [f"{c}:{p}" for c, p in sorted(set(phraya) - set(dnadiff))],
            "dnadiff_only": [f"{c}:{p}" for c, p in sorted(set(dnadiff) - set(phraya))],
            "phraya_raw_variants": raw_phraya, "dnadiff_raw_variants": dn_snps + dn_indels + dn_invalid,
        }
        csp2["phraya_snp_cutoff"] = raw_phraya > CSP2_SNP_CUTOFF
        csp2["dnadiff_snp_cutoff"] = csp2["dnadiff_raw_variants"] > CSP2_SNP_CUTOFF
        out["csp2"] = csp2
        # CSP2 recomputes reference coverage from BED rows passing min_len/min_iden (screenSNPDiffs),
        # so the raw union breadth above is not what it consumes.
        ref_total = float(header.get("Reference_Assembly_Bases", "nan").replace("NA", "nan"))
        phraya_spans = [(r["Ref_Contig"], int(r["Ref_Start"]), int(r["Ref_End"])) for r in covered(bed)
                        if int(r["Ref_Aligned"]) >= args.csp2_min_len and float(r["Perc_Iden"]) >= args.csp2_min_iden]
        dnadiff_spans = csp2_coverage_dnadiff(args.dnadiff, args.csp2_min_len, args.csp2_min_iden)
        csp2["phraya_reference_coverage_pct"] = pct(union_length(phraya_spans), ref_total) if phraya_spans and ref_total == ref_total else 0.0
        csp2["dnadiff_reference_coverage_pct"] = pct(union_length(dnadiff_spans), ref_total) if dnadiff_spans and ref_total == ref_total else 0.0
        csp2["reference_coverage_delta_pp"] = round(csp2["phraya_reference_coverage_pct"] - csp2["dnadiff_reference_coverage_pct"], 4)
        csp2_sets_equal = (not csp2["phraya_only"] and not csp2["dnadiff_only"]
                           and csp2["allele_match"] == len(shared))
        csp2["snps_equal"] = csp2_sets_equal
        gate = {
            "csp2_cutoff_agrees": csp2["phraya_snp_cutoff"] == csp2["dnadiff_snp_cutoff"],
            # Set equality is required unconditionally; for cutoff pairs the kept sets are raw
            # (no csp2 filter applied), so they differ and this correctly reports the disagreement.
            "csp2_snps_equal": csp2_sets_equal,
            "no_duplicates": out["duplicate_ref_positions"] == 0,
            "csp2_coverage": csp2["phraya_snp_cutoff"] and csp2["dnadiff_snp_cutoff"]
                             or abs(csp2["reference_coverage_delta_pp"]) <= args.max_breadth_delta,
            "header_invariants": all(invariants.values()),
        }
    else:
        gate = {
            "recall": out["recall"] is not None and out["recall"] >= args.min_recall,
            "allele_concordance": out["allele_concordance"] is not None and out["allele_concordance"] >= args.min_allele_concordance,
            "extra_positions": out["extra_fraction"] is not None and out["extra_fraction"] <= args.max_extra_fraction,
            "no_duplicates": out["duplicate_ref_positions"] == 0,
            "header_invariants": all(invariants.values()),
        }
    for side in ("reference", "query"):
        delta = out.get(f"{side}_breadth_delta_pp")
        if delta is not None and not args.csp2_filter:
            gate[f"{side}_breadth"] = abs(delta) <= args.max_breadth_delta
    out["gate"] = gate
    out["gate_passed"] = all(gate.values())

    json.dump(out, sys.stdout, indent=2)
    sys.stdout.write("\n")
    if args.csp2_filter:
        c = out["csp2"]
        print(f"csp2 kept phraya {c['phraya_kept']} dnadiff {c['dnadiff_kept']} shared {c['shared']}  "
              f"phraya-only {len(c['phraya_only'])} dnadiff-only {len(c['dnadiff_only'])}  cutoff "
              f"{c['phraya_snp_cutoff']}/{c['dnadiff_snp_cutoff']}  gate {'PASS' if out['gate_passed'] else 'FAIL'}",
              file=sys.stderr)
    else:
        print(f"recall {out['recall']}  allele {out['allele_concordance']}  extra {out['extra_positions']}  "
              f"dups {out['duplicate_ref_positions']}  breadth ref {out['reference_breadth_pct']} "
              f"qry {out['query_breadth_pct']}  gate {'PASS' if out['gate_passed'] else 'FAIL'}", file=sys.stderr)
    return 1 if args.gate and not out["gate_passed"] else 0


if __name__ == "__main__":
    sys.exit(main())
