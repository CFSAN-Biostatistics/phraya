#!/usr/bin/env python3
"""
Aggregate per-replicate timing results into score.py input format.

Scans results/{run_id}/{target}/{aligner}/rep_*/timing.json
Computes mean wall_time_s and peak_rss_gb across replicates
Computes placement accuracy (PA) from SAM files (paftools.js) or .phraya.queries (wgsim)
Counts reads from timing.txt (captured at run time) or FASTQ files as fallback
Outputs JSON matching score.py schema.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from statistics import mean, stdev
from typing import Optional

SCRIPT_DIR = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# Stream Triad
# ---------------------------------------------------------------------------

def read_stream_triad(run_dir: Path) -> float:
    stream_file = run_dir / "stream_triad.txt"
    if not stream_file.exists():
        sys.exit(f"ERROR: {stream_file} not found — run STREAM Triad characterisation first")

    values = []
    for line in stream_file.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        # Accept either single float ("125385.6") or "<label> <float>" format
        parts = line.split()
        for part in reversed(parts):  # last token most likely to be the number
            try:
                values.append(float(part))
                break
            except ValueError:
                continue

    if not values:
        sys.exit(f"ERROR: No valid Triad values in {stream_file}")
    return mean(values)


# ---------------------------------------------------------------------------
# Timing.txt / timing.json parsing
# ---------------------------------------------------------------------------

def parse_timing_txt(path: Path) -> dict:
    """Parse key=value pairs from timing.txt into a dict."""
    result = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if "=" in line:
            k, _, v = line.partition("=")
            result[k.strip()] = v.strip()
    return result


def load_rep_timing(rep_dir: Path) -> Optional[dict]:
    """Load timing data from rep_dir, returning a unified dict or None."""
    # Prefer timing.json (written by benchmark.slurm); fall back to timing.txt
    json_file = rep_dir / "timing.json"
    txt_file = rep_dir / "timing.txt"

    raw = {}
    if json_file.exists():
        try:
            raw = json.loads(json_file.read_text())
        except json.JSONDecodeError as e:
            print(f"WARNING: {json_file}: {e}", file=sys.stderr)
    elif txt_file.exists():
        raw = parse_timing_txt(txt_file)
    else:
        return None

    if "error" in raw:
        print(f"WARNING: {rep_dir} has error: {raw['error']}", file=sys.stderr)
        return None

    # Normalise field names (wrappers write wall_seconds; score.py wants wall_time_s)
    out = {}
    out["wall_time_s"] = float(raw.get("wall_time_s") or raw.get("wall_seconds") or 0)
    out["peak_rss_gb"] = float(raw.get("peak_rss_gb") or 0)
    out["total_reads"] = int(raw.get("total_reads") or 0)
    out["n_aligned"] = int(raw.get("n_aligned") or 0)
    out["unaligned_frac"] = float(raw.get("unaligned_frac") or 0)
    out["n_variants"] = raw.get("n_variants")
    out["alphabet"] = raw.get("alphabet", "dna")  # older runs predate this field
    return out


# ---------------------------------------------------------------------------
# Placement accuracy
# ---------------------------------------------------------------------------

K8_BIN = os.environ.get("K8_BIN", "/nfs/software/apps/micromamba/1.5.8/envs/minimap2-v2.28/bin/k8")
PAFTOOLS_CANDIDATES = [
    os.environ.get("PAFTOOLS_BIN", ""),
    "/nfs/software/apps/micromamba/1.5.8/envs/minimap2-v2.28/bin/paftools.js",
    "/nfs/software/apps/micromamba/1.5.8/envs/quast-v5.2.0/bin/paftools.js",
]

def find_paftools() -> str | None:
    """Returns path to paftools.js if k8 is available to run it."""
    if not Path(K8_BIN).exists():
        return None
    for candidate in PAFTOOLS_CANDIDATES:
        if candidate and Path(candidate).exists():
            return candidate
    return None


SAMTOOLS_BIN = os.environ.get(
    "SAMTOOLS_BIN",
    "/nfs/software/apps/micromamba/1.5.8/envs/samtools-v1.20/bin/samtools",
)


def compute_pa_sam(sam_file: Path, paftools: str, tolerance: int = 10) -> float | None:
    """Placement accuracy from a SAM/BAM file using wgsim read-name encoding.

    Uses sam_accuracy.py which parses wgsim-format QNAME fields directly.
    paftools.js mapeval does not support wgsim read names (only dwgsim format).
    """
    helper = SCRIPT_DIR / "sam_accuracy.py"
    if not helper.exists():
        print(f"WARNING: {helper} not found — PA will be null for SAM aligners", file=sys.stderr)
        return None
    python3 = os.environ.get("PYTHON3_BIN", "python3")
    try:
        result = subprocess.run(
            [python3, str(helper), str(sam_file),
             f"--tolerance={tolerance}", f"--samtools={SAMTOOLS_BIN}"],
            capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            print(f"WARNING: sam_accuracy.py failed for {sam_file}: {result.stderr[:200]}", file=sys.stderr)
            return None
        parts = result.stdout.strip().split("\t")
        return float(parts[0])
    except subprocess.TimeoutExpired:
        print(f"WARNING: PA computation timed out for {sam_file}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"WARNING: PA failed for {sam_file}: {e}", file=sys.stderr)
        return None


def compute_pa_phraya(queries_file: Path, total_reads: int = 0, tolerance: int = 10) -> float | None:
    """Placement accuracy from .phraya.queries using wgsim read-name encoding."""
    helper = SCRIPT_DIR / "phraya_accuracy.py"
    if not helper.exists():
        print(f"WARNING: {helper} not found — PA will be null for phraya", file=sys.stderr)
        return None
    try:
        python3 = os.environ.get("PYTHON3_BIN", "python3")
        cmd = [python3, str(helper), str(queries_file), f"--tolerance={tolerance}"]
        if total_reads > 0:
            cmd.append(f"--total-reads={total_reads}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            print(f"WARNING: phraya_accuracy.py failed: {result.stderr[:200]}", file=sys.stderr)
            return None
        parts = result.stdout.strip().split("\t")
        return float(parts[0])
    except Exception as e:
        print(f"WARNING: PA (phraya) failed for {queries_file}: {e}", file=sys.stderr)
        return None


def compute_pa_protein(queries_file: Path, truth_file: Path, total_queries: int = 0, tolerance: int = 10) -> float | None:
    """Placement accuracy for phraya protein-mode alignment, using
    gen_synthetic_protein.py's own truth.tsv sidecar (ADR-0013) instead of
    wgsim/dwgsim read-name encoding — see phraya_accuracy_protein.py."""
    helper = SCRIPT_DIR / "phraya_accuracy_protein.py"
    if not helper.exists():
        print(f"WARNING: {helper} not found — PA will be null for phraya-protein", file=sys.stderr)
        return None
    if not truth_file.exists():
        print(f"WARNING: {truth_file} not found — PA will be null for phraya-protein", file=sys.stderr)
        return None
    try:
        python3 = os.environ.get("PYTHON3_BIN", "python3")
        cmd = [python3, str(helper), str(queries_file), f"--truth={truth_file}", f"--tolerance={tolerance}"]
        if total_queries > 0:
            cmd.append(f"--total-queries={total_queries}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            print(f"WARNING: phraya_accuracy_protein.py failed: {result.stderr[:200]}", file=sys.stderr)
            return None
        parts = result.stdout.strip().split("\t")
        return float(parts[0])
    except Exception as e:
        print(f"WARNING: PA (phraya-protein) failed for {queries_file}: {e}", file=sys.stderr)
        return None


def compute_iec(phraya_file: Path, truth_file: Path) -> float | None:
    """Indel Event Concordance (see BENCHMARK_EXPANSION.md) — how often a
    read/query's reported CIGAR consolidates its true indel events into the
    right number of I/D ops. Runs `phraya filter --format tsv` and joins
    against a gen_synthetic*.py truth sidecar via compute_indel_recovery.py.
    Returns None (not 0.0) when there's nothing to measure — a target with no
    indel events in its truth data isn't "0% concordant", it's inapplicable."""
    helper = SCRIPT_DIR.parent.parent / "local" / "compute_indel_recovery.py"
    if not helper.exists() or not truth_file.exists() or not phraya_file.exists():
        return None
    try:
        phraya_bin = os.environ.get("PHRAYA_ROOT", "") 
        phraya_bin = str(Path(phraya_bin) / "target" / "release" / "phraya") if phraya_bin else "phraya"
        filter_proc = subprocess.run(
            [phraya_bin, "filter", str(phraya_file), "--format", "tsv"],
            capture_output=True, text=True, timeout=600,
        )
        if filter_proc.returncode != 0:
            print(f"WARNING: phraya filter failed for IEC on {phraya_file}: {filter_proc.stderr[:200]}", file=sys.stderr)
            return None
        python3 = os.environ.get("PYTHON3_BIN", "python3")
        iec_proc = subprocess.run(
            [python3, str(helper), f"--truth={truth_file}"],
            input=filter_proc.stdout, capture_output=True, text=True, timeout=300,
        )
        if iec_proc.returncode != 0:
            print(f"WARNING: compute_indel_recovery.py failed: {iec_proc.stderr[:200]}", file=sys.stderr)
            return None
        for line in iec_proc.stdout.splitlines():
            if line.startswith("iec\t"):
                value = line.split("\t", 1)[1]
                return None if value == "nan" else float(value)
        return None
    except Exception as e:
        print(f"WARNING: IEC computation failed for {phraya_file}: {e}", file=sys.stderr)
        return None


def gap_model_from_aligner(aligner: str) -> str | None:
    """Parses the `--gap-model` axis (ADR-0014) out of an aligner-variant name
    (e.g. `phraya-sensitive-affine` -> "affine"). `None` when the aligner name
    doesn't encode a gap model (non-phraya tools, or `phraya-sensitive` bare —
    this script can't observe phraya's own runtime default, so it doesn't
    guess)."""
    if aligner.endswith("-affine"):
        return "affine"
    if aligner.endswith("-linear"):
        return "linear"
    return None


# ---------------------------------------------------------------------------
# Read counting
# ---------------------------------------------------------------------------

def count_fastq_reads(path: Path) -> int | None:
    try:
        result = subprocess.run(
            f"zcat {path} | wc -l",
            shell=True, capture_output=True, text=True, timeout=300,
        )
        if result.returncode != 0:
            return None
        return int(result.stdout.strip()) // 4
    except Exception:
        return None


# ---------------------------------------------------------------------------
# targets.conf loader
# ---------------------------------------------------------------------------

def load_targets_conf(run_dir: Path, filename: str = "targets.conf") -> tuple[dict, dict]:
    """Returns (size_metric, target_paths) dicts keyed by target id. `size_metric`
    is the 4th pipe-delimited column, whatever it means for `filename`
    (genome GB for targets.conf, proteome residue count for
    targets_protein.conf — both share the same 4-field shape)."""
    candidates = [
        run_dir.parent.parent / "scripts" / "benchmark" / "slurm" / "config" / filename,
        SCRIPT_DIR.parent / "config" / filename,
    ]
    conf_path = next((p for p in candidates if p.exists()), None)
    size_metric: dict[str, float] = {}
    target_paths: dict[str, str] = {}
    if conf_path is None:
        print(f"WARNING: {filename} not found — sizes will be 0", file=sys.stderr)
        return size_metric, target_paths
    for line in conf_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("|")
        if len(parts) >= 4:
            tid = parts[0].strip()
            size_metric[tid] = float(parts[3].strip())
            target_paths[tid] = parts[1].strip()
    return size_metric, target_paths


def load_targets_protein_conf(run_dir: Path) -> tuple[dict, dict]:
    """Returns (proteome_residues, target_paths) from targets_protein.conf."""
    return load_targets_conf(run_dir, filename="targets_protein.conf")


# ---------------------------------------------------------------------------
# Main aggregation
# ---------------------------------------------------------------------------

def aggregate(run_dir_str: str) -> dict:
    run_dir = Path(run_dir_str)
    stream_triad_gbps = read_stream_triad(run_dir)
    genome_sizes, dna_target_paths = load_targets_conf(run_dir)
    proteome_residues, protein_target_paths = load_targets_protein_conf(run_dir)
    # PA for SAM/BAM aligners uses sam_accuracy.py (wgsim read-name parsing),
    # not paftools.js mapeval (which expects dwgsim format, not wgsim format).
    paftools = "unused"  # kept for compatibility with call site below

    data_root = Path.home() / "data-commons" / "test" / "benchmarking" / "alignment"

    # results[aligner][target_id] = list of per-rep dicts
    results: dict[str, dict[str, list]] = {}

    for target_dir in sorted(run_dir.iterdir()):
        if not target_dir.is_dir() or target_dir.name.startswith("."):
            continue
        target_id = target_dir.name

        for aligner_dir in sorted(target_dir.iterdir()):
            if not aligner_dir.is_dir():
                continue
            aligner = aligner_dir.name
            results.setdefault(aligner, {}).setdefault(target_id, [])

            for rep_dir in sorted(aligner_dir.iterdir()):
                if not rep_dir.is_dir() or not rep_dir.name.startswith("rep_"):
                    continue
                timing = load_rep_timing(rep_dir)
                if timing is not None:
                    results[aligner][target_id].append(timing)
                else:
                    print(f"WARNING: No timing in {rep_dir}", file=sys.stderr)

    output = {
        "platform": {
            "stream_triad_gbps": round(stream_triad_gbps / 1000, 2),  # MB/s → GB/s
            "threads": 8,
            "cpu_model": "unknown",
        },
        "aligners": [],
    }

    for aligner, targets in results.items():
        aligner_entry = {"name": aligner, "version": "unknown", "targets": []}

        for target_id, reps in targets.items():
            if not reps:
                print(f"WARNING: No valid reps for {aligner} {target_id}", file=sys.stderr)
                continue

            wall_times = [r["wall_time_s"] for r in reps]
            rss_values = [r["peak_rss_gb"] for r in reps]
            mean_wall = mean(wall_times)
            if len(wall_times) > 1:
                cv = (stdev(wall_times) / mean_wall * 100) if mean_wall > 0 else 0.0
                if cv > 5.0:
                    print(f"WARNING: {aligner} {target_id} CV={cv:.1f}% (>5%)", file=sys.stderr)

            # Read count: prefer what wrappers captured, fall back to FASTQ
            total_reads = max((r["total_reads"] for r in reps), default=0)
            if total_reads == 0 and target_id in dna_target_paths:
                fastq = data_root / dna_target_paths[target_id] / "data" / "reads" / "reads_1_30k.fastq.gz"
                if not fastq.exists():
                    fastq = data_root / dna_target_paths[target_id] / "data" / "reads" / "reads_1.fastq.gz"
                if fastq.exists():
                    count = count_fastq_reads(fastq)
                    if count is not None:
                        total_reads = count

            # Unaligned fraction: mean across reps
            unaligned_fracs = [r["unaligned_frac"] for r in reps if r["unaligned_frac"] > 0]
            mean_unaligned = mean(unaligned_fracs) if unaligned_fracs else None

            # Placement accuracy (use rep_0)
            rep0_dir = run_dir / target_id / aligner / "rep_0"
            pa = None
            is_phraya = aligner.startswith("phraya")
            alphabet = reps[0].get("alphabet", "dna")
            gap_model = gap_model_from_aligner(aligner)

            if is_phraya and alphabet == "protein":
                queries_file = rep0_dir / "alignment.phraya.queries"
                truth_file = data_root / protein_target_paths.get(target_id, "") / "data" / "protein" / "queries.fasta.truth.tsv"
                if queries_file.exists():
                    print(f"Computing PA (phraya-protein) for {aligner} {target_id}...", file=sys.stderr)
                    pa = compute_pa_protein(queries_file, truth_file, total_queries=total_reads)
            elif is_phraya:
                queries_file = rep0_dir / "alignment.phraya.queries"
                if queries_file.exists():
                    print(f"Computing PA (phraya) for {aligner} {target_id}...", file=sys.stderr)
                    pa = compute_pa_phraya(queries_file, total_reads=total_reads)
            elif alphabet == "protein":
                # DIAMOND/BLASTP write tabular (not SAM) output — no PA helper
                # implemented yet for that format. Timing/RSS/unaligned_frac are
                # still computed above; pa stays None rather than a fabricated 0.0.
                pass
            else:
                sam_file = rep0_dir / "alignment.sam"
                if not sam_file.exists():
                    bam_file = rep0_dir / "alignment.bam"
                    if bam_file.exists():
                        sam_file = bam_file
                if sam_file.exists():
                    print(f"Computing PA (sam) for {aligner} {target_id}...", file=sys.stderr)
                    pa = compute_pa_sam(sam_file, paftools)

            # Indel Event Concordance (ADR-0014): phraya-only (needs a .phraya
            # file to `filter --format tsv`), and only where a truth sidecar
            # with indel events exists — protein targets always have one
            # (gen_synthetic_protein.py), DNA SLURM targets don't (dwgsim
            # truth isn't in our sidecar format — see BENCHMARK_EXPANSION.md).
            iec = None
            if is_phraya and alphabet == "protein":
                phraya_file = rep0_dir / "alignment.phraya"
                truth_file = data_root / protein_target_paths.get(target_id, "") / "data" / "protein" / "queries.fasta.truth.tsv"
                iec = compute_iec(phraya_file, truth_file)

            # Variant count for bwa-pipeline (informational)
            n_variants = None
            for r in reps:
                if r.get("n_variants"):
                    try:
                        n_variants = int(r["n_variants"])
                    except (ValueError, TypeError):
                        pass
                    break

            entry = {
                "id": target_id,
                "reads": total_reads,
                "wall_time_s": round(mean_wall, 2),
                "threads": 8,
                "alphabet": alphabet,
                "pa": round(pa, 4) if pa is not None else None,
                "mcs": 0.0,  # MAPQ calibration — not yet implemented; 0.0 so score.py can run
                "peak_rss_gb": round(mean(rss_values), 3) if any(v > 0 for v in rss_values) else None,
                "unaligned_frac": round(mean_unaligned, 4) if mean_unaligned is not None else None,
            }
            if alphabet == "protein":
                residues = proteome_residues.get(target_id, 0.0)
                entry["proteome_residues"] = residues
                # RNT (Residue-Normalized Throughput): BNT's formula
                # (queries / (wall × threads × triad_GBps)) with "reads"
                # swapped for "queries" — same units/interpretation, applied to
                # the axis meaningful for protein workloads (see
                # BENCHMARK_EXPANSION.md; genome-GB-normalized BNT/MEI don't
                # apply to a proteome-scale workload).
                denom = mean_wall * 8 * stream_triad_gbps / 1000  # GB/s
                entry["rnt"] = round(total_reads / denom, 4) if denom > 0 else None
            else:
                entry["genome_size_gb"] = genome_sizes.get(target_id, 0.0)
            if gap_model is not None:
                entry["gap_model"] = gap_model
            if iec is not None:
                entry["iec"] = round(iec, 4)
            if n_variants is not None:
                entry["n_variants"] = n_variants

            aligner_entry["targets"].append(entry)

        if aligner_entry["targets"]:
            output["aligners"].append(aligner_entry)

    return output


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: aggregate_results.py <run_dir>", file=sys.stderr)
        sys.exit(1)
    result = aggregate(sys.argv[1])
    print(json.dumps(result, indent=2))
