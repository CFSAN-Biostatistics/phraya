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
    return out


# ---------------------------------------------------------------------------
# Placement accuracy
# ---------------------------------------------------------------------------

PAFTOOLS_CANDIDATES = [
    os.environ.get("PAFTOOLS_BIN", ""),
    "paftools.js",
    "/nfs/software/apps/micromamba/1.5.8/envs/quast-v5.2.0/bin/paftools.js",
    "/nfs/software/apps/micromamba/1.5.8/envs/minimap2-v2.28/bin/paftools.js",
]

def find_paftools() -> str | None:
    for candidate in PAFTOOLS_CANDIDATES:
        if subprocess.run(["which", candidate], capture_output=True).returncode == 0:
            return candidate
        if Path(candidate).exists():
            return candidate
    return None


SAMTOOLS_BIN = os.environ.get(
    "SAMTOOLS_BIN",
    "/nfs/software/apps/micromamba/1.5.8/envs/samtools-v1.20/bin/samtools",
)


def compute_pa_sam(sam_file: Path, paftools: str, tolerance: int = 10) -> float | None:
    """Placement accuracy from a SAM file via paftools.js mapeval."""
    try:
        sam2paf = subprocess.run(
            f"{SAMTOOLS_BIN} view -F4 {sam_file} | {paftools} sam2paf -",
            shell=True, capture_output=True, text=True, timeout=600,
        )
        if sam2paf.returncode != 0:
            print(f"WARNING: sam2paf failed for {sam_file}", file=sys.stderr)
            return None

        mapeval = subprocess.run(
            [paftools, "mapeval", "-"],
            input=sam2paf.stdout, capture_output=True, text=True, timeout=600,
        )
        if mapeval.returncode != 0:
            print(f"WARNING: mapeval failed for {sam_file}", file=sys.stderr)
            return None

        for line in mapeval.stdout.splitlines():
            if line.startswith("Q"):
                parts = line.split()
                if len(parts) >= 5 and parts[1] == str(tolerance):
                    try:
                        return float(parts[4])
                    except (ValueError, IndexError):
                        continue
        print(f"WARNING: Could not parse PA (d={tolerance}) from mapeval for {sam_file}", file=sys.stderr)
        return None
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

def load_targets_conf(run_dir: Path) -> tuple[dict, dict]:
    """Returns (genome_sizes, target_paths) dicts keyed by target id."""
    candidates = [
        run_dir.parent.parent / "scripts" / "benchmark" / "slurm" / "config" / "targets.conf",
        SCRIPT_DIR.parent / "config" / "targets.conf",
    ]
    conf_path = next((p for p in candidates if p.exists()), None)
    genome_sizes: dict[str, float] = {}
    target_paths: dict[str, str] = {}
    if conf_path is None:
        print("WARNING: targets.conf not found — genome sizes will be 0", file=sys.stderr)
        return genome_sizes, target_paths
    for line in conf_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("|")
        if len(parts) >= 4:
            tid = parts[0].strip()
            genome_sizes[tid] = float(parts[3].strip())
            target_paths[tid] = parts[1].strip()
    return genome_sizes, target_paths


# ---------------------------------------------------------------------------
# Main aggregation
# ---------------------------------------------------------------------------

def aggregate(run_dir_str: str) -> dict:
    run_dir = Path(run_dir_str)
    stream_triad_gbps = read_stream_triad(run_dir)
    genome_sizes, target_paths = load_targets_conf(run_dir)
    paftools = find_paftools()
    if paftools is None:
        print("WARNING: paftools.js not found — PA will be null for SAM-output aligners", file=sys.stderr)

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
            if total_reads == 0 and target_id in target_paths:
                fastq = data_root / target_paths[target_id] / "data" / "reads" / "reads_1_30k.fastq.gz"
                if not fastq.exists():
                    fastq = data_root / target_paths[target_id] / "data" / "reads" / "reads_1.fastq.gz"
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

            if is_phraya:
                queries_file = rep0_dir / "alignment.phraya.queries"
                if queries_file.exists():
                    print(f"Computing PA (phraya) for {aligner} {target_id}...", file=sys.stderr)
                    pa = compute_pa_phraya(queries_file, total_reads=total_reads)
            elif paftools:
                sam_file = rep0_dir / "alignment.sam"
                if not sam_file.exists():
                    # bwa-pipeline uses BAM — convert to SAM on the fly
                    bam_file = rep0_dir / "alignment.bam"
                    if bam_file.exists():
                        sam_file = bam_file  # samtools view handles both
                if sam_file.exists():
                    print(f"Computing PA (paftools) for {aligner} {target_id}...", file=sys.stderr)
                    pa = compute_pa_sam(sam_file, paftools)

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
                "pa": round(pa, 4) if pa is not None else None,
                "mcs": None,  # MAPQ calibration — not yet implemented
                "peak_rss_gb": round(mean(rss_values), 3) if any(v > 0 for v in rss_values) else None,
                "genome_size_gb": genome_sizes.get(target_id, 0.0),
                "unaligned_frac": round(mean_unaligned, 4) if mean_unaligned is not None else None,
            }
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
