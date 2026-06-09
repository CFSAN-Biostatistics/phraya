#!/usr/bin/env python3
"""
Aggregate per-replicate timing results into score.py input format.

Scans results/{run_id}/{target}/{aligner}/rep_*/timing.json
Computes mean±std across 3 replicates per (target, aligner)
Computes placement accuracy (PA) from SAM files via paftools.js
Counts reads from FASTQ files
Outputs JSON matching score.py schema.
"""
import json
import sys
import subprocess
import re
from pathlib import Path
from statistics import mean, stdev


def count_fastq_reads(fastq_path):
    """Count reads in FASTQ file (4 lines per record)."""
    try:
        result = subprocess.run(
            f"zcat {fastq_path} | wc -l",
            shell=True, capture_output=True, text=True, timeout=300
        )
        if result.returncode != 0:
            return None
        line_count = int(result.stdout.strip())
        return line_count // 4
    except (subprocess.TimeoutExpired, ValueError) as e:
        print(f"WARNING: Failed to count reads in {fastq_path}: {e}", file=sys.stderr)
        return None


def compute_placement_accuracy(sam_file, target_id):
    """
    Compute placement accuracy from SAM file using paftools.js.

    Returns PA at d=10bp, or None if computation fails.
    """
    try:
        # Convert SAM → PAF
        sam2paf = subprocess.run(
            f"samtools view -F4 {sam_file} | paftools.js sam2paf -",
            shell=True, capture_output=True, text=True, timeout=300
        )
        if sam2paf.returncode != 0:
            print(f"WARNING: sam2paf failed for {sam_file}", file=sys.stderr)
            return None

        # Run mapeval
        mapeval = subprocess.run(
            ["paftools.js", "mapeval", "-"],
            input=sam2paf.stdout, capture_output=True, text=True, timeout=300
        )
        if mapeval.returncode != 0:
            print(f"WARNING: mapeval failed for {sam_file}", file=sys.stderr)
            return None

        # Parse mapeval output for PA at d=10bp
        # Expected format: lines like "Q 10    12345   11234   0.9123"
        # Where columns are: type, distance, total, correct, accuracy
        for line in mapeval.stdout.splitlines():
            if line.startswith("Q") and "10" in line:
                parts = line.split()
                if len(parts) >= 5 and parts[1] == "10":
                    try:
                        accuracy = float(parts[4])
                        return accuracy
                    except (ValueError, IndexError):
                        continue

        print(f"WARNING: Could not parse PA from mapeval output for {sam_file}", file=sys.stderr)
        return None

    except subprocess.TimeoutExpired:
        print(f"WARNING: PA computation timed out for {sam_file}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"WARNING: PA computation failed for {sam_file}: {e}", file=sys.stderr)
        return None


def aggregate(run_dir):
    """Aggregate timing results from run directory."""
    run_dir = Path(run_dir)

    # Read STREAM Triad (average across nodes)
    stream_file = run_dir / "stream_triad.txt"
    if not stream_file.exists():
        print(f"ERROR: {stream_file} not found", file=sys.stderr)
        print("Run STREAM Triad characterization first", file=sys.stderr)
        sys.exit(1)

    triad_values = []
    for line in stream_file.read_text().splitlines():
        if line.strip():
            parts = line.split()
            if len(parts) >= 2:
                try:
                    triad_values.append(float(parts[1]))
                except ValueError:
                    continue

    if not triad_values:
        print(f"ERROR: No valid Triad values in {stream_file}", file=sys.stderr)
        sys.exit(1)

    stream_triad_gbps = mean(triad_values)

    # Scan results directory for timing.json files
    results = {}  # aligner → {target_id → [timing_data]}

    for target_dir in run_dir.iterdir():
        if not target_dir.is_dir() or target_dir.name.startswith('.') or target_dir.name == 'stream_triad.txt':
            continue
        target_id = target_dir.name

        for aligner_dir in target_dir.iterdir():
            if not aligner_dir.is_dir():
                continue
            aligner = aligner_dir.name

            if aligner not in results:
                results[aligner] = {}
            if target_id not in results[aligner]:
                results[aligner][target_id] = []

            for rep_dir in aligner_dir.iterdir():
                if not rep_dir.is_dir() or not rep_dir.name.startswith('rep_'):
                    continue

                timing_file = rep_dir / "timing.json"
                if not timing_file.exists():
                    print(f"WARNING: {timing_file} not found", file=sys.stderr)
                    continue

                try:
                    timing = json.loads(timing_file.read_text())
                    if "error" not in timing:
                        results[aligner][target_id].append(timing)
                    else:
                        print(f"WARNING: {timing_file} contains error: {timing['error']}", file=sys.stderr)
                except (json.JSONDecodeError, KeyError) as e:
                    print(f"WARNING: Failed to parse {timing_file}: {e}", file=sys.stderr)
                    continue

    # Load targets.conf to get genome sizes
    # Assume targets.conf is in ../config/ relative to run_dir
    targets_conf_path = run_dir.parent.parent / "scripts" / "benchmark" / "slurm" / "config" / "targets.conf"
    genome_sizes = {}
    if targets_conf_path.exists():
        for line in targets_conf_path.read_text().splitlines():
            if line.strip() and not line.startswith('#'):
                parts = line.split('|')
                if len(parts) >= 4:
                    target_id = parts[0].strip()
                    genome_size_gb = float(parts[3].strip())
                    genome_sizes[target_id] = genome_size_gb
    else:
        print(f"WARNING: targets.conf not found at {targets_conf_path}", file=sys.stderr)

    # Aggregate: compute mean wall_time_s and peak_rss_gb per aligner per target
    output = {
        "platform": {
            "stream_triad_gbps": round(stream_triad_gbps, 1),
            "threads": 8,  # From global.json
            "cpu_model": "unknown",  # TODO: capture from SLURM or /proc/cpuinfo
        },
        "aligners": []
    }

    for aligner, targets in results.items():
        aligner_entry = {
            "name": aligner,
            "version": "unknown",  # TODO: capture from module show or --version
            "targets": []
        }

        for target_id, timings in targets.items():
            if not timings:
                print(f"WARNING: No valid timings for {aligner} {target_id}", file=sys.stderr)
                continue

            wall_times = [t["wall_time_s"] for t in timings]
            rss_values = [t["peak_rss_gb"] for t in timings]

            # Compute mean and warn if CV > 5%
            mean_wall = mean(wall_times)
            if len(wall_times) > 1:
                cv = (stdev(wall_times) / mean_wall * 100) if mean_wall > 0 else 0.0
                if cv > 5.0:
                    print(f"WARNING: {aligner} {target_id} has CV={cv:.1f}% (>5%)", file=sys.stderr)

            # Compute PA from first replicate's SAM file
            pa = 0.0
            first_rep_dir = run_dir / target_id / aligner / "rep_0"
            sam_file = first_rep_dir / "alignment.sam"
            if sam_file.exists():
                print(f"Computing PA for {aligner} {target_id}...", file=sys.stderr)
                pa_result = compute_placement_accuracy(str(sam_file), target_id)
                if pa_result is not None:
                    pa = pa_result

            # Count reads from FASTQ (once per target, not per aligner)
            read_count = "unknown"
            # Try to infer data path from run_dir structure
            # This is fragile - better approach would be to pass DATA_ROOT as arg
            # For now, assume standard layout
            data_root = Path.home() / "data-commons" / "test" / "benchmarking" / "alignment"

            # Find target path from targets.conf
            target_path = None
            if targets_conf_path.exists():
                for line in targets_conf_path.read_text().splitlines():
                    if line.strip() and not line.startswith('#'):
                        parts = line.split('|')
                        if len(parts) >= 2 and parts[0].strip() == target_id:
                            target_path = parts[1].strip()
                            break

            if target_path:
                fastq_path = data_root / target_path / "data" / "reads" / "reads_1.fastq.gz"
                if fastq_path.exists():
                    print(f"Counting reads for {target_id}...", file=sys.stderr)
                    count = count_fastq_reads(str(fastq_path))
                    if count is not None:
                        read_count = count

            # Get genome size from targets.conf
            genome_size_gb = genome_sizes.get(target_id, 0.0)
            if genome_size_gb == 0.0:
                print(f"WARNING: genome_size_gb not set for {target_id} in targets.conf", file=sys.stderr)

            aligner_entry["targets"].append({
                "id": target_id,
                "reads": read_count,
                "wall_time_s": round(mean_wall, 2),
                "threads": 8,
                "pa": round(pa, 4) if isinstance(pa, float) else 0.0,
                "mcs": 0.0,  # TODO: MAPQ calibration score (requires MAPQ-stratified PA)
                "peak_rss_gb": round(mean(rss_values), 2),
                "genome_size_gb": genome_size_gb,
            })

        if aligner_entry["targets"]:
            output["aligners"].append(aligner_entry)

    return output


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: aggregate_results.py <run_dir>", file=sys.stderr)
        sys.exit(1)

    result = aggregate(sys.argv[1])
    print(json.dumps(result, indent=2))
