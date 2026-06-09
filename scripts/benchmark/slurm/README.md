# SLURM-Based HPC Aligner Benchmark

Production-grade SLURM benchmarking suite for comparing BWA-MEM2 and minimap2 across 16 benchmark targets with proper hardware normalization (STREAM Triad), 3-replicate averaging, and score.py-compatible output.

## Quick Start

```bash
cd ~/phraya/scripts/benchmark/slurm

# Dry run (preview tasks without submitting)
./run_benchmark.sh --dry-run

# Run on default targets (small + medium: T1-T7b, 9 targets)
./run_benchmark.sh

# Include large targets (T8b 17Gb, T8c 4.3Gb)
./run_benchmark.sh --large

# Use custom targets file
./run_benchmark.sh --targets my_targets.conf
```

## Prerequisites

### Data
- Benchmark dataset pre-staged at `~/data-commons/test/benchmarking/alignment/`
- Run `acquire.py --all` to download references and simulate reads

### HPC Environment
- SLURM workload manager
- Module system with `bwa-mem2`, `minimap2`, `samtools`
- Adjust module names in `config/global.env` if needed

### Phraya (Phase 2)
- Phraya integration **blocked on issue #160** (batch-mode feature)
- Phase 1 benchmarks BWA-MEM2 vs minimap2 only

## Directory Structure

```
phraya/scripts/benchmark/slurm/
├── run_benchmark.sh          # Main orchestrator (user entry point)
├── benchmark.slurm           # SLURM array job script
├── stream.slurm              # STREAM Triad platform characterization
├── config/
│   ├── targets.conf          # Target list (default: small+medium)
│   └── global.env            # Environment variables
├── wrappers/
│   ├── bwa-mem2.sh           # BWA-MEM2 invocation + timing
│   └── minimap2.sh           # minimap2 invocation + timing
└── utils/
    ├── nodelist.sh           # Query sinfo for node rotation
    ├── parse_time.py         # Parse /usr/bin/time -v → JSON
    ├── aggregate_results.py  # Collect replicates → score.py format
    └── aggregate.slurm       # Aggregation SLURM job

results/
└── run_YYYYMMDD_HHMMSS/
    ├── {target}/{aligner}/rep_{0,1,2}/
    │   ├── alignment.sam
    │   ├── timing.txt
    │   └── timing.json
    ├── stream_triad.txt
    └── results.json          # Aggregated for score.py
```

## Architecture

### SLURM Array Jobs
- Single array job with 3D indexing: `task_id → (target_idx, aligner_idx, replicate_idx)`
- Array size: `N_targets × 2_aligners × 3_replicates`
- Example: 9 targets × 2 aligners × 3 reps = 54 tasks

### Node Rotation
- Each replicate runs on different node → averages hardware quirks
- Pre-flight script queries `sinfo` for idle/mixed nodes
- Array task maps `replicate_idx mod 3 → nodelist[idx]`

### Timing Capture
- Wrap aligner invocations with `/usr/bin/time -v 2> timing.txt`
- Parse for "Elapsed (wall clock) time" and "Maximum resident set size"
- Convert to JSON: `{"wall_time_s": float, "peak_rss_gb": float}`

### Hardware Normalization
- STREAM Triad: one-time platform characterization per node type
- Measures sustained main-memory bandwidth (GB/s)
- Used by score.py for BNT (Bandwidth-Normalized Throughput) calculation

## Configuration

### Targets (`config/targets.conf`)

```bash
# Format: TARGET_ID|ORGANISM_PATH|SIZE_CLASS|GENOME_SIZE_GB
T3|mycobacterium_tuberculosis/h37rv|small|0.00441     # 4.4 Mb
T4|staphylococcus_aureus/mrsa252|small|0.00290        # 2.9 Mb
T5|plasmodium_falciparum/3d7|small|0.0233             # 23 Mb
...
T1|homo_sapiens/chr1|medium|0.249                     # 249 Mb
T2|gallus_gallus/grcg7b|medium|1.05                   # 1 Gb
...
# T8b|triticum_aestivum/hexaploid|large|14.5          # 17 Gb (excluded by default)
```

**Note**: The `GENOME_SIZE_GB` field is used for MEI (Memory Efficiency Index) calculation. Values should be canonical genome sizes from reference assemblies. Currently set to 0.0 placeholders — populate with an agent reading manifest.json files.

### Environment (`config/global.env`)

```bash
DATA_ROOT="${DATA_ROOT:-$HOME/data-commons/test/benchmarking/alignment}"
PHRAYA_ROOT="${PHRAYA_ROOT:-$HOME/phraya}"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/phraya/results/benchmark}"

THREADS=8
REPLICATES=3

MODULE_BWA="bwa-mem2"
MODULE_MINIMAP="minimap2"
MODULE_SAMTOOLS="samtools"
```

## Workflow

### 1. STREAM Triad (once per run)
- Characterizes platform memory bandwidth
- Downloads STREAM source to `~/.cache/stream/`
- Compiles with `-O3 -march=native -fopenmp`
- Runs Triad kernel, extracts GB/s
- Output: `stream_triad.txt` (one line per node)

### 2. Index Building (once per target)
- BWA-MEM2: `bwa-mem2 index reference.fasta`
- minimap2: `minimap2 -d reference.mmi reference.fasta`
- flock-protected (one index build per target across all replicates)

### 3. Alignment (array job)
- Each task decodes `$SLURM_ARRAY_TASK_ID` → (target, aligner, replicate)
- Calls wrapper script: `wrappers/$ALIGNER.sh $REF $R1 $R2 $OUT_DIR $THREADS`
- Wrapper uses `/usr/bin/time -v` for timing
- Parses timing.txt → timing.json

### 4. Aggregation (after all tasks complete)
- Scans `results/{run}/{target}/{aligner}/rep_*/timing.json`
- Computes mean±std across 3 replicates
- Warns if CV > 5% (measurement noise)
- **Computes placement accuracy (PA):**
  - Converts first replicate's SAM → PAF: `samtools view -F4 | paftools.js sam2paf`
  - Evaluates with `paftools.js mapeval` (reads true coords from dwgsim names)
  - Extracts PA at d=10bp
- **Counts reads from FASTQ:**
  - Runs `zcat reads_1.fastq.gz | wc -l` to get exact read count
  - Runs on compute node (aggregate.slurm), not login node
- Reads genome sizes from targets.conf (4th field)
- Outputs JSON matching score.py schema

**Requirements:**
- `paftools.js` (from minimap2 suite) in PATH
- `samtools` module loaded
- Aggregation job has sufficient time (0:30:00 default, may need longer for large targets)

### 5. Scoring (manual)
```bash
python ~/data-commons/test/benchmarking/alignment/score.py results/run_*/results.json --sensitivity
```

## Output Format

### `results.json` (score.py input)
```json
{
  "platform": {
    "stream_triad_gbps": 80.0,
    "threads": 8
  },
  "aligners": [
    {
      "name": "bwa-mem2",
      "version": "2.2.1",
      "targets": [
        {
          "id": "T3",
          "reads": "unknown",
          "wall_time_s": 3.1,
          "threads": 8,
          "pa": 0.0,
          "mcs": 0.0,
          "peak_rss_gb": 0.8,
          "genome_size_gb": 0.0
        }
      ]
    }
  ]
}
```

**Computed automatically:**
- ✅ `reads`: counted from FASTQ via `zcat | wc -l`
- ✅ `pa`: placement accuracy at d=10bp via `paftools.js mapeval`
- ✅ `genome_size_gb`: read from targets.conf (4th field)

**Still placeholder:**
- ⚠️ `mcs`: MAPQ calibration score (requires MAPQ-stratified PA, deferred to Phase 2)

## Troubleshooting

### Node Discovery Fails
```
ERROR: No available nodes in partition batch
```
**Solution**: Check `sinfo` output, adjust `SLURM_PARTITION` in `utils/nodelist.sh`, or skip node rotation (hardcode nodelist).

### Module Load Fails
```
WARNING: Failed to load bwa-mem2
```
**Solution**: Run `module avail` to find correct module names, update `config/global.env`.

### Input File Not Found
```
ERROR: Input file not found: .../reference.fasta
```
**Solution**: Run `python ~/data-commons/test/benchmarking/alignment/acquire.py --target T3` for missing target.

### High CV (>5%)
```
WARNING: bwa-mem2 T3 has CV=8.2% (>5%)
```
**Cause**: High variance across replicates (hardware quirks, OS jitter, cache effects).
**Solution**: Review node logs, consider adding more replicates or filtering outliers.

## Verification

### Small-Scale Test
```bash
# Edit targets.conf → only T3 (MTB, 4.4 Mb)
./run_benchmark.sh --dry-run  # Array size: 2 aligners × 3 reps = 6
./run_benchmark.sh

# Monitor
squeue -u $USER
tail -f results/run_*/slurm-*.log

# Validate
ls results/run_*/T3/{bwa-mem2,minimap2}/rep_*/timing.json
cat results/run_*/results.json | jq .
```

### Full Benchmark
- Default targets: 9 targets (T1-T7b, exclude T8b/T8c)
- Array size: 9 × 2 × 3 = 54 tasks
- Expected runtime: ~4-6 hours (dominated by T1 chr1 249Mb, T2 chicken 1Gb)

### Success Criteria
- [ ] All 54 tasks complete without error
- [ ] CV < 5% for each (target, aligner) across 3 replicates
- [ ] `results.json` validates against score.py schema
- [ ] score.py produces BNT, CAS, CBS table

## Normalization Metrics

### BNT (Bandwidth-Normalized Throughput)
```
BNT = reads / (wall_time_s × threads × stream_triad_gbps)
```
Units: reads per thread-GB. Platform-independent speed measure.

### CAS (Composite Accuracy Score)
```
CAS = geometric_mean(PA, MCS, VA)
```
- **PA**: Placement accuracy (fraction within 10bp of true origin)
- **MCS**: MAPQ calibration score (1 − ECE)
- **VA**: Variant accuracy (F1 for SNP+indel, GIAB truth)

### CBS (Combined Benchmark Score)
```
CBS = CAS^α × BNT_norm^(1−α)  (default α=0.7)
```
Single-number leaderboard summary (accuracy-weighted).

### MEI (Memory Efficiency Index)
```
MEI = peak_rss_gb / genome_size_gb
```
Index "bloat factor" (lower is better).

## Known Limitations

### Phase 1 (Current)
- ✅ BWA-MEM2 and minimap2 working end-to-end
- ✅ PA computation automated (paftools.js mapeval)
- ✅ Read counting automated (FASTQ wc -l)
- ✅ Genome sizes from targets.conf (hard-coded)
- ⚠️ Phraya integration **blocked on issue #160** (batch-mode CLI)
- ⚠️ MCS (MAPQ calibration) still placeholder (requires MAPQ stratification)

### Phase 2 (Future)
- Add phraya wrapper once batch-mode ships (#160)
- Implement MCS computation (MAPQ-stratified PA analysis)
- Capture aligner versions from `module show` or `--version`
- Auto-populate genome sizes in targets.conf (agent to read manifests)

## References

- [Normalization spec](~/aligner_benchmark_normalization_spec.md): BNT/CAS/CBS framework
- [AGENTS.md](~/data-commons/test/benchmarking/alignment/AGENTS.md): Step-by-step protocol
- [Issue #160](https://github.com/CFSAN-Biostatistics/phraya/issues/160): Phraya batch-mode feature request
- [Plan file](~/.claude/plans/agile-gliding-stonebraker.md): Implementation design

## License

This benchmarking suite is part of the Phraya project. Same license applies.
