# SLURM-Based HPC Aligner Benchmark

Production-grade SLURM benchmarking suite for comparing BWA-MEM2 and minimap2 across 16 benchmark targets with proper hardware normalization (STREAM Triad), 3-replicate averaging, and score.py-compatible output.

## Reporting rule

Benchmark write-ups (`benchmarking-N.md`, handoff results sections) are produced **only after
the entire suite has run on Reedling2** (the full targets × aligners × 3-replicate SLURM matrix
for the alphabet being reported). Do not substitute local or single-dataset spot checks, even
when they are quicker to run. If Reedling2 is unreachable, the data is missing, or the run
fails, stop and report the failure rather than falling back to a smaller local result.

Note: `config/targets.conf` is the short-read DNA set only (T1-T7b, T8a; T8b/T8c with `--large`), so
a plain `./run_benchmark.sh` is the 9-target matrix. Long-read, protein and assembly-pair targets
live in `targets_longread.conf`, `targets_protein.conf` and `targets_assembly.conf`. Set `DATA_ROOT`
to the staged data (`/nfs/scratch/Justin.Payne/data-commons/test/benchmarking/alignment`).

## Quick Start

```bash
cd ~/phraya/scripts/benchmark/slurm

# Build phraya first (required)
cd ~/phraya && cargo build --release

# Dry run (preview tasks without submitting)
cd ~/phraya/scripts/benchmark/slurm
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

### Phraya
- Build phraya: `cd ~/phraya && cargo build --release`
- Binary location: `~/phraya/target/release/phraya` (configurable via `PHRAYA_ROOT` in global.env)

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
│   ├── minimap2.sh           # minimap2 invocation + timing
│   └── phraya.sh             # phraya invocation + timing (batch-mode)
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
- Array size: `N_targets × 3_aligners × 3_replicates`
- Example: 9 targets × 3 aligners × 3 reps = 81 tasks

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
- phraya: sketches computed during `phraya plan` (no separate index)
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
./run_benchmark.sh --dry-run  # Array size: 3 aligners × 3 reps = 9
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
- Array size: 9 × 3 × 3 = 81 tasks
- Expected runtime: ~4-6 hours (dominated by T1 chr1 249Mb, T2 chicken 1Gb)

### Success Criteria
- [ ] All 81 tasks complete without error (9 targets × 3 aligners × 3 replicates)
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

## Protein-space benchmarking (ADR-0013)

`--alphabet protein` runs a separate aligner set (`diamond-blastp`, `blastp`,
`phraya-protein`, `phraya-protein-sensitive`) against `config/targets_protein.conf`
instead of `config/targets.conf`. Protein wrappers take **4** positional args
(`<proteome.fasta> <query.fasta> <out_dir> <threads>`, no paired reads) —
different from the DNA harness's 5-arg contract, since protein search has no
read pairing. See `../BENCHMARK_EXPANSION.md` for the full design rationale.

```bash
./run_benchmark.sh --alphabet protein
./run_benchmark.sh --alphabet both   # submits DNA and protein as separate array jobs
```

Protein targets are generated once by `../local/gen_synthetic_protein.py` and
staged like STREAM Triad's platform characterization — there's no external
"real proteome + known ground truth" corpus to point at the way `acquire.py`
stages real genomes for DNA targets. `results.json` entries for protein
targets carry `proteome_residues` (not `genome_size_gb`) and `rnt` (not
relevant to score.py's BNT, which is genome-GB-normalized) — Residue-Normalized
Throughput, `queries / (wall_s × threads × stream_triad_GBps)`.

**Known gap**: placement accuracy (`pa`) is only computed for `phraya-protein*`
today, via `phraya_accuracy_protein.py` (reads `gen_synthetic_protein.py`'s own
truth sidecar — no wgsim-style encoding exists for protein). DIAMOND/BLASTP
write tabular output, not SAM, and don't have a PA helper yet; their `pa` field
is `null`, not a fabricated `0.0`.

## Long-read benchmarking

`--alphabet long-read` runs `minimap2-ont` (`-x map-ont`), `minimap2-pb` (`-x map-pb`) and
`phraya-longread` (`phraya align --strategy long-read`) against `config/targets_longread.conf`
on single-end reads (4-arg wrapper contract `<ref> <reads.fq> <out_dir> <threads>`).

```bash
DATA_ROOT=/nfs/scratch/Justin.Payne/data-commons/test/benchmarking/alignment \
  ./run_benchmark.sh --alphabet long-read                 # 8 threads
THREADS=1 DATA_ROOT=... ./run_benchmark.sh --alphabet long-read   # single-thread pass
```

Per target the data must provide `data/reference/reference.fasta`, `data/reads/reads_1.fq` and
the generator's truth sidecar `data/reads/reads_1.fq.truth.tsv` (`gen_synthetic_long.py` format).

- `wall_seconds` is the aligner step only for every tool: `phraya align`, `minimap2` with a prebuilt
  index. `phraya plan` (read + reference sketching) is the analogue of minimap2's index build and is
  reported separately as `plan_time_s`. minimap2 indices are preset-specific
  (`reference.map-ont.mmi`, `reference.map-pb.mmi`) and built atomically.
- PA / precision / `n_correct_any` come from `utils/long_read_accuracy.py` for both phraya
  (`cross_space.phraya.queries`) and minimap2 (`alignment.sam`; best = primary record, retained
  placements = primary + secondary) under one rule: the best placement's start is within a tolerance
  of `true_start` (soft-clip-adjusted for SAM). Because noisy read ends are soft-clipped, a tight start
  test penalises clipping, so two tolerances are reported: `pa` at **2000 bp** (locus-correct; drives CAS)
  and `pa_start_100bp` (start precision). Planted-repeat targets (`data/reference/segdup_truth.bed`) add
  `n_repeat_inside`, `n_repeat_both_copies`, `n_overhang`, `n_overhang_correct`. Targets with no truth by
  design (real reads) get `agree_frac` etc.: agreement of each aligner's best placement with
  `minimap2-ont`'s primary, never PA.
- Run one pass per thread count (`THREADS=8` default, `THREADS=1`), each its own run directory. Long-read
  tasks can exceed the script's 4 h limit at 1 thread: `export SBATCH_TIMELIMIT=12:00:00` before submitting.
- `results.json` long-read entries add `bases_per_s`, `plan_time_s`, `wall_cv_pct`, `precision`,
  `n_truth`, `n_mapped`, `n_correct`, `n_correct_any`, `pa_tolerance_bp`, `pa_start_100bp`, and `sweep_of` (read from the
  target's data-commons `manifest.json`). `data-commons/.../score.py` auto-detects long-read results and
  scores them separately: BNT is bases-based (bases/s per thread per GB/s), CAS is the geometric mean of
  PA (MCS is unimplemented; a null PA, e.g. real reads, is left out of CAS), and `sweep_of` targets go to
  an amortization table instead of the aggregates. It refuses files that mix long-read with other alphabets.
- If any task fails, the dependent aggregation job never starts; run
  `utils/aggregate_results.py <run_dir> > <run_dir>/results.json` by hand (with `DATA_ROOT` set).

## Gap-affine scoring benchmarking (ADR-0014)

`phraya-sensitive-linear` and `phraya-sensitive-affine` (once ADR-0014 ships;
today only `phraya-sensitive-linear` explicitly forces `--gap-model linear`,
since bare `phraya-sensitive` becomes affine-by-default automatically once the
flag exists) run on the **same existing DNA targets**, unchanged —
dwgsim's default read simulation already includes realistic indels, so no new
target axis was needed.

**Indel Event Concordance (IEC)** — whether a read's true indel events are
reported as the same number of CIGAR I/D ops — is fully implemented and
immediately usable at local-bench scale (`../local/compute_indel_recovery.py`,
wired into `../local/run_local_bench.sh` via `INDEL_RATE`). It is **not**
computed at SLURM scale in this pass: that needs dwgsim's true-indel encoding
parsed out of its read names, a different (and unverified in this codebase)
format from the `wgsim` encoding the existing PA helpers already parse —
flagged as follow-up work, not implemented speculatively. `results.json`
entries do carry a `gap_model` field (`linear`/`affine`, parsed from the
aligner-variant name) so linear-vs-affine timing/RSS/PA comparisons work today
even without SLURM-scale IEC.

## Known Limitations

### Phase 1 (Complete)
- ✅ BWA-MEM2, minimap2, and **phraya** working end-to-end
- ✅ Phraya batch-mode integration (issue #160 resolved, commit 7551134)
- ✅ PA computation automated (paftools.js mapeval for BWA/minimap2)
- ✅ Read counting automated (FASTQ wc -l)
- ✅ Genome sizes from targets.conf (hard-coded)

### Phase 2 (Future)
- ⚠️ Phraya PA computation: requires SAM export or .phraya → SAM conversion
- ⚠️ MCS (MAPQ calibration) still placeholder (requires MAPQ stratification)
- Capture aligner versions from `module show` or `--version`
- Auto-populate genome sizes in targets.conf (agent to read manifests)

## References

- [Normalization spec](~/aligner_benchmark_normalization_spec.md): BNT/CAS/CBS framework
- [AGENTS.md](~/data-commons/test/benchmarking/alignment/AGENTS.md): Step-by-step protocol
- [Issue #160](https://github.com/CFSAN-Biostatistics/phraya/issues/160): Phraya batch-mode feature request
- [Plan file](~/.claude/plans/agile-gliding-stonebraker.md): Implementation design

## License

This benchmarking suite is part of the Phraya project. Same license applies.
