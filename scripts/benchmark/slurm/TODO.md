# Benchmark Suite TODO

## Before First Run

### 0. Build Phraya Binary

**Status:** REQUIRED — phraya wrapper expects release binary

**Task:** Build phraya release binary:
```bash
cd ~/phraya
cargo build --release
# Binary: target/release/phraya
```

**Verification:**
```bash
~/phraya/target/release/phraya --version
```

## Before First Run

### 1. Populate Genome Sizes in targets.conf

**Status:** REQUIRED — all entries currently have `0.0` placeholder

**Task:** Launch an agent to read manifest.json files and populate the 4th field (GENOME_SIZE_GB) in config/targets.conf.

**Agent prompt:**
```
Read manifest.json files from ~/data-commons/test/benchmarking/alignment/ for each target listed in config/targets.conf. Extract the canonical genome size from the `length_bp` field, convert to GB (divide by 1e9), and update targets.conf with the correct values.

Targets to process:
- T3: mycobacterium_tuberculosis/h37rv
- T4: staphylococcus_aureus/mrsa252
- T5: plasmodium_falciparum/3d7
- T6: clostridioides_difficile/630
- T7a: candida_albicans/sc5314_haploid
- T1: homo_sapiens/chr1
- T2: gallus_gallus/grcg7b
- T7b: candida_albicans/sc5314_diploid
- T8a: triticum_aestivum/chr3b
- T8b: triticum_aestivum/hexaploid (commented out)
- T8c: triticum_aestivum/aegilops_tauschii_aet_v4 (commented out)

Output format: TARGET_ID|ORGANISM_PATH|SIZE_CLASS|GENOME_SIZE_GB
```

**Expected output example:**
```
T3|mycobacterium_tuberculosis/h37rv|small|0.00441
T4|staphylococcus_aureus/mrsa252|small|0.00290
...
```

### 2. Verify HPC Module Names

**Status:** ASSUMED STABLE — confirm before first run

**Task:** SSH to HPC, run `module avail`, verify:
- `bwa-mem2` module exists
- `minimap2` module exists
- `samtools` module exists

If names differ, update `config/global.env`:
```bash
MODULE_BWA="actual-bwa-name"
MODULE_MINIMAP="actual-minimap-name"
MODULE_SAMTOOLS="actual-samtools-name"
```

### 3. Verify Data Staged

**Status:** ASSUMED PRESENT — verify at least one target

**Task:** Check that benchmark dataset is present:
```bash
ls ~/data-commons/test/benchmarking/alignment/mycobacterium_tuberculosis/h37rv/data/reference/reference.fasta
ls ~/data-commons/test/benchmarking/alignment/mycobacterium_tuberculosis/h37rv/data/reads/reads_1.fastq.gz
```

If missing, run:
```bash
cd ~/data-commons/test/benchmarking/alignment
python acquire.py --target T3  # Or --all for all targets
```

### 4. Test Node Discovery

**Status:** UNTESTED — may fail on non-SLURM systems

**Task:** Verify node discovery works:
```bash
cd ~/phraya/scripts/benchmark/slurm
./utils/nodelist.sh 3
```

**Expected output:** Comma-separated node list (e.g., `node01,node02,node03`)

**If fails:** Check SLURM partition name in `utils/nodelist.sh` line 8:
```bash
PARTITION="${SLURM_PARTITION:-batch}"  # Change "batch" to your partition
```

### 4a. Verify Phraya Thread Control

**Status:** IMPORTANT — phraya doesn't have -t flag yet

**Current approach:** phraya wrapper sets `RAYON_NUM_THREADS=$THREADS` env var

**Verification:**
```bash
cd ~/phraya
RAYON_NUM_THREADS=8 target/release/phraya --help
# Should run without error
```

**Known limitation:** If phraya spawns processes that don't respect RAYON_NUM_THREADS, fair thread comparison with BWA/minimap2 may be invalid.

## Before Full Run

### 5. Small-Scale Validation

**Status:** NOT YET RUN

**Task:** Test with single target (T3) to validate end-to-end:
```bash
cd ~/phraya/scripts/benchmark/slurm

# Edit config/targets.conf: comment out all except T3
# Or create a custom targets file:
cat > config/test_targets.conf <<EOF
T3|mycobacterium_tuberculosis/h37rv|small|0.00441
EOF

# Dry run
./run_benchmark.sh --targets config/test_targets.conf --dry-run
# Expected: 3 aligners × 3 reps = 9 tasks

# Real run
./run_benchmark.sh --targets config/test_targets.conf

# Monitor
squeue -u $USER
tail -f results/run_*/slurm-*.log

# Verify outputs
ls results/run_*/T3/{bwa-mem2,minimap2}/rep_*/timing.json
cat results/run_*/results.json | jq .
```

**Success criteria:**
- All 9 tasks complete without error
- `timing.json` files present for all replicates
- `results.json` contains valid PA and read counts (not 0.0 or "unknown")
- CV < 5% for both aligners

### 6. Verify PA Computation Dependencies

**Status:** UNTESTED

**Task:** Ensure paftools.js is available:
```bash
which paftools.js
# Should return path like /usr/local/bin/paftools.js or via minimap2 module

# If not found, install minimap2 or add to PATH
```

**Task:** Ensure samtools works:
```bash
module load samtools  # Or whatever your module name is
samtools --version
```

## After First Run

### 7. Calibrate Aggregation Timeout

**Status:** DEFAULT 30 MIN — may be too short for large targets

**Current:** `aggregate.slurm` has `#SBATCH --time=0:30:00`

**Check:** After first full run, check aggregation job log:
```bash
tail results/run_*/aggregate.log
```

**If timed out:** Increase time in `utils/aggregate.slurm`:
```bash
#SBATCH --time=1:00:00  # Or 2:00:00 for large targets
```

### 8. Validate score.py Compatibility

**Status:** NOT YET TESTED

**Task:** Feed results.json to score.py:
```bash
python ~/data-commons/test/benchmarking/alignment/score.py results/run_*/results.json --sensitivity
```

**Expected output:** BNT, CAS, CBS table

**If fails:** Check JSON schema mismatch, update aggregate_results.py field names to match score.py expectations.

## Optional Enhancements

### 9. Capture Aligner Versions

**Status:** PLACEHOLDER "unknown"

**Enhancement:** In `benchmark.slurm`, capture versions:
```bash
BWA_VERSION=$(bwa-mem2 version 2>&1 | head -1)
MINIMAP_VERSION=$(minimap2 --version 2>&1 | head -1)
# Write to $OUT_DIR/versions.txt
```

Aggregate in `aggregate_results.py`.

### 10. MCS (MAPQ Calibration Score) Computation

**Status:** DEFERRED TO PHASE 2

**Complexity:** Requires MAPQ-stratified PA analysis (PA per MAPQ bin: 0-9, 10-19, etc.)

**Block:** paftools.js mapeval doesn't output MAPQ bins by default. Need custom parser or switch to different eval tool.

---

## Deployment Checklist (Quick Reference)

**Before first run:**
- [ ] Populate genome sizes in targets.conf (agent task)
- [ ] Verify module names on HPC
- [ ] Verify benchmark data staged (at least T3)
- [ ] Test node discovery: `./utils/nodelist.sh 3`

**Small-scale validation:**
- [ ] Dry run with T3 only
- [ ] Real run with T3 only (9 tasks)
- [ ] Verify results.json has valid PA and read counts

**Full benchmark:**
- [ ] Run on 9 default targets (81 tasks)
- [ ] Check aggregation logs for timeouts
- [ ] Feed results.json to score.py
- [ ] Validate BNT/CAS/CBS output

**Post-run:**
- [ ] Document CV values (measurement noise)
- [ ] Archive results for reproducibility
- [ ] Update targets.conf with validated genome sizes
