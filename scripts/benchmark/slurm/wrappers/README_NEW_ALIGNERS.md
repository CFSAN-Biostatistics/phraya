# New Aligner Wrappers

Added three new aligners to the benchmark harness (issues #211, #212, #213):

## Bowtie2

**Wrapper:** `bowtie2.sh`  
**Tool detection:** Checks `$BOWTIE2_BIN`, falls back to PATH  
**Preset:** `--very-sensitive` (max recall for variant calling)  
**Index:** `.bt2` files (6 total), flock-protected build  
**Output:** SAM

**Installation options:**
```bash
# Conda
conda install -c bioconda bowtie2

# Apptainer/Singularity (HPC without root)
apptainer pull docker://biocontainers/bowtie2:v2.5.1_cv1
# Then set BOWTIE2_BIN to wrapper script calling apptainer exec
```

## minibwa

**Wrapper:** `minibwa.sh`  
**Tool detection:** Checks `$MINIBWA_BIN`, falls back to PATH  
**Algorithm:** Hybrid bwa-mem/minimap2 (3× faster than bwa-mem)  
**Index:** `.l2b` + `.mbw` files, flock-protected build  
**Output:** SAM (not byte-identical to bwa-mem due to minimap2 chaining)

**Installation options:**
```bash
# Build from source
git clone https://github.com/lh3/minibwa
cd minibwa && make
# Dependencies: zlib, SSE4.2/NEON

# Set path
export MINIBWA_BIN=/path/to/minibwa/minibwa
```

## rammap

**Wrapper:** `rammap.sh`  
**Tool detection:** Checks `$RAMMAP_BIN`, falls back to PATH  
**Algorithm:** Pure Rust minimap2 clone with SIMD DP  
**Preset:** `-x sr` (short reads)  
**Index:** None (inline with alignment, like minimap2)  
**Output:** SAM (identical to minimap2)

**Performance note:** ~20-40% slower than minimap2 on short reads, but eliminates C build dependencies.

**Installation options:**
```bash
# Cargo
cargo install --git https://github.com/jwanglab/rammap

# Binary release
wget https://github.com/jwanglab/rammap/releases/download/v1.1.1/rammap-x86_64-linux
chmod +x rammap-x86_64-linux
export RAMMAP_BIN=$PWD/rammap-x86_64-linux

# Apptainer (if available as container)
# Set RAMMAP_BIN to wrapper calling apptainer exec
```

## Integration

All three wrappers follow the existing pattern:
- Accept 5 args: `<ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>`
- Output `timing.txt` with standard fields (wall_seconds, peak_rss_gb, n_aligned, unaligned_frac)
- SAM output → `alignment.sam`
- Flock-protected index builds (bowtie2, minibwa; rammap has no index step)

Benchmark runner and SLURM scripts updated to include all three (9 aligners total now).
