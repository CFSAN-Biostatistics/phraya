#!/bin/bash
# Acquire reference sequences and simulate reads for medium and large benchmark targets.
# Targets: T1 (human chr1), T2 (chicken GRCg7b), T7b (diploid Candida),
#          T8a (wheat chr3b), T8b (hexaploid wheat, --large), T8c (Ae. tauschii, --large)
#
# Usage: acquire_medium_large.sh [--large] [T1 T2 T7b T8a ...]
#   --large  also acquire T8b (17Gb hexaploid wheat) and T8c (4.3Gb Ae. tauschii)

set -euo pipefail

DATA_ROOT="${DATA_ROOT:-$HOME/data-commons/test/benchmarking/alignment}"
TOOLS_ENV="$HOME/phraya/envs/bench-tools"
WGSIM="/nfs/software/apps/micromamba/1.5.8/envs/samtools-v1.20/bin/wgsim"
SAMTOOLS="/nfs/software/apps/micromamba/1.5.8/envs/samtools-v1.20/bin/samtools"
BWA="/nfs/software/apps/micromamba/1.5.8/envs/bwa-v0.7.18/bin/bwa"
MINIMAP2="/nfs/software/apps/micromamba/1.5.8/envs/minimap2-v2.28/bin/minimap2"

DWGSIM="$TOOLS_ENV/bin/dwgsim"
ART="$TOOLS_ENV/bin/art_illumina"
EFETCH="${TOOLS_ENV}/bin/efetch"
# efetch may already exist in system envs
[[ -f "$EFETCH" ]] || EFETCH="/nfs/software/apps/micromamba/1.5.8/envs/entrez-direct-22.7/bin/efetch"
[[ -f "$EFETCH" ]] || EFETCH="/nfs/software/apps/micromamba/1.5.8/envs/seqsero2-v1.2.1/bin/efetch"
DATASETS="$TOOLS_ENV/bin/datasets"

INCLUDE_LARGE=0
TARGETS=()

for arg in "$@"; do
    case "$arg" in
        --large) INCLUDE_LARGE=1 ;;
        T*) TARGETS+=("$arg") ;;
    esac
done

if [[ ${#TARGETS[@]} -eq 0 ]]; then
    TARGETS=(T1 T2 T7b T8a)
    [[ $INCLUDE_LARGE -eq 1 ]] && TARGETS+=(T8b T8c)
fi

# ---- helpers ----

check_tools() {
    for bin in "$DWGSIM" "$ART" "$EFETCH" "$DATASETS"; do
        [[ -f "$bin" ]] || { echo "ERROR: $bin not found — run micromamba install first" >&2; exit 1; }
    done
}

build_indexes() {
    local REF=$1
    [[ -f "${REF}.bwt" ]] || { echo "  Building BWA index..."; $BWA index "$REF"; }
    [[ -f "${REF%.fasta}.mmi" ]] || {
        echo "  Building minimap2 index..."
        $MINIMAP2 -d "${REF%.fasta}.mmi" "$REF"
    }
    [[ -f "${REF}.fai" ]] || $SAMTOOLS faidx "$REF"
}

subsample_30k() {
    local READS_DIR=$1
    local R1="$READS_DIR/reads_1.fastq.gz"
    local R2="$READS_DIR/reads_2.fastq.gz"
    local R1_30K="$READS_DIR/reads_1_30k.fastq.gz"
    local R2_30K="$READS_DIR/reads_2_30k.fastq.gz"
    if [[ ! -f "$R1_30K" ]] || [[ ! -f "$R2_30K" ]]; then
        echo "  Subsampling to 30k pairs..."
        zcat "$R1" | head -120000 | gzip > "$R1_30K"
        zcat "$R2" | head -120000 | gzip > "$R2_30K"
    fi
}

sim_dwgsim() {
    # dwgsim -1 150 -2 150 -C COVERAGE REF reads
    local REF=$1 COVERAGE=$2 READS_DIR=$3 SEED=$4
    mkdir -p "$READS_DIR/tmp_sim"
    $DWGSIM -1 150 -2 150 -C "$COVERAGE" -z "$SEED" "$REF" "$READS_DIR/tmp_sim/reads"
    mv "$READS_DIR/tmp_sim/reads.bwa.read1.fastq.gz" "$READS_DIR/reads_1.fastq.gz"
    mv "$READS_DIR/tmp_sim/reads.bwa.read2.fastq.gz" "$READS_DIR/reads_2.fastq.gz"
    rm -rf "$READS_DIR/tmp_sim"
}

sim_art() {
    # art_illumina -ss HS25 -p -l 150 -f COVERAGE -i REF -o prefix
    local REF=$1 COVERAGE=$2 READS_DIR=$3 SEED=$4
    mkdir -p "$READS_DIR/tmp_sim"
    $ART -ss HS25 -p -l 150 -f "$COVERAGE" -i "$REF" \
        -rs "$SEED" -o "$READS_DIR/tmp_sim/reads" -q
    gzip -c "$READS_DIR/tmp_sim/reads1.fq" > "$READS_DIR/reads_1.fastq.gz"
    gzip -c "$READS_DIR/tmp_sim/reads2.fq" > "$READS_DIR/reads_2.fastq.gz"
    rm -rf "$READS_DIR/tmp_sim"
}

# ---- per-target acquisition ----

acquire_T1() {
    local ORG=homo_sapiens/chr1
    local REF_DIR="$DATA_ROOT/$ORG/data/reference"
    local READS_DIR="$DATA_ROOT/$ORG/data/reads"
    local REF="$REF_DIR/reference.fasta"
    mkdir -p "$REF_DIR" "$READS_DIR"

    if [[ ! -f "$REF" ]]; then
        echo "[T1] Downloading human chr1 (NC_000001.11) via efetch (~750 MB)..."
        $EFETCH -db nucleotide -id NC_000001.11 -format fasta > "$REF"
    fi
    build_indexes "$REF"

    if [[ ! -f "$READS_DIR/reads_1.fastq.gz" ]]; then
        echo "[T1] Simulating reads (dwgsim, 30x, 150bp PE)..."
        sim_dwgsim "$REF" 30 "$READS_DIR" 42
    fi
    subsample_30k "$READS_DIR"
    echo "[T1] Done."
}

acquire_T2() {
    local ORG=gallus_gallus/grcg7b
    local REF_DIR="$DATA_ROOT/$ORG/data/reference"
    local READS_DIR="$DATA_ROOT/$ORG/data/reads"
    local REF="$REF_DIR/reference.fasta"
    mkdir -p "$REF_DIR" "$READS_DIR"

    if [[ ! -f "$REF" ]]; then
        echo "[T2] Downloading chicken GRCg7b (GCF_016699485.2) via ncbi-datasets (~3.2 GB)..."
        TMP="$REF_DIR/tmp_dl"
        mkdir -p "$TMP"
        $DATASETS download genome accession GCF_016699485.2 \
            --include genome --filename "$TMP/ncbi_dataset.zip"
        unzip -o "$TMP/ncbi_dataset.zip" -d "$TMP"
        cat "$TMP/ncbi_dataset/data/GCF_016699485.2/"*.fna > "$REF"
        rm -rf "$TMP"
    fi
    build_indexes "$REF"

    if [[ ! -f "$READS_DIR/reads_1.fastq.gz" ]]; then
        echo "[T2] Simulating reads (dwgsim, 30x, 150bp PE)..."
        sim_dwgsim "$REF" 30 "$READS_DIR" 42
    fi
    subsample_30k "$READS_DIR"
    echo "[T2] Done."
}

acquire_T7b() {
    local ORG=candida_albicans/sc5314_diploid
    local REF_DIR="$DATA_ROOT/$ORG/data/reference"
    local READS_DIR="$DATA_ROOT/$ORG/data/reads"
    local REF="$REF_DIR/reference.fasta"
    mkdir -p "$REF_DIR" "$READS_DIR"

    if [[ ! -f "$REF" ]]; then
        echo "[T7b] Downloading diploid C. albicans SC5314 Assembly 22 from CGD..."
        # Direct download of Assembly 22 chromosomes from CGD
        BASE="http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/chromosomes"
        for CHR in Ca22chr1A_C_albicans_SC5314 Ca22chr1B_C_albicans_SC5314 \
                   Ca22chr2A_C_albicans_SC5314 Ca22chr2B_C_albicans_SC5314 \
                   Ca22chr3A_C_albicans_SC5314 Ca22chr3B_C_albicans_SC5314 \
                   Ca22chr4A_C_albicans_SC5314 Ca22chr4B_C_albicans_SC5314 \
                   Ca22chr5A_C_albicans_SC5314 Ca22chr5B_C_albicans_SC5314 \
                   Ca22chr6A_C_albicans_SC5314 Ca22chr6B_C_albicans_SC5314 \
                   Ca22chr7A_C_albicans_SC5314 Ca22chrRA_C_albicans_SC5314; do
            curl -sf "$BASE/$CHR.fasta.gz" | gzip -d >> "$REF" 2>/dev/null || true
        done
        # If CGD fetch failed, try the assembly zip
        if [[ ! -s "$REF" ]]; then
            echo "[T7b] CGD chromosome fetch failed, trying assembly zip..."
            TMP="$REF_DIR/tmp_dl"; mkdir -p "$TMP"
            curl -L "http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_version_A22-s07-m01-r110_chromosomes.fasta.gz" \
                | gzip -d > "$REF" 2>/dev/null || true
            rm -rf "$TMP"
        fi
        [[ -s "$REF" ]] || { echo "[T7b] ERROR: failed to download reference" >&2; return 1; }
    fi
    build_indexes "$REF"

    if [[ ! -f "$READS_DIR/reads_1.fastq.gz" ]]; then
        echo "[T7b] Simulating reads (dwgsim, 50x, 150bp PE)..."
        sim_dwgsim "$REF" 50 "$READS_DIR" 42
    fi
    subsample_30k "$READS_DIR"
    echo "[T7b] Done."
}

acquire_T8a() {
    local ORG=triticum_aestivum/chr3b
    local REF_DIR="$DATA_ROOT/$ORG/data/reference"
    local READS_DIR="$DATA_ROOT/$ORG/data/reads"
    local REF="$REF_DIR/reference.fasta"
    mkdir -p "$REF_DIR" "$READS_DIR"

    if [[ ! -f "$REF" ]]; then
        echo "[T8a] Downloading wheat chr3B from URGI (~831 MB)..."
        # IWGSC RefSeq v1.0 chr3B from ENA/URGI
        # Primary URL: ftp ENA
        curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/519/105/GCA_900519105.1_IWGSC_v1.0/GCA_900519105.1_IWGSC_v1.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr3B.fna.gz" \
            | gzip -d > "$REF" 2>/dev/null || {
            echo "[T8a] NCBI fetch failed, trying EBI..."
            curl -L "https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/ca/CABB*.fasta.gz" \
                | gzip -d > "$REF" 2>/dev/null || true
        }
        [[ -s "$REF" ]] || { echo "[T8a] ERROR: failed to download chr3B reference" >&2; return 1; }
    fi
    build_indexes "$REF"

    if [[ ! -f "$READS_DIR/reads_1.fastq.gz" ]]; then
        echo "[T8a] Simulating reads (art_illumina HS25, 10x, 150bp PE)..."
        sim_art "$REF" 10 "$READS_DIR" 42
    fi
    subsample_30k "$READS_DIR"
    echo "[T8a] Done."
}

acquire_T8c() {
    local ORG=triticum_aestivum/aegilops_tauschii_aet_v4
    local REF_DIR="$DATA_ROOT/$ORG/data/reference"
    local READS_DIR="$DATA_ROOT/$ORG/data/reads"
    local REF="$REF_DIR/reference.fasta"
    mkdir -p "$REF_DIR" "$READS_DIR"

    if [[ ! -f "$REF" ]]; then
        echo "[T8c] Downloading Ae. tauschii AET v4 (GCA_002575655.1, ~4.3 GB)..."
        TMP="$REF_DIR/tmp_dl"; mkdir -p "$TMP"
        $DATASETS download genome accession GCA_002575655.1 \
            --include genome --filename "$TMP/ncbi_dataset.zip"
        unzip -o "$TMP/ncbi_dataset.zip" -d "$TMP"
        cat "$TMP/ncbi_dataset/data/GCA_002575655.1/"*.fna > "$REF"
        rm -rf "$TMP"
    fi
    build_indexes "$REF"

    if [[ ! -f "$READS_DIR/reads_1.fastq.gz" ]]; then
        echo "[T8c] Simulating reads (art_illumina HS25, 10x, 150bp PE)..."
        sim_art "$REF" 10 "$READS_DIR" 42
    fi
    subsample_30k "$READS_DIR"
    echo "[T8c] Done."
}

acquire_T8b() {
    local ORG=triticum_aestivum/hexaploid
    local REF_DIR="$DATA_ROOT/$ORG/data/reference"
    local READS_DIR="$DATA_ROOT/$ORG/data/reads"
    local REF="$REF_DIR/reference.fasta"
    mkdir -p "$REF_DIR" "$READS_DIR"

    if [[ ! -f "$REF" ]] || [[ ! -s "$REF" ]]; then
        # Remove any truncated file from a previous failed attempt
        rm -f "$REF"
        echo "[T8b] Downloading IWGSC RefSeq v1.0 hexaploid wheat (sequential per chromosome to avoid append races)..."
        BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/519/105/GCA_900519105.1_IWGSC_v1.0/GCA_900519105.1_IWGSC_v1.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA"
        for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D \
                   chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D \
                   chr7A chr7B chr7D; do
            echo "[T8b]   downloading $CHR..."
            curl -sL "$BASE/${CHR}.fna.gz" | gzip -d >> "$REF" || { echo "[T8b] WARNING: $CHR failed, continuing"; }
        done
        # chrUn (unanchored scaffolds)
        curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/519/105/GCA_900519105.1_IWGSC_v1.0/GCA_900519105.1_IWGSC_v1.0_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrUn.fna.gz" \
            | gzip -d >> "$REF" 2>/dev/null || true
        [[ -s "$REF" ]] || { echo "[T8b] ERROR: reference empty after download" >&2; return 1; }
        echo "[T8b] Reference assembled: $(du -sh $REF | cut -f1)"
    fi
    build_indexes "$REF"

    if [[ ! -f "$READS_DIR/reads_1.fastq.gz" ]]; then
        echo "[T8b] Simulating reads (art_illumina HS25, 5x, 150bp PE)..."
        sim_art "$REF" 5 "$READS_DIR" 42
    fi
    subsample_30k "$READS_DIR"
    echo "[T8b] Done."
}

# ---- main ----

check_tools

echo "Acquiring targets: ${TARGETS[*]}"

for TID in "${TARGETS[@]}"; do
    case "$TID" in
        T1)  acquire_T1 &  ;;
        T2)  acquire_T2    ;;  # serial: large download, runs alone
        T7b) acquire_T7b & ;;
        T8a) acquire_T8a & ;;
        T8b) acquire_T8b   ;;  # serial: 17 GB
        T8c) acquire_T8c & ;;
        *) echo "Unknown target: $TID" >&2 ;;
    esac
done
wait
echo "=== All targets acquired ==="
