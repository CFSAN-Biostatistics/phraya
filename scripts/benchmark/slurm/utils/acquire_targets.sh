#!/bin/bash
# Acquire reference sequences and simulate reads for benchmark targets.
# Uses curl (NCBI E-utilities) for single-sequence genomes,
# micromamba+ncbi-datasets-cli for multi-sequence assemblies.
# Simulates reads with wgsim (from samtools env).
#
# Usage: acquire_targets.sh [T3 T4 T5 T6 T7a ...]
#   (defaults to all small/medium small targets except T3 which already has data)

set -euo pipefail

DATA_ROOT="${DATA_ROOT:-$HOME/data-commons/test/benchmarking/alignment}"
WGSIM="/nfs/software/apps/micromamba/1.5.8/envs/samtools-v1.20/bin/wgsim"
MAMBA="/nfs/software/apps/micromamba/1.5.8/bin/micromamba"
DATASETS_ENV="$HOME/phraya/envs/ncbi-datasets"

TARGETS=("${@:-T4 T5 T6 T7a}")

# ---- Tool setup ----

# Install ncbi-datasets-cli if not present
if [[ ! -f "$DATASETS_ENV/bin/datasets" ]]; then
    echo "Installing ncbi-datasets-cli..."
    "$MAMBA" create -p "$DATASETS_ENV" ncbi-datasets-cli -c conda-forge -c bioconda -y 2>&1 | tail -5
fi
DATASETS="$DATASETS_ENV/bin/datasets"

# ---- Per-target acquisition ----

acquire_target() {
    local TID=$1
    case "$TID" in
        T4) ORG=staphylococcus_aureus/mrsa252;  REFSEQ=NC_002952.2;       METHOD=efetch ;;
        T5) ORG=plasmodium_falciparum/3d7;      ASSEMBLY=GCF_000002765.4; METHOD=datasets; COVERAGE=50 ;;
        T6) ORG=clostridioides_difficile/630;   REFSEQ=NC_009089.1;       METHOD=efetch ;;
        T7a) ORG=candida_albicans/sc5314_haploid; ASSEMBLY=GCF_000182965.3; METHOD=datasets; COVERAGE=50 ;;
        *) echo "Unknown target: $TID" >&2; return 1 ;;
    esac

    COVERAGE="${COVERAGE:-100}"
    REF_DIR="$DATA_ROOT/$ORG/data/reference"
    READS_DIR="$DATA_ROOT/$ORG/data/reads"
    mkdir -p "$REF_DIR" "$READS_DIR"
    REF="$REF_DIR/reference.fasta"

    # Download reference
    if [[ ! -f "$REF" ]]; then
        echo "[$TID] Downloading reference ($METHOD)..."
        case "$METHOD" in
            efetch)
                curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${REFSEQ}&rettype=fasta&retmode=text" \
                    > "$REF"
                ;;
            datasets)
                TMP="$REF_DIR/tmp_dl"
                mkdir -p "$TMP"
                "$DATASETS" download genome accession "$ASSEMBLY" \
                    --include genome --filename "$TMP/ncbi_dataset.zip"
                unzip -o "$TMP/ncbi_dataset.zip" -d "$TMP"
                cat "$TMP/ncbi_dataset/data/$ASSEMBLY/"*.fna > "$REF"
                rm -rf "$TMP"
                ;;
        esac
        echo "[$TID] Reference downloaded: $(wc -l < "$REF") lines"
    else
        echo "[$TID] Reference already exists, skipping download"
    fi

    # Build BWA index (needed for bwa-pipeline wrapper flock check)
    if [[ ! -f "${REF}.bwt" ]]; then
        echo "[$TID] Building BWA index..."
        /nfs/software/apps/micromamba/1.5.8/envs/bwa-v0.7.18/bin/bwa index "$REF"
    fi

    # Build minimap2 index
    if [[ ! -f "${REF%.fasta}.mmi" ]]; then
        echo "[$TID] Building minimap2 index..."
        /nfs/software/apps/micromamba/1.5.8/envs/minimap2-v2.28/bin/minimap2 \
            -d "${REF%.fasta}.mmi" "$REF"
    fi

    # Simulate reads with wgsim if not present
    READS_1="$READS_DIR/reads_1.fastq.gz"
    READS_1_30K="$READS_DIR/reads_1_30k.fastq.gz"

    if [[ ! -f "$READS_1_30K" ]]; then
        if [[ ! -f "$READS_1" ]]; then
            echo "[$TID] Simulating reads (${COVERAGE}x coverage, wgsim)..."
            GENOME_SIZE=$(awk '/^>/{next}{s+=length($0)}END{print s}' "$REF")
            N_READS=$(( GENOME_SIZE * COVERAGE / 2 / 150 ))
            "$WGSIM" -N "$N_READS" -1 150 -2 150 -e 0.005 -r 0.001 \
                "$REF" \
                "$READS_DIR/reads_1_full.fastq" "$READS_DIR/reads_2_full.fastq"
            gzip -c "$READS_DIR/reads_1_full.fastq" > "$READS_1"
            gzip -c "$READS_DIR/reads_2_full.fastq" > "$READS_DIR/reads_2.fastq.gz"
            rm "$READS_DIR/reads_1_full.fastq" "$READS_DIR/reads_2_full.fastq"
            echo "[$TID] Full reads: $N_READS pairs"
        fi

        # Subsample to 30k pairs
        echo "[$TID] Subsampling to 30k pairs..."
        zcat "$READS_1" | head -120000 | gzip > "$READS_1_30K"
        zcat "$READS_DIR/reads_2.fastq.gz" | head -120000 | gzip > "$READS_DIR/reads_2_30k.fastq.gz"
        echo "[$TID] 30k subsample done"
    else
        echo "[$TID] 30k reads already exist, skipping"
    fi

    echo "[$TID] Done."
}

for TID in "${TARGETS[@]}"; do
    acquire_target "$TID" &
done
wait
echo "All targets acquired."
