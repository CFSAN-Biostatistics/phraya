#!/bin/bash
# Generate compressed test fixtures
cd "$(dirname "$0")"

# Compress FASTA files
gzip -c sample.fasta > sample.fasta.gz
bzip2 -c sample.fasta > sample.fasta.bz2

# Compress FASTQ files
gzip -c sample.fastq > sample.fastq.gz
bzip2 -c sample.fastq > sample.fastq.bz2

echo "Created compressed fixtures"
