#!/usr/bin/env python3
"""Generate compressed versions of test fixtures."""
import gzip
import bz2
from pathlib import Path

fixtures_dir = Path(__file__).parent

# Compress sample.fasta
with open(fixtures_dir / "sample.fasta", "rb") as f:
    content = f.read()
    with gzip.open(fixtures_dir / "sample.fasta.gz", "wb") as gz:
        gz.write(content)
    with bz2.open(fixtures_dir / "sample.fasta.bz2", "wb") as bz:
        bz.write(content)

# Compress sample.fastq
with open(fixtures_dir / "sample.fastq", "rb") as f:
    content = f.read()
    with gzip.open(fixtures_dir / "sample.fastq.gz", "wb") as gz:
        gz.write(content)
    with bz2.open(fixtures_dir / "sample.fastq.bz2", "wb") as bz:
        bz.write(content)

print("Created compressed fixtures:")
print("  - sample.fasta.gz")
print("  - sample.fasta.bz2")
print("  - sample.fastq.gz")
print("  - sample.fastq.bz2")
