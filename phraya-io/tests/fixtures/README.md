# Test Fixtures for FASTA/FASTQ I/O

This directory contains test fixtures for comprehensive testing of FASTA/FASTQ parsing with compression support.

## Generating Compressed Fixtures

Before running tests, generate the compressed versions of the fixture files:

```bash
cd phraya-io/tests/fixtures
bash generate_fixtures.sh
```

This will create:
- `sample.fasta.gz` - Gzip-compressed FASTA
- `sample.fasta.bz2` - Bzip2-compressed FASTA
- `sample.fastq.gz` - Gzip-compressed FASTQ
- `sample.fastq.bz2` - Bzip2-compressed FASTQ

## Fixture Files

### Valid Files
- `sample.fasta` - 3 sequences in FASTA format
- `sample.fastq` - 3 sequences in FASTQ format with quality scores

### Malformed Files (for error testing)
- `malformed.fasta` - FASTA with missing headers and empty sequences
- `malformed.fastq` - FASTQ with missing quality lines and length mismatches

## Test Coverage

These fixtures support testing all 6 combinations required by acceptance criteria:
1. FASTA plain
2. FASTA gzipped
3. FASTA bzip2
4. FASTQ plain
5. FASTQ gzipped
6. FASTQ bzip2
