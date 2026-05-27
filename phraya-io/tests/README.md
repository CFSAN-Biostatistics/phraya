# FASTA/FASTQ I/O Tests - Issue #2

This directory contains comprehensive acceptance tests for FASTA/FASTQ I/O functionality with auto-detection.

## Test Status: RED (All Tests Failing)

All tests are currently FAILING because the implementation does not exist yet. This is intentional - these are TDD acceptance tests written in the RED phase.

## Test Files

### `fasta_fastq_io.rs`
Integration tests covering all acceptance criteria:
- Parse FASTA (plain, .gz, .bz2) with no quality field
- Parse FASTQ (plain, .gz, .bz2) with quality field populated
- Auto-detection via magic bytes and extension fallback
- All 6 combinations (FASTA/FASTQ × plain/gz/bz2)
- Error handling for malformed files
- Streaming interface (iterator, not Vec)

**Test count**: 22 tests covering:
- 6 basic parsing tests (one per format/compression combo)
- 3 auto-detection tests (magic bytes + extension fallback)
- 2 streaming interface tests
- 4 error handling tests
- 4 edge case tests (empty files, wrapped sequences, quality validation)
- 3 iterator behavior tests

### `api_surface.rs`
API contract tests documenting the expected public interface:
- `parse_sequences()` function signature
- `SequenceReader` iterator behavior
- `SequenceFormat::detect()` and query methods
- `ParseError` error variants
- `Sequence` type structure from phraya-core

**Test count**: 10 tests

## Test Fixtures

See `fixtures/README.md` for details on test data files.

Before running tests, generate compressed fixtures:
```bash
cd phraya-io/tests/fixtures
bash generate_fixtures.sh
```

## Expected Behavior

### Current (RED phase)
```
error[E0432]: unresolved import `phraya_core::Sequence`
error[E0432]: unresolved imports `phraya_io::parse_sequences`, `phraya_io::ParseError`, `phraya_io::SequenceFormat`
```

All tests fail at compile time because types don't exist.

### After Implementation (GREEN phase)
All 32 tests should pass, demonstrating:
- Correct parsing of all 6 format/compression combinations
- Proper auto-detection logic
- Streaming interface with lazy evaluation
- Comprehensive error handling
- Edge case support

## Implementation Requirements

To make these tests pass, implement:

1. **In phraya-core**:
   - `Sequence` struct with fields: `id`, `description`, `data`, `quality`, `pairing_info`
   - Traits: `Clone`, `Debug`, `Serialize`, `Deserialize`

2. **In phraya-io**:
   - `ParseError` enum with variants:
     - `IoError { path, source }`
     - `MalformedEntry { line, reason }`
     - `QualityLengthMismatch { seq_len, qual_len }`
     - `UnsupportedFormat { path, reason }`
   - `SequenceFormat` enum/struct with detection and query methods
   - `SequenceReader` iterator type
   - `parse_sequences()` function returning `Result<SequenceReader, ParseError>`
   - Compression support via flate2 (gzip) and bzip2
   - Auto-detection logic based on magic bytes with extension fallback

## Test Coverage

| Acceptance Criterion | Test Coverage |
|---------------------|---------------|
| Parse FASTA (plain) | `test_parse_fasta_plain` |
| Parse FASTA (.gz) | `test_parse_fasta_gzipped` |
| Parse FASTA (.bz2) | `test_parse_fasta_bz2` |
| Parse FASTQ (plain) | `test_parse_fastq_plain` |
| Parse FASTQ (.gz) | `test_parse_fastq_gzipped` |
| Parse FASTQ (.bz2) | `test_parse_fastq_bz2` |
| Auto-detection (magic bytes) | `test_auto_detection_magic_bytes_gzip`, `test_auto_detection_magic_bytes_bzip2` |
| Auto-detection (extension) | `test_auto_detection_extension_fallback` |
| Streaming interface | `test_streaming_interface`, `test_streaming_large_file_memory_efficiency`, `test_iterator_is_lazy` |
| Error handling | `test_error_malformed_fasta`, `test_error_malformed_fastq`, `test_error_nonexistent_file`, `test_error_unsupported_format` |
| Edge cases | `test_empty_file`, `test_sequence_with_newlines_in_data`, `test_fastq_quality_score_validation`, `test_fasta_multiline_description` |
| Complete matrix | `test_compression_format_combination_matrix` |
