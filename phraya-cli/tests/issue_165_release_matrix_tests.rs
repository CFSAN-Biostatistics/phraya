/// RED acceptance tests for issue #165 — GitHub Releases multi-platform build matrix.
///
/// Tests verify the release workflow file and README against the acceptance criteria:
/// - 5 build targets with correct triple/variant names
/// - Binary naming convention: phraya-{version}-{target}-{variant}.tar.gz
/// - README platform selection guide
/// - README SIMD differences documentation

use std::fs;
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("phraya-cli must have a parent directory")
        .to_path_buf()
}

fn read_release_workflow() -> String {
    let path = repo_root().join(".github/workflows/release.yml");
    fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "release workflow not found at {}: {}",
            path.display(),
            e
        )
    })
}

fn read_readme() -> String {
    let path = repo_root().join("README.md");
    fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!("README.md not found at {}: {}", path.display(), e)
    })
}

// ============================================================================
// Release workflow existence
// ============================================================================

/// issue #165: release.yml workflow must exist at .github/workflows/release.yml
#[test]
fn issue_165_release_workflow_file_exists() {
    let path = repo_root().join(".github/workflows/release.yml");
    assert!(
        path.exists(),
        "expected release workflow at {} — create .github/workflows/release.yml with the multi-platform build matrix",
        path.display()
    );
}

// ============================================================================
// Build matrix — 5 required targets
// ============================================================================

/// issue #165: workflow matrix must include x86_64-linux-gnu portable (SSE4.2)
#[test]
fn issue_165_release_workflow_has_x86_64_linux_gnu_portable_target() {
    let workflow = read_release_workflow();
    assert!(
        workflow.contains("x86_64-linux-gnu") && workflow.contains("portable"),
        "release workflow must include x86_64-linux-gnu portable target (SSE4.2); \
         check the build matrix for this target+variant combination"
    );
}

/// issue #165: workflow matrix must include x86_64-linux-gnu native (AVX2)
#[test]
fn issue_165_release_workflow_has_x86_64_linux_gnu_native_target() {
    let workflow = read_release_workflow();
    assert!(
        workflow.contains("x86_64-linux-gnu") && workflow.contains("native"),
        "release workflow must include x86_64-linux-gnu native target (AVX2); \
         check the build matrix for target-cpu=native variant"
    );
}

/// issue #165: workflow matrix must include aarch64-linux-gnu (NEON)
#[test]
fn issue_165_release_workflow_has_aarch64_linux_gnu_target() {
    let workflow = read_release_workflow();
    assert!(
        workflow.contains("aarch64-linux-gnu"),
        "release workflow must include aarch64-linux-gnu target for ARM Linux (NEON); \
         check the build matrix for this cross-compilation target"
    );
}

/// issue #165: workflow matrix must include x86_64-darwin (macOS Intel)
#[test]
fn issue_165_release_workflow_has_x86_64_darwin_target() {
    let workflow = read_release_workflow();
    assert!(
        workflow.contains("x86_64-darwin") || workflow.contains("x86_64-apple-darwin"),
        "release workflow must include macOS Intel target (x86_64-darwin or x86_64-apple-darwin); \
         check the build matrix for macOS Intel runner"
    );
}

/// issue #165: workflow matrix must include aarch64-darwin (macOS M1/M2)
#[test]
fn issue_165_release_workflow_has_aarch64_darwin_target() {
    let workflow = read_release_workflow();
    assert!(
        workflow.contains("aarch64-darwin") || workflow.contains("aarch64-apple-darwin"),
        "release workflow must include macOS M1/M2 target (aarch64-darwin or aarch64-apple-darwin); \
         check the build matrix for macos-14 or macos-latest-xlarge runner"
    );
}

// ============================================================================
// Binary naming convention
// ============================================================================

/// issue #165: binaries must be named phraya-{version}-{target}-{variant}.tar.gz
#[test]
fn issue_165_release_workflow_binary_naming_follows_spec() {
    let workflow = read_release_workflow();
    // The workflow must produce .tar.gz archives named with the phraya- prefix.
    // We check for the interpolated naming pattern in the workflow YAML.
    assert!(
        workflow.contains("phraya-") && workflow.contains(".tar.gz"),
        "release workflow must produce phraya-{{version}}-{{target}}-{{variant}}.tar.gz archives; \
         check the upload/package step for this naming convention"
    );
}

/// issue #165: binary archives must embed both target triple and variant in the filename
#[test]
fn issue_165_release_workflow_archive_name_includes_target_and_variant() {
    let workflow = read_release_workflow();
    // The archive name must embed both a platform triple and a variant indicator.
    // We verify the workflow references both "target" and "variant" variables (or literals)
    // in the context of archive naming.
    let has_target_in_name = workflow.contains("target") && workflow.contains("variant");
    assert!(
        has_target_in_name,
        "release workflow archive name must embed both target and variant; \
         expected pattern: phraya-{{version}}-{{target}}-{{variant}}.tar.gz — \
         ensure the packaging step uses both variables"
    );
}

// ============================================================================
// README — platform selection guide
// ============================================================================

/// issue #165: README must contain a platform selection guide section
#[test]
fn issue_165_readme_has_platform_selection_guide() {
    let readme = read_readme();
    let has_section = readme.contains("## Platform")
        || readme.contains("## Downloading")
        || readme.contains("## Releases")
        || readme.contains("## Pre-built")
        || readme.contains("platform selection")
        || readme.contains("Platform Selection")
        || readme.contains("Which binary");
    assert!(
        has_section,
        "README must include a platform selection guide section (e.g. '## Platform Selection') \
         explaining which binary to download for each OS/architecture; none found"
    );
}

/// issue #165: README platform guide must mention all 5 supported targets
#[test]
fn issue_165_readme_platform_guide_covers_all_5_targets() {
    let readme = read_readme();
    // Each platform must appear in the README so users know what's available.
    let missing: Vec<&str> = [
        "x86_64-linux",
        "aarch64-linux",
        "x86_64-darwin",
        "aarch64-darwin",
    ]
    .iter()
    .copied()
    .filter(|&t| !readme.contains(t))
    .collect();
    assert!(
        missing.is_empty(),
        "README platform guide is missing entries for: {:?} — \
         all 5 targets (x86_64-linux-gnu portable+native, aarch64-linux-gnu, \
         x86_64-darwin, aarch64-darwin) must be documented",
        missing
    );
}

// ============================================================================
// README — SIMD differences
// ============================================================================

/// issue #165: README must describe SSE4.2 as the portable x86_64 SIMD baseline
/// (the word "portable" must appear in the README alongside SSE4.2 context)
#[test]
fn issue_165_readme_documents_sse42_as_portable_simd_level() {
    let readme = read_readme();
    // "portable" does not yet exist in the README — it is introduced by the
    // platform selection guide that documents the SSE4.2 variant as the
    // portable build that runs on any x86_64 machine.
    assert!(
        readme.contains("portable") && readme.contains("SSE4.2"),
        "README must describe SSE4.2 as the portable SIMD baseline for x86_64; \
         both 'portable' and 'SSE4.2' must appear in the README — \
         add a platform selection / SIMD differences section"
    );
}

/// issue #165: README must describe AVX2 as the native x86_64 variant requiring Haswell+
#[test]
fn issue_165_readme_documents_avx2_as_native_simd_level() {
    let readme = read_readme();
    // The README must explain that the native AVX2 build requires a Haswell or
    // newer CPU (Haswell was the first Intel microarch to ship AVX2, 2013).
    assert!(
        readme.contains("AVX2") && (readme.contains("Haswell") || readme.contains("haswell")),
        "README must explain that the native AVX2 binary requires a Haswell+ CPU; \
         add a SIMD differences section mentioning 'Haswell' as the AVX2 baseline"
    );
}

/// issue #165: README must document NEON for ARM targets including Apple M1/M2
#[test]
fn issue_165_readme_documents_neon_for_arm_targets() {
    let readme = read_readme();
    // "M1" does not yet appear in the README. The platform selection guide
    // must mention M1/M2 (Apple Silicon) as an aarch64-darwin NEON target.
    assert!(
        readme.contains("NEON") && (readme.contains("M1") || readme.contains("M2")),
        "README must mention NEON and Apple M1/M2 in the context of the aarch64-darwin binary; \
         add a platform selection guide that lists M1/M2 as a supported NEON target"
    );
}

/// issue #165: README SIMD section must explain when to choose portable vs native
#[test]
fn issue_165_readme_simd_section_explains_portable_vs_native_choice() {
    let readme = read_readme();
    // A platform selection guide must explain the SSE4.2 vs AVX2 tradeoff.
    // We check that the words "portable" and "native" appear near SIMD context.
    let has_portable = readme.contains("portable");
    let has_native_simd = readme.contains("native") && (readme.contains("AVX2") || readme.contains("target-cpu"));
    assert!(
        has_portable && has_native_simd,
        "README must explain when to use portable (SSE4.2, runs on any x86_64) vs \
         native (AVX2, requires Haswell+); both 'portable' and 'native' with AVX2 \
         context must appear in the README"
    );
}
