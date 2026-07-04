/// Issue #184: CLI strategy flag tests for recall-axis ladder.
/// Tests the --strategy flag parsing and behavior with renamed (exact -> sensitive) and
/// K-capped (K=1/5/∞) anchor selection.
use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};
use phraya_core::types::Sequence;
use phraya_io::plan::{write_plan, PhrayaPlan, UseCase};
use std::collections::HashMap;
use std::path::PathBuf;
use tempfile::TempDir;

/// Helper to get phraya-cli manifest path for cargo run commands (same pattern as
/// issue_181_preset_rename.rs).
fn get_manifest_path() -> PathBuf {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    std::path::Path::new(&manifest_dir).join("Cargo.toml")
}

fn make_plan_with_fasta(fasta_path: &str) -> PhrayaPlan {
    PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec![fasta_path.to_string()],
        "2026-06-06T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    )
}

fn make_seqs() -> (Sequence, Sequence) {
    let mut query = vec![b'A'; 100];
    let mut target = vec![b'A'; 200];
    query[50] = b'T';
    target[50] = b'C';
    (
        Sequence::new(query, None, "read1".to_string(), None),
        Sequence::new(target, None, "ref".to_string(), None),
    )
}

// ============================================================================
// CLI PARSING TESTS: Verify --strategy accepts/rejects correct values
// ============================================================================

/// issue #184: AlignConfig accepts strategy sensitive and has correct radius
#[test]
fn issue_184_cli_align_accepts_strategy_sensitive() {
    let config = AlignConfig::new(Strategy::Sensitive);
    assert_eq!(config.coverage_window_radius, 25);
    assert_eq!(config.strategy, Strategy::Sensitive);
}

/// issue #184: AlignConfig accepts strategy balanced and has correct radius
#[test]
fn issue_184_cli_align_accepts_strategy_balanced() {
    let config = AlignConfig::new(Strategy::Balanced);
    assert_eq!(config.coverage_window_radius, 50);
}

/// issue #184: AlignConfig accepts strategy fast and has correct radius
#[test]
fn issue_184_cli_align_accepts_strategy_fast() {
    let config = AlignConfig::new(Strategy::Fast);
    assert_eq!(config.coverage_window_radius, 150);
}

/// issue #184: `phraya align --strategy exact` is rejected by the real CLI parser
/// (no longer a valid strategy). Spawns the actual binary (pattern from
/// issue_181_preset_rename.rs) rather than simulating the match locally -- a local
/// simulation would pass whether or not phraya-cli/src/main.rs was ever updated.
#[test]
fn issue_184_cli_align_rejects_exact_strategy() {
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "align",
            "nonexistent.phrayaplan",
            "--strategy",
            "exact",
        ])
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        !output.status.success(),
        "phraya align --strategy exact must fail"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("unknown strategy"),
        "stderr must mention 'unknown strategy': {stderr}"
    );
    assert!(
        stderr.contains("sensitive"),
        "stderr must list 'sensitive' as a valid option: {stderr}"
    );
    assert!(
        !stderr.contains("expected fast, balanced, or exact"),
        "stderr must not still advertise 'exact' as valid: {stderr}"
    );
}

/// issue #184: an invalid --strategy value returns the updated error message via the
/// real CLI parser.
#[test]
fn issue_184_cli_align_rejects_invalid_strategy_with_sensitive_option() {
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "align",
            "nonexistent.phrayaplan",
            "--strategy",
            "invalid_strategy",
        ])
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        !output.status.success(),
        "phraya align --strategy invalid_strategy must fail"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("unknown strategy"),
        "stderr must mention 'unknown strategy': {stderr}"
    );
    assert!(
        stderr.contains("sensitive"),
        "stderr must list 'sensitive' as valid option: {stderr}"
    );
    assert!(
        !stderr.contains("expected fast, balanced, or exact"),
        "stderr must not still advertise 'exact' as valid: {stderr}"
    );
}

/// issue #184: `phraya align --strategy sensitive` is accepted by the real CLI parser.
#[test]
fn issue_184_cli_align_accepts_sensitive_strategy_string() {
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "align",
            "nonexistent.phrayaplan",
            "--strategy",
            "sensitive",
        ])
        .output()
        .expect("Failed to execute phraya align");

    // The strategy parses fine; the run still fails, but for a *different* reason
    // (missing plan file), never "unknown strategy".
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("unknown strategy"),
        "stderr must not reject 'sensitive' as an unknown strategy: {stderr}"
    );
}

/// issue #184: default AlignConfig is still balanced
#[test]
fn issue_184_cli_align_uses_default_strategy_without_flag() {
    let config = AlignConfig::default();
    assert_eq!(
        config.coverage_window_radius, 50,
        "default strategy must be balanced (±50bp)"
    );
    assert_eq!(config.strategy, Strategy::Balanced);
}

/// issue #184: sensitive strategy produces ±25bp coverage window (same as old exact)
#[test]
fn issue_184_cli_sensitive_produces_narrow_coverage_window() {
    let config = AlignConfig::new(Strategy::Sensitive);
    assert_eq!(
        config.coverage_window_radius, 25,
        "sensitive strategy must use ±25bp window (same as old exact)"
    );
}

// ============================================================================
// ALIGNMENT BEHAVIOR TESTS: Verify strategies work correctly via align_task_with_config
// ============================================================================

/// issue #184: all three strategies produce different coverage windows
#[test]
fn issue_184_cli_different_strategies_produce_different_windows() {
    let dir = TempDir::new().unwrap();
    let tmp = dir.path().join("plan.phrayaplan");
    let plan = make_plan_with_fasta(tmp.to_str().unwrap());
    write_plan(&tmp, &plan).unwrap();
    let plan = phraya_io::plan::read_plan(&tmp).unwrap();

    let (query, target) = make_seqs();

    let result_fast =
        align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Fast))
            .expect("fast alignment should succeed");
    let result_balanced = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Balanced),
    )
    .expect("balanced alignment should succeed");
    let result_sensitive = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Sensitive),
    )
    .expect("sensitive alignment should succeed");

    assert!(
        !result_fast.variants.is_empty(),
        "fast should produce variants"
    );
    assert!(
        !result_balanced.variants.is_empty(),
        "balanced should produce variants"
    );
    assert!(
        !result_sensitive.variants.is_empty(),
        "sensitive should produce variants"
    );

    let fast_window = result_fast.variants[0].local_coverage().len();
    let balanced_window = result_balanced.variants[0].local_coverage().len();
    let sensitive_window = result_sensitive.variants[0].local_coverage().len();

    // Verify the hierarchy:
    // fast (±150) > balanced (±50) > sensitive (±25)
    assert!(
        fast_window > balanced_window,
        "fast strategy window ({fast_window}) must be larger than balanced ({balanced_window})"
    );
    assert!(
        balanced_window > sensitive_window,
        "balanced strategy window ({balanced_window}) must be larger than sensitive ({sensitive_window})"
    );
}

/// issue #184: sensitive strategy produces the same window as old exact strategy
#[test]
fn issue_184_cli_sensitive_window_matches_old_exact() {
    let config = AlignConfig::new(Strategy::Sensitive);
    // The window radius should match the old exact strategy (±25bp).
    assert_eq!(
        config.coverage_window_radius, 25,
        "sensitive must use the same window radius as the old exact strategy (±25bp)"
    );
}

/// issue #184: override coverage_window decouples radius from strategy
#[test]
fn issue_184_cli_coverage_window_override_is_orthogonal() {
    // Verify the base radii are strategy-dependent
    assert_eq!(AlignConfig::new(Strategy::Fast).coverage_window_radius, 150);
    assert_eq!(
        AlignConfig::new(Strategy::Sensitive).coverage_window_radius,
        25
    );

    // Override allows decoupling
    let cfg_fast = AlignConfig::new(Strategy::Fast).with_coverage_window_radius(25);
    assert_eq!(cfg_fast.coverage_window_radius, 25);
    assert_eq!(
        cfg_fast.strategy,
        Strategy::Fast,
        "override must not change the strategy"
    );

    let cfg_sensitive = AlignConfig::new(Strategy::Sensitive).with_coverage_window_radius(150);
    assert_eq!(cfg_sensitive.coverage_window_radius, 150);
    assert_eq!(
        cfg_sensitive.strategy,
        Strategy::Sensitive,
        "override must not change the strategy"
    );
}

// ============================================================================
// ANCHOR CAP (K) TESTS: Verify K=1/5/∞ behavior
// ============================================================================

/// issue #184: fast strategy (K=1) produces single-anchor placements on multi-mapping loci
#[test]
fn issue_184_cli_fast_strategy_k1_single_placement() {
    let dir = TempDir::new().unwrap();
    let tmp = dir.path().join("plan.phrayaplan");
    let plan = make_plan_with_fasta(tmp.to_str().unwrap());
    write_plan(&tmp, &plan).unwrap();
    let plan = phraya_io::plan::read_plan(&tmp).unwrap();

    // Create a 50bp repetitive read that matches at two tandem positions.
    let unit = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // ~56bp
    let query = Sequence::new(unit.to_vec(), None, "read1".to_string(), None);

    // Target with tandem repeats
    let mut target_bases = unit.to_vec();
    target_bases.extend_from_slice(unit);
    target_bases.extend_from_slice(unit);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);

    let result = align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Fast))
        .expect("fast alignment should succeed");

    assert_eq!(
        result.query_positions.len(),
        1,
        "Fast (K=1) must report exactly 1 placement on tandem repeats, got {}",
        result.query_positions.len()
    );
}

/// issue #184: balanced strategy (K=5) preserves multi-mapping on ambiguous loci
#[test]
fn issue_184_cli_balanced_strategy_k5_preserves_multimapping() {
    let dir = TempDir::new().unwrap();
    let tmp = dir.path().join("plan.phrayaplan");
    let plan = make_plan_with_fasta(tmp.to_str().unwrap());
    write_plan(&tmp, &plan).unwrap();
    let plan = phraya_io::plan::read_plan(&tmp).unwrap();

    let unit = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let query = Sequence::new(unit.to_vec(), None, "read1".to_string(), None);

    let mut target_bases = unit.to_vec();
    target_bases.extend_from_slice(unit);
    target_bases.extend_from_slice(unit);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);

    let result = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Balanced),
    )
    .expect("balanced alignment should succeed");

    assert!(
        result.query_positions.len() >= 2,
        "Balanced (K=5) must preserve multi-mapping, got {} placements",
        result.query_positions.len()
    );
}

/// issue #184: sensitive strategy (K=∞) reports all seed-derived anchors
#[test]
fn issue_184_cli_sensitive_strategy_kinf_all_anchors() {
    let dir = TempDir::new().unwrap();
    let tmp = dir.path().join("plan.phrayaplan");
    let plan = make_plan_with_fasta(tmp.to_str().unwrap());
    write_plan(&tmp, &plan).unwrap();
    let plan = phraya_io::plan::read_plan(&tmp).unwrap();

    let unit = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let query = Sequence::new(unit.to_vec(), None, "read1".to_string(), None);

    // Multiple tandem repeats to ensure K=∞ reports all
    let mut target_bases = Vec::new();
    for _ in 0..4 {
        target_bases.extend_from_slice(unit);
    }
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);

    let result = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Sensitive),
    )
    .expect("sensitive alignment should succeed");

    assert!(
        result.query_positions.len() >= 3,
        "Sensitive (K=∞) must report all seed-derived anchors (≥3 for 4× tandem), got {}",
        result.query_positions.len()
    );
}

// ============================================================================
// FACTORY METHOD TESTS: Verify AlignConfig builder API
// ============================================================================

/// issue #184: AlignConfig::fast() factory method works
#[test]
fn issue_184_cli_config_factory_fast() {
    let config = AlignConfig::fast();
    assert_eq!(config.strategy, Strategy::Fast);
    assert_eq!(config.coverage_window_radius, 150);
}

/// issue #184: AlignConfig::balanced() factory method works
#[test]
fn issue_184_cli_config_factory_balanced() {
    let config = AlignConfig::balanced();
    assert_eq!(config.strategy, Strategy::Balanced);
    assert_eq!(config.coverage_window_radius, 50);
}

/// issue #184: AlignConfig::sensitive() factory method works (replaces exact())
#[test]
fn issue_184_cli_config_factory_sensitive() {
    let config = AlignConfig::sensitive();
    assert_eq!(config.strategy, Strategy::Sensitive);
    assert_eq!(config.coverage_window_radius, 25);
}

// ============================================================================
// NOMENCLATURE: real "exact" removal is exercised by
// issue_184_cli_align_rejects_exact_strategy above (subprocess-based -- the CLI
// genuinely refuses "exact" at runtime, not just a type-level presence check).
// ============================================================================
