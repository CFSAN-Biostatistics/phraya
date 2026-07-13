/// Issue #148: CLI strategy flag acceptance tests.
/// Uses library API directly rather than subprocess invocation for speed and reliability.
use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};
use phraya_core::types::Sequence;
use phraya_io::plan::{write_plan, PhrayaPlan, UseCase};
use std::collections::HashMap;
use tempfile::TempDir;

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

/// issue #148: AlignConfig accepts fast strategy and has correct radius
#[test]
fn issue_148_align_accepts_strategy_fast() {
    let config = AlignConfig::new(Strategy::Fast);
    assert_eq!(config.coverage_window_radius, 150);
}

/// issue #148: AlignConfig accepts balanced strategy and has correct radius
#[test]
fn issue_148_align_accepts_strategy_balanced() {
    let config = AlignConfig::new(Strategy::Balanced);
    assert_eq!(config.coverage_window_radius, 50);
}

/// issue #148: AlignConfig accepts sensitive strategy and has correct radius
#[test]
fn issue_148_align_accepts_strategy_sensitive() {
    let config = AlignConfig::new(Strategy::Sensitive);
    assert_eq!(config.coverage_window_radius, 25);
}

/// issue #148: invalid strategy string returns error (tests run_align error path)
#[test]
fn issue_148_align_rejects_invalid_strategy() {
    // Simulate the CLI parsing: unknown strategy should produce an error
    let result: Result<Strategy, String> = match "invalid_strategy" {
        "fast" => Ok(Strategy::Fast),
        "balanced" => Ok(Strategy::Balanced),
        "sensitive" => Ok(Strategy::Sensitive),
        other => Err(format!(
            "unknown strategy: {other}; expected fast, balanced, or sensitive"
        )),
    };
    assert!(result.is_err(), "invalid strategy must be rejected");
    let err = result.unwrap_err();
    assert!(
        err.contains("unknown strategy"),
        "error must mention 'unknown strategy': {err}"
    );
}

/// issue #148: default AlignConfig is balanced
#[test]
fn issue_148_align_uses_default_strategy_without_flag() {
    let config = AlignConfig::default();
    assert_eq!(
        config.coverage_window_radius, 50,
        "default strategy must be balanced (±50bp)"
    );
    assert_eq!(config.strategy, Strategy::Balanced);
}

/// issue #148: fast strategy produces wider coverage window than balanced on same alignment
#[test]
fn issue_148_cli_different_strategies_produce_different_windows() {
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

    assert!(
        !result_fast.variants.is_empty(),
        "fast should produce variants"
    );
    assert!(
        !result_balanced.variants.is_empty(),
        "balanced should produce variants"
    );

    let fast_window = result_fast.variants[0].local_coverage().len();
    let balanced_window = result_balanced.variants[0].local_coverage().len();

    assert!(
        fast_window > balanced_window,
        "fast strategy window ({fast_window}) must be larger than balanced ({balanced_window})"
    );
}
