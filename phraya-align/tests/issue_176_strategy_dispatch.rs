//! Issue #176: alignment strategy ladder (Sensitive / Balanced / Fast).
//!
//! New semantics:
//!   - Sensitive: seeded WFA, all anchors (K=∞) — the canonical reference path.
//!   - Balanced: Myers fitting (≤500bp) with WFA fallback, top 5 anchors (K=5) — exact results, faster engine.
//!   - Fast:     low-sensitivity (seed subsampling + divergence cutoff, K=1).
//!
//! Because Myers and WFA compute identical edit distances, Sensitive and Balanced must agree
//! on the variants they call for a uniquely-mapping read. This file pins that invariant.

use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};
use phraya_core::types::Sequence;
use phraya_io::plan::{PhrayaPlan, UseCase};
use std::collections::HashMap;

fn make_plan() -> PhrayaPlan {
    PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["test".to_string()],
        "2026-06-01T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    )
}

/// A 150bp read with two SNPs against a 300bp reference window. Sensitive (WFA) and Balanced
/// (Myers) must report the same variant positions and reference/alt bases.
#[test]
fn issue_176_sensitive_and_balanced_call_identical_variants() {
    // Deterministic pseudo-random reference so seeds actually anchor the read.
    let mut state: u64 = 0xD1B54A32D192ED03;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) as u32
    };
    let bases = [b'A', b'C', b'G', b'T'];
    let read: Vec<u8> = (0..150).map(|_| bases[(next() % 4) as usize]).collect();

    // Target = read region (with 2 SNPs) + a divergent tail to ~2× length.
    let mut region = read.clone();
    region[40] = if region[40] == b'A' { b'C' } else { b'A' };
    region[110] = if region[110] == b'G' { b'T' } else { b'G' };
    let tail: Vec<u8> = (0..160).map(|_| bases[(next() % 4) as usize]).collect();
    let mut target_bases = region;
    target_bases.extend_from_slice(&tail);

    let query = Sequence::new(read, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let sensitive = align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Sensitive))
        .expect("sensitive alignment should succeed");
    let balanced =
        align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Balanced))
            .expect("balanced alignment should succeed");

    let key = |v: &phraya_core::types::VariantObservation| {
        (v.position(), v.ref_base(), {
            let mut alleles: Vec<(u8, u32)> =
                v.all_alleles().iter().map(|(b, c)| (*b, *c)).collect();
            alleles.sort();
            alleles
        })
    };
    let mut sensitive_vars: Vec<_> = sensitive.variants.iter().map(key).collect();
    let mut balanced_vars: Vec<_> = balanced.variants.iter().map(key).collect();
    sensitive_vars.sort();
    balanced_vars.sort();

    assert_eq!(
        sensitive_vars, balanced_vars,
        "Sensitive (WFA) and Balanced (Myers) must call identical variants for a uniquely-mapping read"
    );
}

/// The strategy preset sets a default coverage-window radius, but it can be overridden
/// orthogonally — choosing an algorithm should not lock the annotation width.
#[test]
fn issue_176_coverage_window_radius_override_is_orthogonal() {
    // Preset default still applies when not overridden.
    assert_eq!(AlignConfig::new(Strategy::Fast).coverage_window_radius, 150);
    assert_eq!(AlignConfig::new(Strategy::Sensitive).coverage_window_radius, 25);

    // Override decouples the radius from the strategy.
    let cfg = AlignConfig::new(Strategy::Fast).with_coverage_window_radius(10);
    assert_eq!(cfg.coverage_window_radius, 10);
    assert_eq!(cfg.strategy, Strategy::Fast, "override must not change the strategy");
}

/// Deterministic pseudo-random DNA of a given length.
fn random_dna(seed: u64, len: usize) -> Vec<u8> {
    let mut state = seed;
    let bases = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[((state >> 33) % 4) as usize]
        })
        .collect()
}

/// Fast sacrifices sensitivity for speed: a read that only aligns at high divergence is
/// dropped under Fast (divergence cutoff) but still called under Balanced.
#[test]
fn issue_176_fast_drops_divergent_read_that_balanced_keeps() {
    let read = random_dna(0x1234_5678, 100);

    // Target: the read with a dense block of 40 substitutions (≈40% divergence over the
    // read), followed by an unrelated tail. The flanking matches still seed an anchor at 0.
    let mut region = read.clone();
    for i in 50..90 {
        region[i] = match region[i] {
            b'A' => b'C',
            b'C' => b'A',
            b'G' => b'T',
            _ => b'G',
        };
    }
    let tail = random_dna(0x9999, 120);
    let mut target_bases = region;
    target_bases.extend_from_slice(&tail);

    let query = Sequence::new(read, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let balanced =
        align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Balanced));
    let fast = align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Fast));

    assert!(
        balanced.is_some(),
        "Balanced should align the divergent read (no cutoff)"
    );
    assert!(
        fast.is_none(),
        "Fast should drop the read whose best alignment exceeds the divergence cutoff"
    );
}

/// Fast must still produce the same primary variant calls as Balanced for a clean,
/// uniquely-mapping read well within the divergence cutoff.
#[test]
fn issue_176_fast_matches_balanced_on_clean_read() {
    let read = random_dna(0xABCD_EF01, 150);

    // Two isolated SNPs, low divergence — comfortably under Fast's cutoff.
    let mut region = read.clone();
    region[37] = if region[37] == b'A' { b'C' } else { b'A' };
    region[120] = if region[120] == b'G' { b'T' } else { b'G' };
    let tail = random_dna(0x5555, 170);
    let mut target_bases = region;
    target_bases.extend_from_slice(&tail);

    let query = Sequence::new(read, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let balanced =
        align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Balanced))
            .expect("balanced should align");
    let fast = align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Fast))
        .expect("fast should align a clean read");

    let positions = |r: &phraya_align::executor::AlignmentResult| {
        let mut p: Vec<(u32, u8)> =
            r.variants.iter().map(|v| (v.position(), v.ref_base())).collect();
        p.sort();
        p
    };
    assert_eq!(
        positions(&fast),
        positions(&balanced),
        "Fast and Balanced must agree on variant positions for a clean uniquely-mapping read"
    );
}

/// Fast's seed subsampling under-reports multi-mapping: against a tandem-duplicated
/// target, Sensitive records both equally-good positions while Fast keeps only the best
/// single anchor. This pins the documented sensitivity tradeoff.
#[test]
fn issue_176_fast_underreports_multimapping_vs_sensitive() {
    let unit = random_dna(0x0FAC_E001, 80);
    // Tandem duplication: the read matches at target offset 0 and offset 80 equally well.
    let mut target_bases = unit.clone();
    target_bases.extend_from_slice(&unit);
    target_bases.extend_from_slice(&random_dna(0x7777, 60));

    let query = Sequence::new(unit, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let sensitive = align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Sensitive))
        .expect("sensitive should align");
    let fast = align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Fast))
        .expect("fast should align");

    assert!(
        sensitive.query_positions.len() >= 2,
        "Sensitive should record both copies of the tandem duplication (multi-mapping), got {}",
        sensitive.query_positions.len()
    );
    assert_eq!(
        fast.query_positions.len(),
        1,
        "Fast keeps only the single best-voted anchor (multi-mapping under-reported), got {}",
        fast.query_positions.len()
    );
}
