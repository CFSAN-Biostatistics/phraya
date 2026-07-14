//! RED acceptance tests for issue #200: `align_read` reuses a read's sketch by content
//! hash (`PhrayaPlan::read_sketches`), instead of always recomputing it from scratch.
//!
//! `get_effective_sketch` (phraya-align/src/executor.rs) already reuses a *target's*
//! stored sketch by looking it up in `plan.kmer_index`/`dense_kmer_index` keyed by
//! sequence ID. Issue #200 asks for the same reuse for *reads*, but keyed by content
//! hash (`PhrayaPlan::read_sketches: HashMap<u64, MinimizerSketch>`, added in #200's
//! first slice, PR #225) rather than by ID — a read surviving into a later pipeline
//! stage isn't guaranteed to keep the same ID, but its content hash is stable.
//!
//! Proof technique: store a deliberately wrong/empty `MinimizerSketch` in
//! `plan.read_sketches` under the read's actual content hash (computed via
//! `phraya_io::plan::read_content_hash`). If `align_read` reuses it instead of
//! recomputing a correct sketch from the read's real bases, seeding against the target
//! finds zero seeds (the empty sketch shares no minimizers with anything) and the read
//! fails to place with a real, non-trivial edit distance — a directly observable
//! difference from correctly recomputing the sketch, which would place the read cleanly
//! against a matching target.

use phraya_align::executor::{align_read, AlignConfig, Strategy, TargetContext};
use phraya_core::types::{sketch_sequence_default, MinimizerSketch, Sequence};
use phraya_io::plan::{read_content_hash, PhrayaPlan, UseCase};
use std::collections::HashMap;

fn make_plan() -> PhrayaPlan {
    PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["test".to_string()],
        "2026-07-13T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    )
}

fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

/// A read whose sketch is present in `plan.read_sketches` under its own content hash,
/// but whose ID does not appear in `plan.kmer_index` at all, must still have its sketch
/// reused — content hash, not sequence ID, is the lookup key for read sketches.
///
/// The stored sketch here is correct (matches the read's real content), so if reuse
/// works, the read places cleanly. This is the positive-control half of the proof; see
/// `read_sketch_reuse_prefers_stored_over_recompute` for the half that distinguishes
/// "reused" from "recomputed anyway."
#[test]
fn read_sketch_is_looked_up_by_content_hash_not_sequence_id() {
    let target_bases = diverse_dna(2_000, 7);
    let target = Sequence::new(target_bases.clone(), None, "ref".to_string(), None);
    let mut plan = make_plan();

    let read_bases = target_bases[500..650].to_vec();
    let read = Sequence::new(
        read_bases.clone(),
        None,
        "read_with_untracked_id".to_string(),
        None,
    );

    // Deliberately do NOT add anything to plan.kmer_index for this read's ID — the only
    // way align_read can find a stored sketch is via read_sketches keyed by content hash.
    let hash = read_content_hash(&read_bases);
    let correct_sketch = sketch_sequence_default(&read);
    plan.read_sketches.insert(hash, correct_sketch);

    let config = AlignConfig::new(Strategy::Balanced);
    let ctx = TargetContext::build(&target, &plan, config.strategy);

    let result = align_read(&ctx, &read, &plan, &config, None);
    assert!(
        result.is_some(),
        "a read with a correct sketch stored under its content hash must place, \
         even though its sequence ID has no entry in plan.kmer_index"
    );
}

/// If `align_read` truly reuses the stored read sketch (rather than silently falling
/// back to recomputing from the read's real bases), poisoning the stored sketch with an
/// empty one must cause the read to fail to place — seeding against an empty sketch
/// finds zero seeds. This is the test that actually distinguishes "reuse happened" from
/// "recompute happened anyway and the stored sketch was ignored."
#[test]
fn read_sketch_reuse_prefers_stored_over_recompute() {
    let target_bases = diverse_dna(2_000, 11);
    let target = Sequence::new(target_bases.clone(), None, "ref".to_string(), None);
    let mut plan = make_plan();

    let read_bases = target_bases[800..950].to_vec();
    let read = Sequence::new(read_bases.clone(), None, "read_poisoned".to_string(), None);

    let hash = read_content_hash(&read_bases);
    let empty_sketch = MinimizerSketch {
        minimizers: vec![],
        k: 21,
        w: 11,
    };
    plan.read_sketches.insert(hash, empty_sketch);

    let config = AlignConfig::new(Strategy::Balanced);
    let ctx = TargetContext::build(&target, &plan, config.strategy);

    // Sanity check: without poisoning, this exact read/target pair places cleanly.
    let mut clean_plan = make_plan();
    let clean_ctx = TargetContext::build(&target, &clean_plan, config.strategy);
    let clean_result = align_read(&clean_ctx, &read, &clean_plan, &config, None);
    assert!(
        clean_result.is_some(),
        "control check failed: this read/target pair should place without any stored sketch"
    );
    clean_plan.read_sketches.clear(); // silence unused-mut if the above changes

    let poisoned_result = align_read(&ctx, &read, &plan, &config, None);
    assert!(
        poisoned_result.is_none(),
        "align_read must prefer the stored (poisoned, empty) read sketch over recomputing \
         a correct one — if this assertion fails, sketch reuse-by-content-hash isn't wired \
         up and align_read is recomputing from the read's real bases regardless of \
         plan.read_sketches"
    );
}

/// A read whose content hash has no entry in `plan.read_sketches` at all must fall back
/// to recomputing its sketch (via `sketch_sequence_default` or equivalent) rather than
/// failing — the reuse path is an optimization, not a requirement.
#[test]
fn read_without_stored_sketch_falls_back_to_recompute() {
    let target_bases = diverse_dna(2_000, 23);
    let target = Sequence::new(target_bases.clone(), None, "ref".to_string(), None);
    let plan = make_plan(); // read_sketches empty

    let read_bases = target_bases[300..450].to_vec();
    let read = Sequence::new(read_bases, None, "read_no_cache".to_string(), None);

    let config = AlignConfig::new(Strategy::Balanced);
    let ctx = TargetContext::build(&target, &plan, config.strategy);

    let result = align_read(&ctx, &read, &plan, &config, None);
    assert!(
        result.is_some(),
        "a read with no stored sketch must still place by recomputing its sketch fresh"
    );
}
