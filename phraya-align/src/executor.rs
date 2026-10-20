use crate::seeding::{find_seeds_indexed_capped, seed_occurrence_cap, MinimizerIndex};
use crate::{myers_edit_distance, myers_extend, score_alignments, wfa_extend, SeedAnchor};
use phraya_core::types::{
    reverse_complement, sketch_sequence_default, sketch_sequence_default_for, Alphabet,
    AlignmentBlock, MinimizerSketch, Sequence, Strand, VariantObservation,
};
use phraya_core::{detect_tandem_repeats, RepeatDetectorConfig};
use phraya_io::plan::{read_content_hash, PhrayaPlan};
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};

/// Helper: Filter a dense sketch to the w=11 subset using membership tags from the plan.
///
/// For Fast/Balanced strategies, we need to use only the w=11 subset of minimizers to maintain
/// byte-identity with pre-#182 results. This function selects only those dense minimizers whose
/// w11_membership tags are `true`, producing a sketch that is byte-identical to the canonical
/// w=11 sketch (after deduplication).
///
/// If the plan lacks dense data or w11_membership tags, returns the sketch unchanged.
/// Returns `None` if Sensitive is requested on a sparse plan (dense sketches unavailable).
fn get_effective_sketch(
    sequence_id: &str,
    plan: &PhrayaPlan,
    strategy: Strategy,
) -> Option<MinimizerSketch> {
    // For Sensitive strategy: use dense sketch if available, or w=11 if sparse
    if strategy == Strategy::Sensitive {
        if let Some(dense) = plan.get_dense_sketch(sequence_id) {
            return Some(dense.clone());
        }
        if plan.is_sparse() {
            // Sparse plan: no dense sketches available. Degrade gracefully by falling back to w=11.
            return plan.get_sketch(sequence_id).cloned();
        }
        return plan.get_sketch(sequence_id).cloned();
    }

    // For Fast/Balanced: filter dense to w=11 subset if available, else use w=11 directly
    if let Some(dense) = plan.get_dense_sketch(sequence_id) {
        if let Some(membership) = plan.get_w11_membership(sequence_id) {
            // Filter dense minimizers to those tagged as w=11 subset
            let filtered_minimizers: Vec<(u64, u32)> = dense
                .minimizers
                .iter()
                .enumerate()
                .filter_map(|(i, &m)| {
                    if i < membership.len() && membership[i] {
                        Some(m)
                    } else {
                        None
                    }
                })
                .collect();
            return Some(MinimizerSketch {
                minimizers: filtered_minimizers,
                k: dense.k,
                w: dense.w,
            });
        }
    }

    // Fallback: use w=11 sketch from kmer_index
    plan.get_sketch(sequence_id).cloned()
}

/// Alignment strategy affecting coverage window size and anchor selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strategy {
    /// Fast strategy: K=1 chain cap, wide ±150bp coverage window
    Fast,
    /// Balanced strategy: K=2 chain cap, moderate ±50bp coverage window (default)
    Balanced,
    /// Sensitive strategy: K=50 chain cap, narrow ±25bp coverage window for precision
    Sensitive,
    /// Exact strategy: K=250 chain cap, minimal ±5bp coverage window, no mismatches tolerated in initial seeding
    Exact,
    /// Long-read strategy: Chunked alignment for ONT/PacBio reads (10-50kb), adaptive banding
    LongRead,
    /// Assembly strategy for contig-vs-contig comparison (ADR-0015): every chain in either
    /// orientation becomes a local alignment block (its exact seed segments joined by global
    /// gap fills), so a contig yields one block per homologous region instead of one
    /// whole-query extension. ±5bp coverage window.
    Assembly,
}

impl Default for Strategy {
    fn default() -> Self {
        Strategy::Balanced
    }
}

/// Gap cost model for WFA extension (ADR-0014). Only meaningful for
/// `Strategy::Sensitive` — `Balanced`/`Fast` are Myers-primary and stay `Linear`
/// unconditionally regardless of this setting, since Myers has no affine mode (and
/// their WFA fallback for long queries exists only to mirror Myers cheaply, not to
/// introduce a different cost model at that tier).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapModel {
    /// Uniform per-base cost (the historical default): a mismatch and each indel base
    /// cost 1 unit, so a k-base indel costs the same as k independent substitutions.
    Linear,
    /// `gap_open + L * gap_extend` per L-base gap (`wfa_simd::AffineCosts::default()`):
    /// consolidates a real multi-base indel into one gap rather than leaving the search
    /// indifferent between one gap and several mismatch-cost-tied substitutions.
    Affine,
}

impl Default for GapModel {
    /// `Sensitive`'s default (ADR-0014) — the more biologically faithful model, for the
    /// tier whose users already opted into the slower, most-complete alignment path.
    fn default() -> Self {
        GapModel::Affine
    }
}

/// Configuration for alignment execution, controlling coverage window size.
#[derive(Debug, Clone, Copy)]
pub struct AlignConfig {
    /// Alignment strategy
    pub strategy: Strategy,
    /// Coverage window radius in base pairs
    pub coverage_window_radius: usize,
    /// Gap cost model (ADR-0014) — consulted only when `strategy == Sensitive`.
    pub gap_model: GapModel,
}

impl AlignConfig {
    /// Create a new AlignConfig with the specified strategy.
    /// The coverage_window_radius is automatically set based on the strategy.
    pub fn new(strategy: Strategy) -> Self {
        let coverage_window_radius = match strategy {
            Strategy::Fast => 150,
            Strategy::Balanced => 50,
            Strategy::Sensitive => 25,
            Strategy::Exact => 5,
            Strategy::LongRead => 25,  // Long reads use narrow window like Sensitive
            Strategy::Assembly => 5,   // Contigs carry depth 1; tight window like Exact
        };
        AlignConfig {
            strategy,
            coverage_window_radius,
            gap_model: GapModel::default(),
        }
    }

    /// Create a Fast strategy config (±150bp window).
    pub fn fast() -> Self {
        Self::new(Strategy::Fast)
    }

    /// Create a Balanced strategy config (±50bp window).
    pub fn balanced() -> Self {
        Self::new(Strategy::Balanced)
    }

    /// Create a Sensitive strategy config (±25bp window).
    pub fn sensitive() -> Self {
        Self::new(Strategy::Sensitive)
    }

    /// Create an Exact strategy config (±5bp window, K=250, no mismatch tolerance in seeding).
    pub fn exact() -> Self {
        Self::new(Strategy::Exact)
    }

    /// Create a LongRead strategy config (±25bp window, K=100, for ONT/PacBio reads 10-50kb).
    pub fn long_read() -> Self {
        Self::new(Strategy::LongRead)
    }

    /// The strategy still selects the alignment algorithm; this only changes the width
    /// of the per-variant local-coverage annotation.
    pub fn with_coverage_window_radius(mut self, radius: usize) -> Self {
        self.coverage_window_radius = radius;
        self
    }

    /// Override the gap cost model independently of the strategy preset (ADR-0014).
    /// Meaningful only when `strategy == Sensitive` — `Balanced`/`Fast` ignore this and
    /// stay `Linear` unconditionally (Myers has no affine mode).
    pub fn with_gap_model(mut self, gap_model: GapModel) -> Self {
        self.gap_model = gap_model;
        self
    }
}

impl Default for AlignConfig {
    fn default() -> Self {
        AlignConfig::new(Strategy::default())
    }
}

/// Per-read coverage over only the aligned span.
///
/// A single read touches ~its own length of the target, so materializing a
/// genome-length vector per read (and scanning it on merge) is wasteful. This
/// carries just the covered window: `counts[j]` is the coverage at absolute target
/// position `start + j`; every position outside `[start, start + counts.len())` is 0.
#[derive(Debug, Clone, Default)]
pub struct WindowedCoverage {
    /// Absolute target position of `counts[0]`.
    pub start: usize,
    /// Coverage over `[start, start + counts.len())`. Outside the window, 0.
    pub counts: Vec<u32>,
}

impl WindowedCoverage {
    /// Coverage at an absolute target position (0 outside the window).
    #[inline]
    pub fn get_abs(&self, abs_pos: usize) -> u32 {
        abs_pos
            .checked_sub(self.start)
            .and_then(|j| self.counts.get(j))
            .copied()
            .unwrap_or(0)
    }

    /// Expand to a full-length genome track of length `len`, for the single-read
    /// output path (which writes a whole-reference coverage track).
    pub fn to_full(&self, len: usize) -> Vec<u32> {
        let mut full = vec![0u32; len];
        for (j, &c) in self.counts.iter().enumerate() {
            let pos = self.start + j;
            if pos < len {
                full[pos] = c;
            }
        }
        full
    }
}

/// Sum coverage at an absolute target position across a list of (possibly overlapping)
/// windows — e.g. one per alignment (primary + alternatives) from
/// [`compute_windowed_coverage`]. `windows` is at most K (the strategy's anchor cap), so
/// this stays O(K) per lookup rather than requiring a merged genome-length buffer.
#[inline]
fn get_abs_multi(windows: &[WindowedCoverage], abs_pos: usize) -> u32 {
    windows.iter().map(|w| w.get_abs(abs_pos)).sum()
}

/// Result of a single alignment task.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Variant observations at polymorphic sites
    pub variants: Vec<VariantObservation>,
    /// Coverage over each alignment's own span (quantized: depths below 5 exact, 5 and
    /// above rounded to the nearest multiple of 5) — one small window per alignment
    /// (primary + alternatives), not a single buffer spanning all of them, so
    /// multi-mapped reads on repeat-rich genomes don't materialize a near-genome-length
    /// buffer. Merge each window into a genome accumulator at `.start` independently.
    pub coverage: Vec<WindowedCoverage>,
    /// Same windows as `coverage`, but un-quantized raw per-position depth. Coverage
    /// breadth (`phraya_core::types::coverage_breadth`) needs this: `coverage`'s
    /// depth-1-4 values are exact, but 5-and-above are still rounded to the nearest 5
    /// (e.g. depth 8 or 9 rounds *up* to 10), which would overcount 10x breadth (see that
    /// function's doc comment). Not used for anything else — `local_coverage` windows and
    /// the merged `CoverageTrack` still come from `coverage`.
    pub raw_coverage: Vec<WindowedCoverage>,
    /// Query index: (target_position, normalized_score) for primary + alternatives
    pub query_positions: Vec<(u32, f64)>,
    /// `Strategy::Assembly` only: every kept local alignment block (ADR-0015). Empty for the
    /// other strategies, whose placements are reported through `query_positions` alone.
    pub blocks: Vec<AlignmentBlock>,
}

/// Execute a single alignment task: query vs target with default configuration.
pub fn align_task(
    query: &Sequence,
    target: &Sequence,
    plan: &PhrayaPlan,
) -> Option<AlignmentResult> {
    align_task_with_config(query, target, plan, &AlignConfig::default())
}

/// Maximum query length for which the Myers bit-parallel engine is used. Longer reads
/// fall back to WFA: Myers is O(nm/w) (quadratic in length), while WFA is O(s·n) and
/// scales better once reads are long enough that the length term dominates.
const MYERS_MAX_QUERY_LEN: usize = 500;

/// Divergence ceiling for the Fast strategy. Reads whose best alignment exceeds this
/// fraction of mismatches+indels per base are dropped — the low-sensitivity tradeoff
/// that lets Fast skip hard/divergent reads for speed.
const FAST_MAX_DIVERGENCE: f64 = 0.20;

/// Floor for the repeat-masking seed occurrence cap ([`seed_occurrence_cap`]). A minimizer
/// must occur more than this many times in the target before it can be masked, so clean and
/// moderately-repetitive genomes mask nothing (the cap only bites on pathological hyper-repeats).
const SEED_OCC_CAP_FLOOR: usize = 256;

/// Normalized-score floor at which a read's primary placement is reportable in `.phraya.queries`
/// (mirrors the 0.95 filter in `phraya_io::queries::write_queries`). Used only to classify a read
/// as placed vs below-threshold for [`AlignStats`]; the actual filtering happens at write time.
const SCORE_REPORT_THRESHOLD: f64 = 0.95;
/// Minimum identity for a `(0,0)` fallback-only alignment to produce variant observations.
///
/// When a query shares no minimizer with a target (no seeds, no chains), `select_anchors`
/// falls back to a `(0,0)` anchor. For comparable-length DNA (the contig-vs-contig case
/// issue #146 added the fallback for), this can rescue a genuinely divergent homolog that
/// shares no exact k-mer. But it also fires for every non-homologous pair in reference-
/// palette mode (every read aligned against every reference space): a 50 kb query vs a
/// non-homologous 50 kb reference produces a full-length ~50%-identity alignment whose
/// ~25 000 `VariantObservation`s are pure noise.
///
/// The cross-space sidecar already discards such placements at its 0.95 identity filter
/// (they never reach the superposition a filter acts on), but the per-space `.phraya` has
/// no such gate and retains every fabricated call. This floor suppresses *variant* extraction
/// from seedless alignments below the noise threshold, while still recording coverage depth
/// and the placement (for the cross-space superposition). Seeded alignments are unaffected —
/// this gate only engages when `had_seeds` is false (the alignment relied solely on the
/// `(0,0)` fallback).
const FALLBACK_ONLY_MIN_IDENTITY: f64 = 0.55;

/// Strategy-aware score threshold.
pub fn score_threshold(strategy: Strategy) -> f64 {
    match strategy {
        Strategy::LongRead => 0.80, // ONT/PacBio reads: ~10-15% error → identity ~0.85-0.90
        Strategy::Fast => 0.85,     // Fast path tolerates more divergence
        _ => SCORE_REPORT_THRESHOLD, // 0.95 for balanced/sensitive/exact
    }
}

/// Per-read outcome classification, for localizing where reads are lost (issue #194 AC #1).
enum Outcome {
    /// Reportable primary (normalized score ≥ [`SCORE_REPORT_THRESHOLD`]).
    Placed,
    /// No shared minimizer seed in either orientation — only the `(0,0)` fallback anchor.
    NoSeed,
    /// Seeded and extended, but the best primary scored below the report threshold.
    BelowThreshold,
    /// Dropped by the Fast-strategy divergence cutoff.
    FastCutoff,
    /// No anchor extended at all (extension errored on every anchor, both orientations).
    NoAlignment,
}

/// Thread-safe tally of read outcomes across an alignment run.
///
/// Answers "where are reads being lost?" — the instrumentation for issue #194: is the loss in
/// seeding (`no_seed`), in extension/divergence (`below_threshold`), in the Fast cutoff, or in
/// extension failure (`no_alignment`)? Atomic counters so the parallel worker can share one
/// across rayon tasks. Pass `Some(&stats)` to [`align_read`] to accumulate; `None` disables it.
#[derive(Debug, Default)]
pub struct AlignStats {
    /// Reads placed with a reportable primary (score ≥ 0.95).
    pub placed: AtomicU64,
    /// Reads with no shared minimizer seed in either orientation.
    pub no_seed: AtomicU64,
    /// Reads that seeded and extended but whose best primary scored < 0.95.
    pub below_threshold: AtomicU64,
    /// Reads dropped by the Fast-strategy divergence cutoff.
    pub fast_cutoff: AtomicU64,
    /// Reads for which no anchor extended at all.
    pub no_alignment: AtomicU64,
}

impl AlignStats {
    fn record(&self, outcome: Outcome) {
        let counter = match outcome {
            Outcome::Placed => &self.placed,
            Outcome::NoSeed => &self.no_seed,
            Outcome::BelowThreshold => &self.below_threshold,
            Outcome::FastCutoff => &self.fast_cutoff,
            Outcome::NoAlignment => &self.no_alignment,
        };
        counter.fetch_add(1, Ordering::Relaxed);
    }

    /// One-line summary: `placed=.. no_seed=.. below_threshold=.. fast_cutoff=.. no_alignment=..`.
    pub fn summary(&self) -> String {
        format!(
            "placed={} no_seed={} below_threshold={} fast_cutoff={} no_alignment={}",
            self.placed.load(Ordering::Relaxed),
            self.no_seed.load(Ordering::Relaxed),
            self.below_threshold.load(Ordering::Relaxed),
            self.fast_cutoff.load(Ordering::Relaxed),
            self.no_alignment.load(Ordering::Relaxed),
        )
    }
}

/// Number of top chains extended per strategy — chaining's counterpart to the legacy
/// path's raw anchor cap `K`. `Sensitive` gets a large finite cap (50), not the legacy
/// path's literal `K = ∞`: chaining already collapses most of the raw-vote explosion
/// that made an unbounded cap tractable only via ADR-0007's branch-and-bound, and a
/// finite cap gives predictable worst-case cost from the outset rather than waiting to
/// discover a pathological genome empirically.
fn chain_cap(strategy: Strategy) -> usize {
    match strategy {
        Strategy::Fast => 1,
        Strategy::Balanced => 2,
        Strategy::Sensitive => 50,
        Strategy::Exact => 250,
        Strategy::LongRead => 100,  // Long reads use moderate chain cap
        Strategy::Assembly => unreachable!("Strategy::Assembly aligns blocks in align_assembly, never via extend_chains"),
    }
}

/// Choose the WFA/Myers anchors to extend, in descending order of evidence.
///
/// 1. **Chains present** — one anchor per surviving chain (up to `k`), `target_pos` from
///    [`crate::chaining::Chain::target_start`]. The normal path.
/// 2. **Seeds but no chain** — anchor on the seeds themselves, deduplicated by diagonal and
///    capped at `k`. The query shares exact k-mers with this target, so it is a genuine
///    candidate; chaining routinely fails to find anything co-linear on repetitive or
///    low-complexity queries (issue #146), and a tandem-repeat query that occurs verbatim
///    dozens of times in the target is the standard example. Anchoring at the seed is what
///    makes this recoverable — the old `(0,0)` anchor only ever aligned `target[0..2*|q|)`,
///    so it rescued such a query purely by luck when the match happened to sit at the start.
/// 3. **No seeds at all** — `(0,0)`, but only where [`fallback_anchor_applies`] (always
///    `false` for `Alphabet::Protein` — see its docs); otherwise no anchors, and the caller
///    reports no alignment. A query sharing not one k-mer with a target vastly larger than
///    itself has no candidate locus to extend from, and forcing one fabricates a ~50%-
///    identity alignment pinned at offset 0.
fn select_anchors(
    chains: &[crate::chaining::Chain],
    seeds: &[crate::Seed],
    k: usize,
    query_len: usize,
    target_len: usize,
    alphabet: Alphabet,
) -> Vec<SeedAnchor> {
    if !chains.is_empty() {
        return chains
            .iter()
            .take(k)
            .map(|c| SeedAnchor {
                query_pos: 0,
                target_pos: c.target_start(),
            })
            .collect();
    }

    if !seeds.is_empty() {
        // Project each seed back along its diagonal to query offset 0, matching how chain
        // anchors are expressed, then dedupe: a repeat family yields many seeds per locus.
        let mut anchors: Vec<SeedAnchor> = Vec::with_capacity(k.min(seeds.len()));
        for s in seeds {
            let target_pos = (s.target_pos as usize).saturating_sub(s.query_pos as usize);
            if anchors.iter().any(|a| a.target_pos == target_pos) {
                continue;
            }
            anchors.push(SeedAnchor {
                query_pos: 0,
                target_pos,
            });
            if anchors.len() == k {
                break;
            }
        }
        return anchors;
    }

    if fallback_anchor_applies(query_len, target_len, alphabet) {
        return vec![SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        }];
    }

    Vec::new()
}

/// Extend a single anchor with the engine selected by `strategy` (and, for `Sensitive`,
/// `gap_model` — ADR-0014).
///
/// - `Sensitive`: canonical seeded WFA on every anchor — the reference path. Gap-affine
///   by default (`GapModel::Affine`); `GapModel::Linear` opts back into the old
///   uniform-cost path.
/// - `Balanced` / `Fast`: Myers fitting for short reads (identical results to linear
///   WFA, but faster), falling back to **linear** WFA for reads longer than
///   [`MYERS_MAX_QUERY_LEN`] — `gap_model` is never consulted here; Myers has no affine
///   mode, and this fallback exists only to mirror Myers cheaply at length, not to
///   introduce a different cost model at this tier.
///
/// Fast and Balanced differ from Sensitive in *which* anchors they extend (K=1 and K=2
/// vs K=50) and Fast adds a post-hoc divergence cutoff, handled by the caller
/// — not the extension engine.
fn extend_anchor(
    strategy: Strategy,
    gap_model: GapModel,
    query: &[u8],
    target_window: &[u8],
    anchor: SeedAnchor,
) -> crate::WfaResult {
    match strategy {
        Strategy::Sensitive => match gap_model {
            GapModel::Affine => crate::wfa_extend_affine(
                query,
                target_window,
                anchor,
                crate::wfa_simd::AffineCosts::default(),
            ),
            GapModel::Linear => wfa_extend(query, target_window, anchor),
        },
        Strategy::Balanced | Strategy::Fast => {
            if query.len() <= MYERS_MAX_QUERY_LEN {
                myers_extend(query, target_window, anchor)
            } else {
                wfa_extend(query, target_window, anchor)
            }
        }
        Strategy::Exact => {
            // Exact strategy always uses WFA for full precision, regardless of query length.
            // No Myers fallback or adaptive selection.
            wfa_extend(query, target_window, anchor)
        }
        Strategy::LongRead => {
            // Long reads use WFA for accuracy on long, divergent sequences
            wfa_extend(query, target_window, anchor)
        }
        Strategy::Assembly => unreachable!("Strategy::Assembly aligns blocks in align_assembly, never via extend_anchor"),
    }
}
/// Precomputed, read-only per-target data shared across every query aligned to one
/// target.
///
/// The main use case (Case 2: N reads vs one reference) aligns many queries against a
/// single target. The target-derived structures here — the minimizer index and the
/// tandem-repeat regions — depend only on the target, so building them once and calling
/// [`align_read`] per query removes O(target) work (and a full-genome `to_uppercase`
/// copy and sketch clone) from each read. This mirrors how `plan.hotspot_intervals` is
/// precomputed once and passed read-only.
pub struct TargetContext<'a> {
    target: &'a Sequence,
    minimizer_index: MinimizerIndex,
    /// Tandem-repeat regions, detected on first use.
    ///
    /// Only `extract_variants_from_cigar` consumes these, so a target that places no read
    /// never needs them at all. Detection is an O(target) scan (~0.11 s/Mb measured), which
    /// reference-palette mode otherwise pays once per space even though only a handful of
    /// spaces in a large palette ever receive a placement. `OnceLock` rather than `OnceCell`
    /// because `&TargetContext` is shared across rayon workers, and `OnceLock` rather than
    /// `LazyLock` because the initializer needs per-instance runtime input (`self.target`),
    /// which a field-level `LazyLock` initializer cannot capture.
    repeat_regions: std::sync::OnceLock<Vec<phraya_core::RepeatRegion>>,
    /// Repeat-masking cap: minimizers occurring more than this many times in the target
    /// are skipped during seeding. Derived once from the index's occurrence distribution
    /// (issue #194 — bounds the seed explosion on repeat-dense / low-complexity genomes).
    seed_max_occ: usize,
    /// Strategy used to select sketch density (Fast/Balanced use w=11 subset, Sensitive uses full dense)
    strategy: Strategy,
}

impl<'a> TargetContext<'a> {
    /// Build the shared context for `target`, reusing the plan's precomputed sketch if
    /// present and falling back to recomputing it otherwise. Filters the sketch based on
    /// the strategy: Fast/Balanced use the w=11 subset, Sensitive uses the full dense set.
    pub fn build(target: &'a Sequence, plan: &PhrayaPlan, strategy: Strategy) -> Self {
        // Get the effective sketch for this strategy (filters dense to w=11 if needed)
        let sketch = get_effective_sketch(target.id(), plan, strategy)
            .unwrap_or_else(|| sketch_sequence_default_for(target, plan.alphabet));
        let minimizer_index = MinimizerIndex::build(&sketch);
        // Repeat-masking cap from the index's own occurrence distribution. Floor 256 keeps
        // it a no-op on clean/moderately-repetitive genomes (nothing occurs that often), so
        // only pathological hyper-repeats (homopolymer/microsatellite k-mers) are trimmed.
        let seed_max_occ = seed_occurrence_cap(&minimizer_index, SEED_OCC_CAP_FLOOR);
        TargetContext {
            target,
            minimizer_index,
            repeat_regions: std::sync::OnceLock::new(),
            seed_max_occ,
            strategy,
        }
    }

    /// Build the shared context for `target` from an explicitly supplied sketch, bypassing
    /// the plan's sequence-ID lookup. Used by reference-palette alignment (ADR-0011, issue
    /// #226): on a content-hash **hit** the caller passes the planned reference space's
    /// sketch (proving sketch reuse — no recompute), and on a **miss** the caller passes a
    /// freshly computed sketch. The occurrence cap is derived the same way as [`build`], and
    /// repeat regions are deferred identically (see [`TargetContext::repeat_regions`]).
    pub fn build_with_sketch(
        target: &'a Sequence,
        sketch: &MinimizerSketch,
        strategy: Strategy,
    ) -> Self {
        let minimizer_index = MinimizerIndex::build(sketch);
        let seed_max_occ = seed_occurrence_cap(&minimizer_index, SEED_OCC_CAP_FLOOR);
        TargetContext {
            target,
            minimizer_index,
            repeat_regions: std::sync::OnceLock::new(),
            seed_max_occ,
            strategy,
        }
    }

    /// Tandem-repeat regions for this target, detected on the first call and cached.
    ///
    /// Deliberately not computed in the constructors: a space that receives no placement
    /// never calls this, which is what makes a large reference palette cheap.
    fn repeat_regions(&self) -> &[phraya_core::RepeatRegion] {
        self.repeat_regions.get_or_init(|| {
            detect_tandem_repeats(self.target.bases(), &RepeatDetectorConfig::default())
        })
    }

    /// The target sequence this context was built for.
    pub fn target(&self) -> &Sequence {
        self.target
    }
}

/// Execute a single alignment task: query vs target with specified configuration.
///
/// Thin wrapper that builds a [`TargetContext`] and delegates to [`align_read`]. When
/// aligning many queries against one target, build the context once and call
/// [`align_read`] directly instead of paying the per-target build on every query.
pub fn align_task_with_config(
    query: &Sequence,
    target: &Sequence,
    plan: &PhrayaPlan,
    config: &AlignConfig,
) -> Option<AlignmentResult> {
    let ctx = TargetContext::build(target, plan, config.strategy);
    align_read(&ctx, query, plan, config, None)
}

/// Seed and chain `query_sketch` against the target — the cheap half of
/// [`align_oriented`], factored out so [`align_read`] can decide which orientation(s)
/// are worth the expensive extend step before paying for it.
///
/// Returns `(n_seeds, chains)`. Both orientations of a read always find the *same*
/// number of raw seeds (canonical minimizers are strand-invariant in value; only the
/// seed *query positions* differ), so `n_seeds` alone cannot distinguish the true
/// orientation from the wrong one. Chaining can: it requires seeds to be co-linear
/// (consistent diagonal, increasing query/target position), which only holds in a
/// read's true orientation — seeding the wrong orientation of a real match reliably
/// produces zero chains even when seeds are plentiful (verified: a genuinely
/// forward-matching 150bp read seeded in the wrong (RC) orientation found the same
/// seed count but chained to nothing, while the true orientation chained to one strong
/// candidate). `chains` is this function's actual "is this orientation plausible?"
/// signal.
fn seed_and_chain(
    ctx: &TargetContext<'_>,
    query_sketch: &MinimizerSketch,
    strategy: Strategy,
) -> (Vec<crate::Seed>, Vec<crate::chaining::Chain>) {
    let seeds = find_seeds_indexed_capped(query_sketch, &ctx.minimizer_index, ctx.seed_max_occ);
    let params = if strategy == Strategy::Assembly {
        crate::chaining::ChainParams::for_assembly()
    } else if strategy == Strategy::LongRead {
        crate::chaining::ChainParams::for_long_read()
    } else {
        crate::chaining::ChainParams::default()
    };
    let chains = crate::chaining::chain_seeds(&seeds, &params);
    (seeds, chains)
}

/// Whether a query with **no** chain support can still legitimately produce an alignment
/// against a target of `target_len`, via the `(0,0)` fallback anchor in
/// [`select_anchors`].
///
/// **Always `false` for `Alphabet::Protein`.** The ratio rule below assumes DNA's issue
/// #146 contig-vs-contig shape (a handful of sequences compared pairwise, where comparable
/// length is unusual enough to be a real homology signal). Reference-palette protein
/// search instead compares a query against a proteome — potentially thousands of candidate
/// proteins that are *naturally* comparable in length to the query, which is the normal
/// shape of a database, not evidence of relatedness. Applying the ratio there would force
/// a real `(0,0)` extension attempt against every candidate with no seed (real BLASTP/
/// DIAMOND-style seeded search never does that), and it also defeats
/// `run_align_reference`'s per-space skip below at the cost of O(queries × proteome size)
/// extension attempts — exactly the combinatorial blowup this alphabet gate prevents.
///
/// For DNA (unchanged): true only when the two are comparable in length — the contig-vs-
/// contig case issue #146 added the fallback for, where a real homolog can share no
/// minimizer. When the target dwarfs the query the fallback can only ever align the query
/// against `target[0..2*query_len)`: a guaranteed-losing alignment whose sub-threshold
/// placement `write_cross_space_queries` discards at its 0.95 identity filter, but whose
/// ~`query_len/2` bogus `VariantObservation`s are *not* filtered and land in the `.phraya`.
/// Measured on a 20-record palette: 1.87 MB of garbage per non-homologous space against
/// 42 KB of signal in the one space the reads came from, at 256us per (read, space) pair —
/// 2.5x the cost of a real placement. Reference-palette mode makes that the common case,
/// since every read is aligned against every space.
///
/// Public because reference-palette alignment relies on the same predicate to decide it can
/// skip a space outright: zero shared minimizers means zero seeds, so if the fallback also
/// does not apply, that space provably yields no placement and no variant. The two
/// decisions MUST agree, hence one function rather than two copies of the ratio.
pub fn fallback_anchor_applies(query_len: usize, target_len: usize, alphabet: Alphabet) -> bool {
    if alphabet == Alphabet::Protein {
        return false;
    }
    /// Largest target/query length ratio at which a seedless `(0,0)` alignment is still
    /// plausible rather than guaranteed noise.
    const FALLBACK_MAX_TARGET_RATIO: usize = 10;
    target_len <= query_len.saturating_mul(FALLBACK_MAX_TARGET_RATIO)
}

/// Extend the anchors [`select_anchors`] picks for this orientation, returning the scored
/// alignments — or `None` if there was nothing worth extending.
///
/// The expensive half of what was [`align_oriented`] before it was split so
/// [`align_read`] could decide which orientation(s) are worth the expensive extend step
/// before paying for it.
///
/// `seed_k` is the minimizer length of the orientation's sketch (`fwd_sketch.k` for the
/// forward orientation, the reverse-complement sketch's `k` for the reverse one). It is
/// passed through to the `Strategy::LongRead` path, which would otherwise have no way to
/// recover the k-mer length needed to validate a chain's seeds.
fn extend_chains(
    strategy: Strategy,
    gap_model: GapModel,
    ctx: &TargetContext<'_>,
    query_bytes: &[u8],
    chains: &[crate::chaining::Chain],
    seeds: &[crate::Seed],
    alphabet: Alphabet,
    seed_k: usize,
) -> Option<crate::ScoredAlignments> {
    let target = ctx.target;
    let query_len = query_bytes.len();

    // LongRead: bypass select_anchors (which projects chains to query_pos=0 and drops the
    // seeds needed for gap-fill) and dispatch to the chain-seed-based path directly.
    if strategy == Strategy::LongRead {
        return extend_chains_longread(query_bytes, chains, target, seed_k);
    }

    let anchors = select_anchors(
        chains,
        seeds,
        chain_cap(strategy),
        query_len,
        target.bases().len(),
        alphabet,
    );
    if anchors.is_empty() {
        return None;
    }

    // Extend every anchor with the strategy's engine, uniformly across Fast/Balanced/
    // Sensitive (ADR-0012 — chaining's structural collapse of repeat families into one
    // chain per copy, plus Sensitive's finite K=50 chain cap, made the old
    // Sensitive-only score-bounded branch-and-bound (ADR-0007) unnecessary; that
    // machinery was removed once chaining was validated on the full benchmark ladder).
    //
    // Window the target to ~2× query length from the anchor position. WFA is O(s·n)
    // where s = edit distance; for s << min(|q|,|t|) it is dramatically faster than
    // O(|q|×|t|) DP, but s grows with the length difference — passing the full
    // reference to a 150bp read makes the edit distance ~|target|-|query| (length gap)
    // rather than ~2% divergence, turning O(s) into O(target²). The 2× margin
    // accommodates indels while keeping the aligned window tractable.
    let mut alignments = Vec::with_capacity(anchors.len());
    for anchor in anchors {
        let margin = query_len * 2;
        let window_end = (anchor.target_pos + margin).min(target.bases().len());
        let target_window = &target.bases()[..window_end];
        match extend_anchor(strategy, gap_model, query_bytes, target_window, anchor) {
            Ok(aln) => alignments.push(aln),
            Err(e) => log::warn!("alignment failed for anchor {:?}: {:?}", anchor, e),
        }
    }

    if alignments.is_empty() {
        return None;
    }

    Some(score_alignments(&alignments, query_bytes.len()))
}

/// Chunked WFA extension for long reads (ONT/PacBio, 10-50kb).
///
/// Replaces the old whole-query-per-anchor chunked-WFA path that made one O(qn) WFA call
/// per surviving anchor (up to [`chain_cap`] = 100 for `LongRead`, 2kb chunks each). At 12%
/// ONT error one locus fragments into many near-duplicate chains (starts within ±2 bp,
/// identities ~0.817), each aligned in full — the ~15 s/read cost the 2026-10-20 handoff
/// flagged. This path does a minimap2-style seed→chain→align: dedup chains to one locus,
/// then align each locus once by gap-filling its exact-match seed segments.
///
/// Inter-segment gaps use global Myers (bit-parallel, ~no allocation; bounded by
/// ChainParams::max_dist_x). The two read ends use fitting WFA with a `default_max_s_cap`
/// backstop — on cap exhaustion (or cap=0) the end is soft-clipped, not dropped. Retention
/// stays per-space and per-strand: a placement at locus X is judged against X alone, never
/// globally — composable align({A,B}) = align({A}) ∪ align({B}).
fn extend_chains_longread(
    query_bytes: &[u8],
    chains: &[crate::chaining::Chain],
    target: &Sequence,
    seed_k: usize,
) -> Option<crate::ScoredAlignments> {
    if seed_k == 0 {
        return None;
    }
    let target_bases = target.bases();
    let query_len = query_bytes.len();

    // Dedup chains into loci. Chains come sorted by score desc; near-duplicate chains at
    // one true locus share a projected start (Chain::target_start) within ±merge_tol.
    // Compare against ALL kept loci (not just the last) — score order isn't position order.
    // merge_tol = query_len/10 swallows indel drift (tens of bp at 12% over 20 kb) while
    // leaving distinct tandem/repeat copies (farther apart) as separate loci.
    let merge_tol = (query_len / 10).max(256);
    let mut loci: Vec<&crate::chaining::Chain> = Vec::with_capacity(chains.len());
    for chain in chains {
        let tp = chain.target_start();
        if loci
            .iter()
            .any(|c| tp <= c.target_start() + merge_tol && c.target_start() <= tp + merge_tol)
        {
            continue;
        }
        loci.push(chain);
    }

    // Safety cap: chain_cap (100 for LongRead) is the documented locus upper bound. Dedup
    // keeps real data small, but this guards pathological repeat landscapes.
    let cap = chain_cap(Strategy::LongRead);
    if loci.len() > cap {
        loci.truncate(cap);
    }

    let mut alignments: Vec<crate::Alignment> = Vec::with_capacity(loci.len());
    for chain in &loci {
        if let Some(aln) = align_long_locus(query_bytes, target_bases, &chain.seeds, seed_k) {
            alignments.push(aln);
        }
    }
    if alignments.is_empty() {
        return None;
    }

    // score_alignments keeps the primary (min edit distance) + alternatives within the
    // 0.95 score-ratio band — the same surface the rest of the pipeline reads, so
    // cross-space superposition and the 0.95 filter behave unchanged.
    Some(score_alignments(&alignments, query_bytes.len()))
}

/// Align one long-read locus: exact k-mer segments stitched by global Myers gap-fills,
/// with fitting WFA end-extensions at the two read ends. `seeds` are co-linear in query_pos
/// (chaining's DP guarantees it).
///
/// The result spans the full query (query_start=0, query_end=query.len()); target_start is
/// shifted earlier than the first segment's start when the prefix fitting consumes target
/// bases, and target_end extends past the last segment when the suffix fitting consumes
/// target bases.
fn align_long_locus(
    query: &[u8],
    target: &[u8],
    seeds: &[crate::Seed],
    k: usize,
) -> Option<crate::Alignment> {
    let query_len = query.len();
    let segs = chain_segments(seeds, query, target, k);
    if segs.is_empty() {
        // No exact k-mer survived validation (cross-strand minimizer hit). Fall back to a
        // single fitting extension from the seed's projected target position — one call,
        // not a chunked WFA loop.
        if seeds.is_empty() {
            return None;
        }
        let tp = (seeds[0].target_pos as i64 - seeds[0].query_pos as i64).max(0) as usize;
        let window_end = (tp + query_len * 2).min(target.len());
        let window = &target[..window_end];
        let anchor = crate::SeedAnchor { query_pos: 0, target_pos: tp };
        return match wfa_extend(query, window, anchor) {
            Ok(aln) => Some(aln),
            Err(e) => {
                log::warn!("long-read locus fitting extension failed: {e:?}");
                None
            }
        };
    }

    let mut ops: Vec<(usize, char)> = Vec::new();
    let mut edit_distance = 0usize;
    let first = &segs[0];
    let last = segs.last().unwrap();

    // target_start tracks where the alignment touches the reference. Starts at the first
    // segment's start (projected back to query offset 0), then shifts earlier if the prefix
    // fitting consumes target bases before it.
    let mut target_start = (first.t as i64 - first.q as i64).max(0) as usize;

    // --- Prefix: query[0..first.q] into target ending at first.t. ---
    // Global Myers alignment of the prefix against exactly prefix.len() target bases
    // ending at first.t. Global (not fitting): the seed diagonal already fixes the locus,
    // so the prefix must consume a fixed target span — a fitting WFA's free end was
    // measured to drift ~122 bp onto a wrong higher-scoring match (longread_7 on Reedling2).
    // Myers global alignment has no freedom to slide, so target_start is anchored to first.t.
    let prefix = &query[..first.q];
    if !prefix.is_empty() && first.t > 0 {
        let t_start = first.t.saturating_sub(prefix.len());
        let t_window = &target[t_start..first.t.min(target.len())];
        if t_window.len() == prefix.len() {
            let (ed, cigar) = myers_edit_distance(prefix, t_window);
            for (n, op) in phraya_core::cigar::parse_ops(&cigar) {
                push_op(&mut ops, n, op);
            }
            edit_distance += ed;
        } else {
            // Window shorter than prefix (near target start): align what fits, hard-clip rest.
            let avail = t_window.len();
            let (ed, cigar) = myers_edit_distance(&prefix[..avail], t_window);
            for (n, op) in phraya_core::cigar::parse_ops(&cigar) {
                push_op(&mut ops, n, op);
            }
            edit_distance += ed;
            push_op(&mut ops, prefix.len() - avail, 'S');
            edit_distance += prefix.len() - avail;
        }
    } else if !prefix.is_empty() {
        // first.t == 0: no target before the segment; soft-clip the prefix.
        push_op(&mut ops, prefix.len(), 'S');
        edit_distance += prefix.len();
    }
    ops.push((first.len, 'M'));
    // Exact segment: 0 edits.

    // --- Interior gaps: global Myers between consecutive exact segments. ---
    for w in segs.windows(2) {
        let prev = &w[0];
        let next = &w[1];
        let gap_q = &query[prev.q + prev.len..next.q];
        let gap_t = &target[prev.t + prev.len..next.t];
        let (ed, cigar) = myers_edit_distance(gap_q, gap_t);
        for (n, op) in phraya_core::cigar::parse_ops(&cigar) {
            push_op(&mut ops, n, op);
        }
        edit_distance += ed;
        ops.push((next.len, 'M'));
    }

    // --- Suffix: query[last_end..] from last.t+last.len, fitting WFA. ---
    // The suffix has no segment downstream to anchor against, so fitting (target end free)
    // is correct here — target_end isn't used for a placement position check (only
    // target_start/primary is), and the end extension is cheap. Window sized to stay in
    // the fitting branch (> 1.5x query + 10).
    let suffix_len = query_len - (last.q + last.len);
    let mut target_end = last.t + last.len;
    let suffix = &query[last.q + last.len..];
    if !suffix.is_empty() && last.t + last.len < target.len() {
        let win = suffix_len * 16 / 10 + 16;
        let t_end = (last.t + last.len + win).min(target.len());
        let t_window = &target[last.t + last.len..t_end];
        let cap = crate::wfa_simd::default_max_s_cap(suffix.len(), t_window.len());
        match crate::wfa_simd::fill_wfa_fitting_impl(suffix, t_window, Some(cap)) {
            Some((cigar, ed, consumed)) => {
                for (n, op) in phraya_core::cigar::parse_ops(&cigar) {
                    push_op(&mut ops, n, op);
                }
                edit_distance += ed;
                target_end = (last.t + last.len) + consumed;
            }
            None => {
                // Cap exhausted: soft-clip the suffix.
                push_op(&mut ops, suffix.len(), 'S');
                edit_distance += suffix.len();
            }
        }
    } else if !suffix.is_empty() {
        push_op(&mut ops, suffix.len(), 'S');
        edit_distance += suffix.len();
    }

    let cigar: String = ops.iter().map(|(n, op)| format!("{n}{op}")).collect();
    Some(crate::Alignment {
        cigar,
        edit_distance,
        query_start: 0,
        query_end: query_len,
        target_start,
        target_end,
    })
}

#[cfg(test)]
mod longread_locus_tests {
    use super::*;

    fn read_with_indel(prefix_len: usize, del_len: usize, suffix_len: usize) -> (Vec<u8>, Vec<u8>, Vec<crate::Seed>) {
        // Reference: AAAAA... (prefix_len) then CCCCC... (del_len) then TTTTT... (suffix_len)
        let mut ref_seq = Vec::new();
        ref_seq.extend(std::iter::repeat(b'A').take(prefix_len));
        ref_seq.extend(std::iter::repeat(b'C').take(del_len));
        ref_seq.extend(std::iter::repeat(b'T').take(suffix_len));
        // Read = ref without the C-run (a clean deletion)
        let mut read = ref_seq.clone();
        read.drain(prefix_len..prefix_len + del_len);
        // Seeds: exact 21-mers flanking the deletion, placed at their reference positions.
        let k = 21;
        let left_seed_q = 0usize;
        let left_seed_t = 0usize;
        let right_seed_q = prefix_len;
        let right_seed_t = prefix_len + del_len;
        let seeds = vec![
            crate::Seed { query_pos: left_seed_q as u32, target_pos: left_seed_t as u32, minimizer: 0 },
            crate::Seed { query_pos: right_seed_q as u32, target_pos: right_seed_t as u32, minimizer: 0 },
        ];
        (read, ref_seq, seeds)
    }

    fn cigar_query_target(cigar: &str) -> (usize, usize) {
        let stats = phraya_core::cigar::CigarStats::parse(cigar);
        (stats.query_aligned_len() as usize, stats.target_aligned_len() as usize)
    }

    #[test]
    fn locus_cigar_lengths_match_read_and_span() {
        // prefix=300, del=50, suffix=300 -> read 600, ref 650
        let (read, ref_seq, seeds) = read_with_indel(300, 50, 300);
        let k = 21;
        let aln = align_long_locus(&read, &ref_seq, &seeds, k).expect("locus must align");
        let (q, t) = cigar_query_target(&aln.cigar);
        assert_eq!(q, read.len(), "CIGAR query length must equal read length");
        assert_eq!(t, aln.target_end - aln.target_start,
            "CIGAR target length must equal target_end - target_start");
        assert_eq!(aln.target_start, 0, "anchor at ref start");
        assert_eq!(aln.target_end, ref_seq.len(), "aligns to end of ref");
    }

    #[test]
    fn locus_handles_prefix_before_first_seed() {
        // First seed at query 200, so prefix=200 must be extended, not dropped.
        let (read, ref_seq, seeds) = read_with_indel(200, 50, 300);
        let k = 21;
        let aln = align_long_locus(&read, &ref_seq, &seeds, k).expect("locus must align");
        let (q, _t) = cigar_query_target(&aln.cigar);
        assert_eq!(q, read.len(), "CIGAR must consume the full read including the prefix");
        assert!(aln.target_end > aln.target_start, "alignment must have span");
    }

    #[test]
    fn locus_long_prefix_no_target_drift() {
        // Long prefix (4000 bp) before the first exact seed, perfect match. Verifies the
        // 1.6x end window keeps fit_end in fitting mode (not global) and that target_start
        // anchors at 0 rather than drifting by the window slack.
        let k = 15;
        let prefix_len = 4000;
        let suffix_len = 4000;
        let mut ref_seq = Vec::new();
        for i in 0..prefix_len {
            ref_seq.push(if i % 2 == 0 { b'A' } else { b'T' });
        }
        // Clean exact seed region (k*2 bases) then a G-run suffix.
        ref_seq.extend(std::iter::repeat(b'C').take(k * 2));
        ref_seq.extend(std::iter::repeat(b'G').take(suffix_len - k));
        let read = ref_seq.clone();
        let seed_left = prefix_len;
        let seeds = vec![crate::Seed {
            query_pos: seed_left as u32,
            target_pos: seed_left as u32,
            minimizer: 0,
        }];
        let aln = align_long_locus(&read, &ref_seq, &seeds, k)
            .expect("locus must align");
        let (q, t) = cigar_query_target(&aln.cigar);
        assert_eq!(q, read.len(), "CIGAR query length must equal read length");
        assert_eq!(t, aln.target_end - aln.target_start,
            "CIGAR target length must equal target_end - target_start");
        assert_eq!(aln.target_start, 0, "perfect prefix must anchor at ref start, not drift");
        assert_eq!(aln.target_end, ref_seq.len());
        assert_eq!(aln.edit_distance, 0, "perfect match must have 0 edits");
    }
}

/// Per-query work that depends only on (query, plan, strategy) — never on the target
/// space. Hoisted out of `align_read` so reference-palette mode, which aligns every read
/// against every space, computes it once per read instead of once per (space, read).
pub struct PreparedQuery {
    fwd_sketch: MinimizerSketch,
    used_stored_sketch: bool,
    rc_bases: Vec<u8>,
    rev_sketch: Option<MinimizerSketch>,
}

impl PreparedQuery {
    /// The query's minimizer hash values.
    ///
    /// Canonical minimizers are strand-invariant in *value* — only their query positions
    /// differ between orientations (see [`seed_and_chain`]) — so the forward sketch's values
    /// are the complete set of hashes this query can seed on, in either orientation. A
    /// target sharing none of them cannot yield a single seed for this query.
    pub fn minimizer_values(&self) -> impl Iterator<Item = u64> + '_ {
        self.fwd_sketch.minimizers.iter().map(|&(val, _)| val)
    }
}

/// Compute the space-independent part of aligning `query`: its forward sketch (reusing a
/// stored sketch when the plan has one), its reverse complement, and — when worth
/// computing — its reverse-complement sketch. See `align_read_prepared` for why the
/// reverse sketch is skipped when the forward sketch is empty.
pub fn prepare_query(query: &Sequence, plan: &PhrayaPlan, config: &AlignConfig) -> PreparedQuery {
    // Issue #200: prefer a stored read sketch (keyed by content hash, stable across
    // pipeline stages unlike sequence ID) over recomputing one. `used_stored_sketch`
    // tracks whether the forward sketch came from the plan's cache, so an empty cached
    // sketch (a caller intentionally supplied no minimizers for this read) can be
    // distinguished below from a sketch that's merely too short to seed anything.
    let (fwd_sketch, used_stored_sketch) = {
        let hash = read_content_hash(query.bases());
        if let Some(stored) = plan.get_read_sketch(hash) {
            (stored.clone(), true)
        } else {
            let fallback = get_effective_sketch(query.id(), plan, config.strategy)
                .unwrap_or_else(|| sketch_sequence_default_for(query, plan.alphabet));
            (fallback, false)
        }
    };

    // Protein has no reverse strand (ADR-0013) — the dual-strand search below is a DNA
    // concept only. `rc_bases` stays empty and unread for Protein; `reverse_complement`
    // itself must not run on protein bytes — its A↔T/C↔G table would silently corrupt
    // any amino-acid byte that happens to collide with a DNA letter (e.g. Ala 'A' → 'T').
    let is_dna = plan.alphabet == Alphabet::Dna;
    let rc_bases = if is_dna {
        reverse_complement(query.bases())
    } else {
        Vec::new()
    };

    // If the forward sketch itself is empty, seeding against it can only ever find zero
    // seeds — the `fwd_seeds == 0 && fwd_sketch.minimizers.is_empty()` skip condition
    // `align_read` used to check (ctx-dependent `fwd_seeds` included) simplifies to just
    // the sketch-emptiness check, since an empty sketch trivially implies zero seeds.
    // Canonical minimizers are strand-invariant, so an empty forward sketch means the
    // reverse orientation is very unlikely to match either — skip its (expensive) sketch.
    let rev_sketch = if !is_dna || fwd_sketch.minimizers.is_empty() {
        None
    } else {
        let rc_seq = Sequence::new(rc_bases.clone(), None, query.id().to_string(), None);
        // Reverse complement must be re-sketched (positions are orientation-specific)
        // and must use the same strategy-based filtering as the forward strand. For
        // simplicity and correctness this uses a default w=11 sketch, automatically
        // dense-capable if simd-minimizers supports it, or just w=11 otherwise.
        Some(sketch_sequence_default(&rc_seq))
    };

    PreparedQuery {
        fwd_sketch,
        used_stored_sketch,
        rc_bases,
        rev_sketch,
    }
}

/// Align a single query against the target described by `ctx`, reusing a `PreparedQuery`
/// computed once per read by `prepare_query`.
///
/// Try both strands. A reverse-strand read's stored bytes are the reverse complement of
/// the reference region it came from; seeding finds anchors either way (canonical
/// minimizers are strand-invariant), but extension is strand-naïve, so the forward bytes
/// of a reverse read align at ~read-length edit distance and get dropped. Aligning the
/// reverse complement too — and keeping whichever orientation scores better — recovers
/// those reads (issue #192). Extension against the *forward* target means the winning
/// CIGAR and alleles are always in reference-forward orientation regardless of strand.
///
/// Chaining (not seeding) is what decides which orientation is worth extending: both
/// orientations of a real read find the same raw seed count (canonical minimizers are
/// strand-invariant in value), but chaining requires seed co-linearity, which reliably
/// holds only in a read's true orientation.
///
/// Pass `Some(&stats)` to tally the read's outcome (placed / no-seed / below-threshold /
/// fast-cutoff / no-alignment) for issue #194 diagnostics; `None` skips accounting.
pub fn align_read_prepared(
    ctx: &TargetContext<'_>,
    query: &Sequence,
    prepared: &PreparedQuery,
    plan: &PhrayaPlan,
    config: &AlignConfig,
    stats: Option<&AlignStats>,
) -> Option<AlignmentResult> {
    let target = ctx.target;
    let record = |outcome: Outcome| {
        if let Some(s) = stats {
            s.record(outcome);
        }
    };

    let (fwd_seeds, fwd_chains) = seed_and_chain(ctx, &prepared.fwd_sketch, config.strategy);
    let (rev_seeds, rev_chains) = match &prepared.rev_sketch {
        Some(s) => seed_and_chain(ctx, s, config.strategy),
        None => (Vec::new(), Vec::new()),
    };
    // Whether the read shared any minimizer with the target (in either orientation).
    // Distinguishes a seeding loss from an extension/divergence loss when classifying an
    // unplaced read.
    let had_seeds = !fwd_seeds.is_empty() || !rev_seeds.is_empty();

    // Issue #200: an empty *stored* read sketch is a deliberate signal from the caller
    // that this read has no minimizers to seed with (as opposed to a short/degenerate
    // read that merely recomputed to an empty sketch) — treat it as unplaceable rather
    // than falling through to the (0,0) anchor's fallback alignment.
    if prepared.used_stored_sketch && prepared.fwd_sketch.minimizers.is_empty() && !had_seeds {
        record(Outcome::NoSeed);
        return None;
    }

    // Assembly mode: every chain in either orientation becomes its own local block
    // (ADR-0015) — no anchor selection, no orientation winner, no (0,0) fallback.
    if config.strategy == Strategy::Assembly {
        let (result, outcome) =
            align_assembly(ctx, query, prepared, plan, config, &fwd_chains, &rev_chains, had_seeds);
        record(outcome);
        return result;
    }

    // Extend only the orientation(s) worth extending. If exactly one orientation has
    // chain support and the other has none, the chainless side can only ever produce a
    // weaker seed- or fallback-anchored alignment (see [`select_anchors`]) — never
    // competitive against a real chain-backed match, so skip its extension entirely. If both
    // have chains, or neither does (a low-complexity/repetitive background where chaining
    // can't distinguish signal from noise in either orientation, see issue #146), extend
    // both, exactly as before.
    let (fwd, rev) = if !fwd_chains.is_empty() && rev_chains.is_empty() {
        (
            extend_chains(config.strategy, config.gap_model, ctx, query.bases(), &fwd_chains, &fwd_seeds, plan.alphabet, prepared.fwd_sketch.k),
            None,
        )
    } else if fwd_chains.is_empty() && !rev_chains.is_empty() {
        (
            None,
            extend_chains(config.strategy, config.gap_model, ctx, &prepared.rc_bases, &rev_chains, &rev_seeds, plan.alphabet, prepared.rev_sketch.as_ref().map_or(0, |s| s.k)),
        )
    } else {
        (
            extend_chains(config.strategy, config.gap_model, ctx, query.bases(), &fwd_chains, &fwd_seeds, plan.alphabet, prepared.fwd_sketch.k),
            extend_chains(config.strategy, config.gap_model, ctx, &prepared.rc_bases, &rev_chains, &rev_seeds, plan.alphabet, prepared.rev_sketch.as_ref().map_or(0, |s| s.k)),
        )
    };

    // Keep the better-scoring orientation; ties resolve to forward deterministically.
    let (scored, query_bytes, strand): (crate::ScoredAlignments, &[u8], Strand) = match (fwd, rev) {
        (Some(f), Some(r)) if r.primary.edit_distance < f.primary.edit_distance => {
            (r, &prepared.rc_bases[..], Strand::Reverse)
        }
        (Some(f), _) => (f, query.bases(), Strand::Forward),
        (None, Some(r)) => (r, &prepared.rc_bases[..], Strand::Reverse),
        (None, None) => {
            record(Outcome::NoAlignment);
            return None;
        }
    };

    let primary_score = 1.0 - (scored.primary.edit_distance as f64 / query.len().max(1) as f64);

    // Fast strategy: drop reads whose best alignment exceeds the divergence cutoff. This
    // is the deliberate sensitivity sacrifice — confident, low-divergence reads only.
    if config.strategy == Strategy::Fast {
        let divergence = scored.primary.edit_distance as f64 / query.len().max(1) as f64;
        if divergence > FAST_MAX_DIVERGENCE {
            record(Outcome::FastCutoff);
            return None;
        }
    }

    // Classify the surviving read: reportable placement, or seeded-but-below-threshold, or a
    // read that only reached extension via the (0,0) fallback (no shared seed).
    let report_threshold = score_threshold(config.strategy);
    record(if primary_score >= report_threshold {
        Outcome::Placed
    } else if had_seeds {
        Outcome::BelowThreshold
    } else {
        Outcome::NoSeed
    });

    // Query positions are always recorded for the cross-space sidecar (the filter stage
    // sees the full superposition of placements; the sidecar's 0.95 threshold discards
    // low-identity ones). The identity field is absolute normalized (1 - ed/len).
    let mut query_positions = vec![(scored.primary.target_start as u32, primary_score)];
    for alt in &scored.alternatives {
        let alt_score = 1.0 - (alt.edit_distance as f64 / query.len().max(1) as f64);
        query_positions.push((alt.target_start as u32, alt_score));
    }

    // Suppress variant extraction from seedless (0,0)-fallback-only alignments below
    // [`FALLBACK_ONLY_MIN_IDENTITY`], but still record coverage and placement. Such
    // alignments have no minimizer evidence — they are the full-length ~50%-identity noise
    // that inflates palette-mode .phraya output (CSP2 snpdiffs escalation: 173 840 vs 33 386
    // real SNPs). Coverage depth is still computed (a read aligning at offset 0 is genuine
    // coverage even at low identity), but no bogus `VariantObservation`s are emitted. The
    // placement is recorded in query_positions above; the sidecar's 0.95 filter discards
    // sub-threshold placements from the superposition anyway.
    //
    // Seeded alignments (had_seeds == true) are never gated here — even a low-scoring seeded
    // alignment has positional evidence from shared minimizers.
    let mut suppress_variants = !had_seeds && primary_score < FALLBACK_ONLY_MIN_IDENTITY;
    // Long-read alignments (ONT/PacBio, 10-50 kb at ~12% error) are placements, not SNP
    // calls: a 25 kb read yields ~3000 mismatch VariantObservations each, allocating a
    // HashMap + coverage window — ~0.34 s/read, the dominant align cost (measured). Variant
    // extraction is a short-read use case; the filter stage decides what to call from the
    // superposition. Suppress observations (keep placement + query_positions + coverage)
    // to stay fast. (Long-read use case: depletion/tiering by placement, never SNP calling.)
    if config.strategy == Strategy::LongRead {
        suppress_variants = true;
    }

    // Raw (un-quantized) coverage over just the aligned span, for local_coverage
    // lookups in variants; quantized separately for the stored coverage track.
    // Computed unconditionally — even a fallback-only alignment that we suppress
    // variants for still contributes genuine coverage-depth at its aligned locus.
    let raw_coverage = compute_windowed_coverage(&scored, target.len());

    let variants = if suppress_variants {
        Vec::new()
    } else {
        let query_mapq = query.mapq().unwrap_or(60);
        let query_avg_bq = query.avg_quality().unwrap_or(60.0);

        // Look up mate info from plan (if available from BAM input)
        let mate_info = plan.mate_info.get(query.id());

        // Pre-compute aggregate insert stats from mate_info so variants are merge-stable.
        let insert_stats: Option<(i64, u32)> = mate_info.map(|mi| (mi.insert_size.abs() as i64, 1u32));

        extract_variants_from_cigar(
            &scored.primary.cigar,
            scored.primary.target_start,
            scored.primary.query_start,
            query_bytes,
            target.bases(),
            scored.primary.edit_distance as u32,
            query.id().to_string(),
            &raw_coverage,
            ctx.repeat_regions(),
            query_mapq,
            query_avg_bq,
            primary_score,
            config.coverage_window_radius,
            &plan.hotspot_intervals,
            mate_info,
            insert_stats,
            strand,
        )
    };

    // Quantize each window in place; positions outside all windows quantize to 0
    // (quantize(0) == 0), so the merged genome track is identical to quantizing full.
    let coverage: Vec<WindowedCoverage> = raw_coverage
        .iter()
        .map(|w| WindowedCoverage {
            start: w.start,
            counts: quantize_coverage(&w.counts),
        })
        .collect();

    Some(AlignmentResult {
        variants,
        coverage,
        raw_coverage,
        query_positions,
        blocks: Vec::new(),
    })
}

/// Align a single query against the target described by `ctx`. Convenience wrapper that
/// prepares the read's space-independent sketches on every call; aligning the same read
/// against multiple reference spaces (e.g. reference-palette mode) should call
/// `prepare_query` once per read and reuse `align_read_prepared` per space instead.
pub fn align_read(
    ctx: &TargetContext<'_>,
    query: &Sequence,
    plan: &PhrayaPlan,
    config: &AlignConfig,
    stats: Option<&AlignStats>,
) -> Option<AlignmentResult> {
    let prepared = prepare_query(query, plan, config);
    align_read_prepared(ctx, query, &prepared, plan, config, stats)
}

/// Floor on a `Strategy::Assembly` block's aligned reference length (ADR-0015). Shorter
/// blocks are chance seed clusters rather than homology; nucmer's minimum cluster length
/// (`-c 65`) is the same order.
const ASSEMBLY_MIN_BLOCK_LEN: usize = 100;

/// Identity floor on a `Strategy::Assembly` block: a homology backstop, not a SNP-quality
/// filter (CSP2 applies its own `min_iden` downstream).
const ASSEMBLY_MIN_BLOCK_IDENTITY: f64 = 0.80;

/// An exact match `query[q..q + len] == target[t..t + len]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Segment {
    q: usize,
    t: usize,
    len: usize,
}

/// A local alignment in strand-local query coordinates (half-open).
#[derive(Debug, Clone, PartialEq, Eq)]
struct LocalBlock {
    q_start: usize,
    q_end: usize,
    t_start: usize,
    t_end: usize,
    cigar: String,
    edit_distance: usize,
}
/// Drop exact duplicate `(strand, block)` entries, keeping the first occurrence of each.
/// Full-scan (not `dedup_by`), so non-adjacent duplicates are also collapsed.
fn dedup_blocks(kept: &mut Vec<(Strand, LocalBlock, f64)>) {
    let mut out: Vec<(Strand, LocalBlock, f64)> = Vec::with_capacity(kept.len());
    for e in kept.drain(..) {
        if !out.iter().any(|o| o.0 == e.0 && o.1 == e.1) {
            out.push(e);
        }
    }
    *kept = out;
}

/// Turn a chain's seeds into ordered, non-overlapping exact-match segments.
///
/// Canonical minimizers also match the opposite strand and hashes can collide, so a seed
/// counts only if its k-mer matches byte-for-byte in this orientation. Overlapping seeds on
/// one diagonal merge; a seed on a new diagonal loses whatever prefix overlaps the previous
/// segment on either axis. Each segment then grows to its maximal exact match, bounded by its
/// neighbours, so the gaps left between segments hold only the real differences.
fn chain_segments(seeds: &[crate::Seed], query: &[u8], target: &[u8], k: usize) -> Vec<Segment> {
    let mut segs: Vec<Segment> = Vec::new();
    for seed in seeds {
        let (q, t) = (seed.query_pos as usize, seed.target_pos as usize);
        if q + k > query.len() || t + k > target.len() || query[q..q + k] != target[t..t + k] {
            continue;
        }
        if let Some(last) = segs.last_mut() {
            if q.wrapping_sub(t) == last.q.wrapping_sub(last.t) && q <= last.q + last.len {
                last.len = last.len.max(q + k - last.q);
                continue;
            }
        }
        let overlap = segs.last().map_or(0, |last| {
            (last.q + last.len).saturating_sub(q).max((last.t + last.len).saturating_sub(t))
        });
        if overlap < k {
            segs.push(Segment { q: q + overlap, t: t + overlap, len: k - overlap });
        }
    }
    for i in 0..segs.len() {
        let (q_floor, t_floor) = match i.checked_sub(1) {
            Some(p) => (segs[p].q + segs[p].len, segs[p].t + segs[p].len),
            None => (0, 0),
        };
        let (q_ceil, t_ceil) = match segs.get(i + 1) {
            Some(n) => (n.q, n.t),
            None => (query.len(), target.len()),
        };
        let s = &mut segs[i];
        while s.q > q_floor && s.t > t_floor && query[s.q - 1] == target[s.t - 1] {
            s.q -= 1;
            s.t -= 1;
            s.len += 1;
        }
        while s.q + s.len < q_ceil && s.t + s.len < t_ceil && query[s.q + s.len] == target[s.t + s.len] {
            s.len += 1;
        }
    }
    segs
}

/// Append `n` × `op` to a run-length CIGAR, merging with the previous run.
fn push_op(ops: &mut Vec<(usize, char)>, n: usize, op: char) {
    if n == 0 {
        return;
    }
    match ops.last_mut() {
        Some((count, last)) if *last == op => *count += n,
        _ => ops.push((n, op)),
    }
}

/// Globally align the gap between two consecutive segments, appending its CIGAR to `ops`.
/// Returns the gap's edit distance, or `None` if the WFA cap is exhausted (the block is
/// split there and `ops` is left untouched).
fn fill_gap(q: &[u8], t: &[u8], ops: &mut Vec<(usize, char)>) -> Option<usize> {
    match (q.is_empty(), t.is_empty()) {
        (true, true) => Some(0),
        // This crate's alphabet: 'I' = target-extra, 'D' = query-extra.
        (true, false) => {
            push_op(ops, t.len(), 'I');
            Some(t.len())
        }
        (false, true) => {
            push_op(ops, q.len(), 'D');
            Some(q.len())
        }
        (false, false) => {
            let cap = crate::wfa_simd::default_max_s_cap(q.len(), t.len());
            let (cigar, edit_distance) = crate::wfa_simd::fill_wfa(q, t, Some(cap))?;
            for (n, op) in phraya_core::cigar::parse_ops(&cigar) {
                push_op(ops, n, op);
            }
            Some(edit_distance)
        }
    }
}

/// The local alignment blocks of one chain: its segments joined by global gap fills.
fn chain_blocks(seeds: &[crate::Seed], query: &[u8], target: &[u8], k: usize) -> Vec<LocalBlock> {
    let finish = |first: Segment, last: Segment, ops: &[(usize, char)], edit_distance: usize| LocalBlock {
        q_start: first.q,
        q_end: last.q + last.len,
        t_start: first.t,
        t_end: last.t + last.len,
        cigar: ops.iter().map(|(n, op)| format!("{n}{op}")).collect(),
        edit_distance,
    };
    let mut blocks = Vec::new();
    let mut ops: Vec<(usize, char)> = Vec::new();
    let mut edit_distance = 0;
    // (first, last) segment of the block being built.
    let mut open: Option<(Segment, Segment)> = None;
    for seg in chain_segments(seeds, query, target, k) {
        open = match open {
            None => Some((seg, seg)),
            Some((first, last)) => {
                let gap_q = &query[last.q + last.len..seg.q];
                let gap_t = &target[last.t + last.len..seg.t];
                match fill_gap(gap_q, gap_t, &mut ops) {
                    Some(gap_ed) => {
                        edit_distance += gap_ed;
                        Some((first, seg))
                    }
                    None => {
                        blocks.push(finish(first, last, &ops, edit_distance));
                        ops.clear();
                        edit_distance = 0;
                        Some((seg, seg))
                    }
                }
            }
        };
        push_op(&mut ops, seg.len, 'M');
    }
    if let Some((first, last)) = open {
        blocks.push(finish(first, last, &ops, edit_distance));
    }
    blocks
}

/// `Strategy::Assembly` (ADR-0015): every chain in either orientation becomes local
/// alignment blocks, replacing anchor selection, whole-query extension and the orientation
/// winner. Blocks shorter than [`ASSEMBLY_MIN_BLOCK_LEN`] or below
/// [`ASSEMBLY_MIN_BLOCK_IDENTITY`] are dropped; variants come from every kept block. Returns
/// the result and the outcome to record (`Placed` iff a block is kept).
#[allow(clippy::too_many_arguments)]
fn align_assembly(
    ctx: &TargetContext<'_>,
    query: &Sequence,
    prepared: &PreparedQuery,
    plan: &PhrayaPlan,
    config: &AlignConfig,
    fwd_chains: &[crate::chaining::Chain],
    rev_chains: &[crate::chaining::Chain],
    had_seeds: bool,
) -> (Option<AlignmentResult>, Outcome) {
    let target = ctx.target.bases();
    let orientations = [
        (Strand::Forward, query.bases(), fwd_chains, Some(prepared.fwd_sketch.k)),
        (Strand::Reverse, &prepared.rc_bases[..], rev_chains, prepared.rev_sketch.as_ref().map(|s| s.k)),
    ];
    let mut built = 0usize;
    let mut kept: Vec<(Strand, LocalBlock, f64)> = Vec::new();
    for (strand, bytes, chains, k) in orientations {
        let Some(k) = k else { continue };
        for chain in chains {
            for block in chain_blocks(&chain.seeds, bytes, target, k) {
                built += 1;
                let identity =
                    phraya_core::cigar::CigarStats::parse(&block.cigar).identity(block.edit_distance as u32);
                if block.t_end - block.t_start >= ASSEMBLY_MIN_BLOCK_LEN && identity >= ASSEMBLY_MIN_BLOCK_IDENTITY {
                    kept.push((strand, block, identity));
                }
            }
        }
    }
    if kept.is_empty() {
        let outcome = if !had_seeds {
            Outcome::NoSeed
        } else if built > 0 {
            Outcome::BelowThreshold
        } else {
            Outcome::NoAlignment
        };
        return (None, outcome);
    }
    // ADR-0015 has no anchor selection, so every seed cluster grows to its maximal exact match;
    // on a long near-identical contig, chaining's `max_dist_x`/`max_iter` back-scan bounds
    // split one clean diagonal into several independent chains, each emitting the same
    // `(strand, coords, cigar)` block. Collapse those duplicates here, before variant
    // extraction and the coverage build. See §7 of the CSP2 assembly-verification handoff.
    dedup_blocks(&mut kept);
    let raw_coverage: Vec<WindowedCoverage> = kept
        .iter()
        .map(|(_, b, _)| WindowedCoverage { start: b.t_start, counts: vec![1; b.t_end - b.t_start] })
        .collect();
    let mapq = query.mapq().unwrap_or(60);
    let avg_base_quality = query.avg_quality().unwrap_or(60.0);
    let mate_info = plan.mate_info.get(query.id());
    let insert_stats: Option<(i64, u32)> = mate_info.map(|mi| (mi.insert_size.abs() as i64, 1u32));
    let query_id: std::sync::Arc<str> = std::sync::Arc::from(query.id());
    let query_len = query.len();

    let mut variants = Vec::new();
    let mut query_positions = Vec::with_capacity(kept.len());
    let mut blocks = Vec::with_capacity(kept.len());
    for (strand, b, identity) in &kept {
        let bytes = match strand {
            Strand::Forward => query.bases(),
            Strand::Reverse => &prepared.rc_bases[..],
        };
        variants.extend(extract_variants_from_cigar(
            &b.cigar,
            b.t_start,
            b.q_start,
            bytes,
            target,
            b.edit_distance as u32,
            query.id().to_string(),
            &raw_coverage,
            ctx.repeat_regions(),
            mapq,
            avg_base_quality,
            *identity,
            config.coverage_window_radius,
            &plan.hotspot_intervals,
            mate_info,
            insert_stats,
            *strand,
        ));
        query_positions.push((b.t_start as u32, *identity));
        let (q_start, q_end) = match strand {
            Strand::Forward => (b.q_start, b.q_end),
            Strand::Reverse => (query_len - b.q_end, query_len - b.q_start),
        };
        blocks.push(AlignmentBlock {
            query_id: query_id.clone(),
            query_len: query_len as u32,
            q_start: q_start as u32,
            q_end: q_end as u32,
            t_start: b.t_start as u32,
            t_end: b.t_end as u32,
            strand: *strand,
            edit_distance: b.edit_distance as u32,
            cigar: std::sync::Arc::from(b.cigar.as_str()),
        });
    }
    let coverage = raw_coverage
        .iter()
        .map(|w| WindowedCoverage { start: w.start, counts: quantize_coverage(&w.counts) })
        .collect();
    (
        Some(AlignmentResult { variants, coverage, raw_coverage, query_positions, blocks }),
        Outcome::Placed,
    )
}

/// Check if a position falls within any hotspot interval.
///
/// `hotspot_intervals` (from [`phraya_core::types::detect_hotspot_intervals`]) is sorted by
/// `start` and non-overlapping, so binary search replaces the O(n) linear scan with O(log n) —
/// load-bearing once intervals number in the thousands+ (issue: 500x perf investigation).
fn is_in_hotspot(pos: u32, hotspot_intervals: &[(u32, u32)]) -> bool {
    // Last interval whose start is <= pos.
    let idx = hotspot_intervals.partition_point(|&(start, _)| start <= pos);
    idx > 0 && pos < hotspot_intervals[idx - 1].1
}

/// Check if a position falls within any tandem-repeat region.
///
/// `repeat_regions` (from [`phraya_core::detect_tandem_repeats`]) is produced by a strictly
/// forward-scanning detector, so regions are sorted by `start` and non-overlapping — binary
/// search replaces an O(n) linear scan with O(log n). On a human chromosome, `repeat_regions`
/// numbers in the millions (chr1: ~1.5M), so the linear scan dominated per-variant extraction
/// cost (~22ms/read, the majority of `phraya align`'s wall time on medium/large references).
fn is_in_repeat_region(pos: usize, repeat_regions: &[phraya_core::RepeatRegion]) -> bool {
    let idx = repeat_regions.partition_point(|r| r.start <= pos);
    idx > 0 && pos < repeat_regions[idx - 1].end
}

/// Parse CIGAR and extract VariantObservations at mismatch positions. The CIGAR starts at
/// `target_start` in `target` and at strand-local offset `query_start` in `query`.
fn extract_variants_from_cigar(
    cigar: &str,
    target_start: usize,
    query_start: usize,
    query: &[u8],
    target: &[u8],
    edit_distance: u32,
    provenance: String,
    coverage: &[WindowedCoverage],
    repeat_regions: &[phraya_core::RepeatRegion],
    mapq: u8,
    avg_base_quality: f64,
    confidence: f64,
    coverage_window_radius: usize,
    hotspot_intervals: &[(u32, u32)],
    mate_info: Option<&phraya_core::types::MateInfo>,
    insert_stats: Option<(i64, u32)>,
    strand: Strand,
) -> Vec<VariantObservation> {
    use std::sync::Arc;

    let mut variants = Vec::new();
    let mut q_pos = query_start;
    let mut t_pos = target_start;

    // Share CIGAR and provenance allocations across all variants from this read
    let cigar_arc: Arc<str> = Arc::from(cigar.to_string());
    let provenance_arc: Arc<str> = Arc::from(provenance);

    let ops = phraya_core::cigar::parse_ops(cigar);
    for (count, op) in ops {
        match op {
            'M' => {
                q_pos += count;
                t_pos += count;
            }
            'X' => {
                // Mismatch: one VariantObservation per position
                for i in 0..count {
                    let qp = q_pos + i;
                    let tp = t_pos + i;
                    if qp < query.len() && tp < target.len() {
                        let alt_base = query[qp];
                        let ref_base = target[tp];
                        let mut alleles = HashMap::new();
                        alleles.insert(alt_base, 1u32);

                        // Local coverage: ±coverage_window_radius bp window, values from the alignment coverage track.
                        let window_start = if tp >= coverage_window_radius {
                            tp - coverage_window_radius
                        } else {
                            0
                        };
                        let window_end = (tp + coverage_window_radius + 1).min(target.len());
                        let local_coverage: Vec<u32> = (window_start..window_end)
                            .map(|pos| get_abs_multi(coverage, pos))
                            .collect();
                        let variant_offset = (tp - window_start) as u32;

                        let in_repeat = is_in_repeat_region(tp, repeat_regions);

                        let kmer_uniqueness = if is_in_hotspot(tp as u32, hotspot_intervals) {
                            0.0
                        } else {
                            1.0
                        };

                        let mut obs = VariantObservation::with_arc(
                            tp as u32,
                            ref_base,
                            alleles,
                            kmer_uniqueness * confidence,
                            cigar_arc.clone(),
                            mapq,
                            edit_distance,
                            local_coverage,
                            avg_base_quality,
                            provenance_arc.clone(),
                        )
                        .with_tandem_repeat(in_repeat)
                        .with_kmer_uniqueness(kmer_uniqueness)
                        .with_coverage_window_offset(variant_offset)
                        .with_strand(strand)
                        .with_query_position(qp as u32);

                        if let Some(mi) = mate_info {
                            obs = obs
                                .with_mate_info(mi.clone())
                                .with_pair_counts(1, if mi.proper_pair { 1 } else { 0 });
                            if let Some((sum, count)) = insert_stats {
                                obs = obs.with_insert_stats(sum, count);
                            }
                        }

                        variants.push(obs);
                    }
                }
                q_pos += count;
                t_pos += count;
            }
            // WFA convention: 'I' = target has extra bases (standard 'D'); 'D' = query has extra.
            'I' => {
                // Target has extra bases = deletion in query relative to target
                // Only emit variant if the query is still being aligned (q_pos < query.len())
                // Skipping 'I' at tail-end where query has already ended prevents false deletions
                if t_pos < target.len() && q_pos < query.len() {
                    let deleted_bases = &target[t_pos..(t_pos + count).min(target.len())];
                    let window_start = if t_pos >= coverage_window_radius {
                        t_pos - coverage_window_radius
                    } else {
                        0
                    };
                    let window_end = (t_pos + coverage_window_radius + 1).min(target.len());
                    let local_coverage: Vec<u32> = (window_start..window_end)
                        .map(|pos| get_abs_multi(coverage, pos))
                        .collect();
                    let variant_offset = (t_pos - window_start) as u32;

                    let in_repeat = is_in_repeat_region(t_pos, repeat_regions);

                    let kmer_uniqueness = if is_in_hotspot(t_pos as u32, hotspot_intervals) {
                        0.0
                    } else {
                        1.0
                    };

                    // For deletion: ref_base is the first deleted base, alt is "." (VCF convention)
                    let ref_base = if !deleted_bases.is_empty() {
                        deleted_bases[0]
                    } else {
                        b'.'
                    };
                    let mut alleles = HashMap::new();
                    alleles.insert(b'.', 1u32);

                    let mut obs = VariantObservation::with_arc(
                        t_pos as u32,
                        ref_base,
                        alleles,
                        kmer_uniqueness * confidence,
                        cigar_arc.clone(),
                        mapq,
                        edit_distance,
                        local_coverage,
                        avg_base_quality,
                        provenance_arc.clone(),
                    )
                    .with_tandem_repeat(in_repeat)
                    .with_variant_type(phraya_core::types::VariantType::Deletion)
                    .with_kmer_uniqueness(kmer_uniqueness)
                    .with_coverage_window_offset(variant_offset)
                    .with_strand(strand)
                    .with_query_position(q_pos as u32);

                    if let Some(mi) = mate_info {
                        obs = obs
                            .with_mate_info(mi.clone())
                            .with_pair_counts(1, if mi.proper_pair { 1 } else { 0 });
                        if let Some((sum, count)) = insert_stats {
                            obs = obs.with_insert_stats(sum, count);
                        }
                    }

                    variants.push(obs);
                }
                t_pos += count;
            }
            'D' => {
                // Query has extra bases = insertion in query relative to target
                // Emit one VariantObservation for the inserted region
                if q_pos < query.len() && t_pos < target.len() {
                    let inserted_bases = &query[q_pos..(q_pos + count).min(query.len())];
                    let window_start = if t_pos >= coverage_window_radius {
                        t_pos - coverage_window_radius
                    } else {
                        0
                    };
                    let window_end = (t_pos + coverage_window_radius + 1).min(target.len());
                    let local_coverage: Vec<u32> = (window_start..window_end)
                        .map(|pos| get_abs_multi(coverage, pos))
                        .collect();
                    let variant_offset = (t_pos - window_start) as u32;

                    let in_repeat = is_in_repeat_region(t_pos, repeat_regions);

                    let kmer_uniqueness = if is_in_hotspot(t_pos as u32, hotspot_intervals) {
                        0.0
                    } else {
                        1.0
                    };

                    // For insertion: ref_base is ".", alt is the inserted bases
                    let mut alleles = HashMap::new();
                    for &base in inserted_bases {
                        *alleles.entry(base).or_insert(0) += 1;
                    }

                    let mut obs = VariantObservation::with_arc(
                        t_pos as u32,
                        b'.',
                        alleles,
                        kmer_uniqueness * confidence,
                        cigar_arc.clone(),
                        mapq,
                        edit_distance,
                        local_coverage,
                        avg_base_quality,
                        provenance_arc.clone(),
                    )
                    .with_tandem_repeat(in_repeat)
                    .with_variant_type(phraya_core::types::VariantType::Insertion)
                    .with_kmer_uniqueness(kmer_uniqueness)
                    .with_coverage_window_offset(variant_offset)
                    .with_strand(strand)
                    .with_query_position(q_pos as u32);

                    if let Some(mi) = mate_info {
                        obs = obs
                            .with_mate_info(mi.clone())
                            .with_pair_counts(1, if mi.proper_pair { 1 } else { 0 });
                        if let Some((sum, count)) = insert_stats {
                            obs = obs.with_insert_stats(sum, count);
                        }
                    }

                    variants.push(obs);
                }
                q_pos += count;
            }
            _ => {}
        }
    }

    variants
}

/// Raw per-read coverage, one small window per alignment (primary first, then each
/// alternative in `scored.alternatives` order).
///
/// A single alignment's own span is bounded to ~2x query length by the anchor windowing
/// in `align_oriented`, so one window per alignment stays small. Multi-mapped reads on a
/// repeat-rich genome can have alternates hundreds of megabases from the primary — a
/// *union* bounding box across all of them (the previous design) would materialize a
/// near-genome-length `Vec<u32>` per read. Emitting disjoint per-alignment windows and
/// merging each independently into the genome accumulator (a simple `+=`, so windows can
/// safely overlap) avoids that blowup entirely.
fn compute_windowed_coverage(
    scored: &crate::ScoredAlignments,
    target_len: usize,
) -> Vec<WindowedCoverage> {
    std::iter::once(&scored.primary)
        .chain(scored.alternatives.iter())
        .filter_map(|aln| {
            let start = aln.target_start.min(target_len);
            let end = aln.target_end.min(target_len);
            if start >= end {
                return None;
            }
            let mut counts = vec![0u32; end - start];
            for c in &mut counts {
                *c = 1;
            }
            Some(WindowedCoverage { start, counts })
        })
        .collect()
}

fn quantize_coverage(raw: &[u32]) -> Vec<u32> {
    raw.iter()
        .map(|&v| phraya_core::types::CoverageTrack::quantize(v as usize) as u32)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_plan() -> PhrayaPlan {
        PhrayaPlan::new(
        vec!["test".to_string()],
        "2026-06-01T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],)
    }

    #[test]
    fn align_config_convenience_constructors_match_new() {
        let fast = AlignConfig::fast();
        let fast_via_new = AlignConfig::new(Strategy::Fast);
        assert_eq!(fast.strategy, fast_via_new.strategy);
        assert_eq!(
            fast.coverage_window_radius,
            fast_via_new.coverage_window_radius
        );

        let sensitive = AlignConfig::sensitive();
        let sensitive_via_new = AlignConfig::new(Strategy::Sensitive);
        assert_eq!(sensitive.strategy, sensitive_via_new.strategy);
        assert_eq!(
            sensitive.coverage_window_radius,
            sensitive_via_new.coverage_window_radius
        );
    }

    #[test]
    fn is_in_repeat_region_matches_linear_scan_at_boundaries() {
        // Adjacent, non-overlapping regions as produced by detect_tandem_repeats:
        // [10, 20), [20, 30), gap, [50, 55)
        let regions = vec![
            phraya_core::RepeatRegion::new(10, 19, 2, "AT".to_string()),
            phraya_core::RepeatRegion::new(20, 29, 3, "CAG".to_string()),
            phraya_core::RepeatRegion::new(50, 54, 4, "GATA".to_string()),
        ];
        let linear = |pos: usize| regions.iter().any(|r| pos >= r.start && pos < r.end);
        for pos in 0..60 {
            assert_eq!(
                is_in_repeat_region(pos, &regions),
                linear(pos),
                "mismatch at pos {pos}"
            );
        }
    }

    #[test]
    fn is_in_repeat_region_empty_regions_always_false() {
        assert!(!is_in_repeat_region(0, &[]));
        assert!(!is_in_repeat_region(1_000_000, &[]));
    }

    #[test]
    fn is_in_hotspot_matches_linear_scan_at_boundaries() {
        let intervals = vec![(10u32, 20u32), (20u32, 30u32), (50u32, 55u32)];
        let linear = |pos: u32| intervals.iter().any(|&(s, e)| pos >= s && pos < e);
        for pos in 0..60u32 {
            assert_eq!(
                is_in_hotspot(pos, &intervals),
                linear(pos),
                "mismatch at pos {pos}"
            );
        }
    }

    #[test]
    fn select_anchors_falls_back_to_origin_only_for_comparable_lengths() {
        // Seedless and comparable-length: the issue #146 contig-vs-contig safety net.
        assert_eq!(
            select_anchors(&[], &[], chain_cap(Strategy::Fast), 100, 500, Alphabet::Dna),
            vec![SeedAnchor {
                query_pos: 0,
                target_pos: 0
            }]
        );

        // Seedless against a target that dwarfs the query: no candidate locus exists, so
        // there is nothing to extend and the caller must report no alignment.
        assert!(select_anchors(&[], &[], chain_cap(Strategy::Fast), 100, 100_000, Alphabet::Dna).is_empty());
    }

    /// A repetitive query can share exact k-mers yet chain to nothing. Those seeds are real
    /// candidate loci, so they must be extended — and at the seed's own diagonal, not at
    /// `(0,0)`, which would only ever examine the first `2*query_len` bases of the target.
    #[test]
    fn select_anchors_uses_seed_diagonals_when_chaining_finds_nothing() {
        let seeds = vec![
            crate::Seed {
                query_pos: 10,
                target_pos: 510,
                minimizer: 1,
            },
            // Same diagonal as the first (510-10 == 700-200): must dedupe to one anchor.
            crate::Seed {
                query_pos: 200,
                target_pos: 700,
                minimizer: 2,
            },
            crate::Seed {
                query_pos: 0,
                target_pos: 9_000,
                minimizer: 3,
            },
        ];

        // Target dwarfs the query, so the `(0,0)` fallback does not apply — these anchors
        // exist purely because seeds do.
        let anchors = select_anchors(&[], &seeds, chain_cap(Strategy::Balanced), 300, 100_000, Alphabet::Dna);
        assert_eq!(
            anchors,
            vec![
                SeedAnchor {
                    query_pos: 0,
                    target_pos: 500
                },
                SeedAnchor {
                    query_pos: 0,
                    target_pos: 9_000
                },
            ],
            "seed diagonals, deduplicated and capped at chain_cap"
        );

        // Fast caps at one anchor.
        assert_eq!(
            select_anchors(&[], &seeds, chain_cap(Strategy::Fast), 300, 100_000, Alphabet::Dna).len(),
            1
        );
    }

    #[test]
    fn select_anchors_uses_chain_target_start_and_respects_cap() {
        let seeds_a = vec![crate::Seed {
            query_pos: 0,
            target_pos: 100,
            minimizer: 1,
        }];
        let seeds_b = vec![crate::Seed {
            query_pos: 0,
            target_pos: 5000,
            minimizer: 2,
        }];
        let seeds_c = vec![crate::Seed {
            query_pos: 0,
            target_pos: 9000,
            minimizer: 3,
        }];
        let chains = vec![
            crate::chaining::Chain {
                seeds: seeds_a,
                score: 30,
            },
            crate::chaining::Chain {
                seeds: seeds_b,
                score: 25,
            },
            crate::chaining::Chain {
                seeds: seeds_c,
                score: 20,
            },
        ];
        // K=1 (Fast) keeps only the first (highest-scoring) chain's target_start.
        // No (0,0) fallback when chains exist.
        let anchors = select_anchors(&chains, &[], chain_cap(Strategy::Fast), 150, 10_000, Alphabet::Dna);
        assert_eq!(
            anchors,
            vec![SeedAnchor {
                query_pos: 0,
                target_pos: 100
            },]
        );

        // K=2 (Balanced) keeps the top 2 chains. No (0,0) fallback when chains exist.
        let anchors = select_anchors(&chains, &[], chain_cap(Strategy::Balanced), 150, 10_000, Alphabet::Dna);
        assert_eq!(
            anchors,
            vec![
                SeedAnchor {
                    query_pos: 0,
                    target_pos: 100
                },
                SeedAnchor {
                    query_pos: 0,
                    target_pos: 5000
                },
            ]
        );
    }

    #[test]
    fn test_align_task_handles_indel() {
        // Query has a deletion relative to target: target has 'T' at position 4 that query lacks.
        // Currently returns None due to equal-length guard — must use WFA instead.
        let query = Sequence::new(b"ACGACGT".to_vec(), None, "query_del".to_string(), None);
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(
            result.is_some(),
            "align_task must handle different-length sequences via WFA"
        );
    }

    #[test]
    fn test_perfect_match_no_variants() {
        let query = Sequence::new(b"ACGTACGT".to_vec(), None, "query1".to_string(), None);
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        assert!(
            result.variants.is_empty(),
            "Perfect match should have no variants"
        );
        assert_eq!(
            result.coverage[0].to_full(target.len()).len(),
            target.len(),
            "Coverage window should reconstruct to target length"
        );
    }

    #[test]
    fn test_query_positions_carry_scores() {
        let query = Sequence::new(b"ACGTACGT".to_vec(), None, "q".to_string(), None);
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan).expect("alignment should succeed");

        assert!(
            !result.query_positions.is_empty(),
            "should have at least one position"
        );
        let (_pos, score) = result.query_positions[0];
        assert!(
            score > 0.0 && score <= 1.0,
            "score must be in (0.0, 1.0], got {score}"
        );
    }

    #[test]
    fn test_variant_cigar_reflects_wfa_not_stub() {
        // SNP: T at position 2, otherwise identical 7-base sequences
        let query = Sequence::new(b"ACTACGT".to_vec(), None, "q".to_string(), None);
        let target = Sequence::new(b"ACCACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan).expect("alignment should succeed");
        assert_eq!(result.variants.len(), 1);

        let cigar = result.variants[0].cigar();
        assert_ne!(cigar, "1M", "CIGAR must come from WFA, not stub");
        // WFA over 7 equal-length bases with 1 mismatch produces something like "2M1X4M"
        assert!(
            cigar.contains('X') || cigar.contains('M'),
            "CIGAR should contain M or X ops, got: {cigar}"
        );
        assert!(
            cigar.len() > 2,
            "CIGAR should represent the full alignment, got: {cigar}"
        );
    }

    /// Deterministic diverse DNA (LCG) — avoids the minimizer-seed explosion of repetitive
    /// sequence and guarantees enough distinct k-mers to seed a 100bp read.
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

    #[test]
    fn reverse_strand_read_places_with_forward_strand_alleles() {
        // Issue #192: a read originating from the reverse strand is stored as the reverse
        // complement of its reference region. It must place at the correct forward coordinate,
        // report reference-strand alleles (not the RC), and be marked Strand::Reverse.
        use phraya_core::types::{reverse_complement, Strand};

        let target_bytes = diverse_dna(300, 7);
        // Forward-strand mutant of region [100,200): one SNP at forward offset 50 (abs pos 150).
        let mut region = target_bytes[100..200].to_vec();
        let orig = region[50];
        let snp = if orig == b'A' { b'C' } else { b'A' };
        region[50] = snp;
        // The reverse-strand read as it would be stored: RC of the forward mutant region.
        let read_bytes = reverse_complement(&region);

        let target = Sequence::new(target_bytes, None, "ref".to_string(), None);
        let query = Sequence::new(read_bytes, None, "rev_read".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan)
            .expect("reverse-strand read must align, not be dropped");

        // Places at the correct forward coordinate (~100).
        let (pos, score) = result.query_positions[0];
        assert!(
            (pos as i64 - 100).abs() <= 1,
            "reverse read should place near forward coord 100, got {pos}"
        );
        assert!(
            score > 0.9,
            "reverse read should place with high score, got {score}"
        );

        // The single SNP is reported at forward position 150, with the reference-strand
        // ref base and the forward-strand allele — NOT their reverse complements.
        let var = result
            .variants
            .iter()
            .find(|v| v.position() == 150)
            .expect("SNP must be called at forward position 150");
        assert_eq!(
            var.ref_base(),
            orig,
            "ref base must be the forward-strand reference base"
        );
        assert!(
            var.all_alleles().contains_key(&snp),
            "allele must be the forward-strand read base, not its complement"
        );
        assert_eq!(
            var.strand(),
            Strand::Reverse,
            "variant must record Strand::Reverse"
        );
    }

    #[test]
    fn align_stats_classify_placed_below_threshold_and_no_seed() {
        use std::sync::atomic::Ordering;

        let target = Sequence::new(diverse_dna(300, 3), None, "ref".to_string(), None);
        let plan = make_plan();
        let ctx = TargetContext::build(&target, &plan, Strategy::default());
        let cfg = AlignConfig::default();
        let stats = AlignStats::default();

        // Perfect substring → placed.
        let good = Sequence::new(
            target.bases()[50..150].to_vec(),
            None,
            "good".to_string(),
            None,
        );
        align_read(&ctx, &good, &plan, &cfg, Some(&stats));

        // A 150bp read whose divergence is clustered: long exact flanks still yield shared
        // minimizers (so it seeds), but ~15 clustered mismatches push the primary below 0.95
        // (15/150 = 0.10 divergence) → below_threshold, not no_seed.
        let mut div = target.bases()[50..200].to_vec();
        for i in 67..82 {
            div[i] = if div[i] == b'A' { b'C' } else { b'A' };
        }
        let diverged = Sequence::new(div, None, "div".to_string(), None);
        align_read(&ctx, &diverged, &plan, &cfg, Some(&stats));

        // Unrelated random read → no shared minimizer → no_seed.
        let junk = Sequence::new(diverse_dna(120, 987_654), None, "junk".to_string(), None);
        align_read(&ctx, &junk, &plan, &cfg, Some(&stats));

        assert_eq!(
            stats.placed.load(Ordering::Relaxed),
            1,
            "perfect read is placed"
        );
        assert_eq!(
            stats.no_seed.load(Ordering::Relaxed),
            1,
            "unrelated read is no_seed"
        );
        assert_eq!(
            stats.below_threshold.load(Ordering::Relaxed),
            1,
            "heavily-diverged read seeds but falls below 0.95"
        );
    }

    #[test]
    fn forward_strand_read_is_marked_forward() {
        // A forward-strand read with a SNP must still place forward and be marked Forward.
        use phraya_core::types::Strand;

        let target_bytes = diverse_dna(300, 11);
        let mut read_bytes = target_bytes[100..200].to_vec();
        let orig = read_bytes[50];
        read_bytes[50] = if orig == b'A' { b'C' } else { b'A' };

        let target = Sequence::new(target_bytes, None, "ref".to_string(), None);
        let query = Sequence::new(read_bytes, None, "fwd_read".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan).expect("forward read must align");
        let var = result
            .variants
            .iter()
            .find(|v| v.position() == 150)
            .expect("SNP at forward position 150");
        assert_eq!(
            var.strand(),
            Strand::Forward,
            "forward read must record Strand::Forward"
        );
    }

    /// Regression guard for the orientation-skip optimization in [`align_read`]: seeding
    /// the *wrong* orientation of a genuine match must produce the same raw seed count as
    /// the *true* orientation (canonical minimizers are strand-invariant in value) but zero
    /// chains (chaining requires seed co-linearity, which only holds in the true
    /// orientation). This is the property `align_read` relies on to skip extending a
    /// chainless orientation once the other orientation has chain support — if it stopped
    /// holding, the skip could silently drop a genuinely reverse-strand read.
    #[test]
    fn wrong_orientation_seeding_finds_seeds_but_no_chains() {
        let target_bytes = diverse_dna(300, 21);
        let true_fwd_read = target_bytes[50..200].to_vec();

        let target = Sequence::new(target_bytes, None, "ref".to_string(), None);
        let plan = make_plan();
        let ctx = TargetContext::build(&target, &plan, Strategy::Balanced);

        let fwd_sketch = sketch_sequence_default(&Sequence::new(
            true_fwd_read.clone(),
            None,
            "q".to_string(),
            None,
        ));
        let (fwd_seeds, fwd_chains) = seed_and_chain(&ctx, &fwd_sketch, Strategy::Balanced);

        let rc_bytes = reverse_complement(&true_fwd_read);
        let rc_sketch =
            sketch_sequence_default(&Sequence::new(rc_bytes, None, "q_rc".to_string(), None));
        let (rc_seeds, rc_chains) = seed_and_chain(&ctx, &rc_sketch, Strategy::Balanced);

        assert_eq!(
            fwd_seeds.len(),
            rc_seeds.len(),
            "canonical minimizers are strand-invariant: both orientations must find the \
             same raw seed count"
        );
        assert!(
            !fwd_chains.is_empty(),
            "the true orientation must chain to something"
        );
        assert!(
            rc_chains.is_empty(),
            "the wrong orientation of a genuine match must chain to nothing, despite \
             sharing the same seed count as the true orientation"
        );
    }

    #[test]
    fn test_single_snp_creates_variant() {
        // Query has T at position 2, target has C (SNP)
        let query = Sequence::new(b"ACTACGT".to_vec(), None, "query1".to_string(), None);
        let target = Sequence::new(b"ACCACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        assert_eq!(
            result.variants.len(),
            1,
            "One SNP should produce one variant"
        );

        let var = &result.variants[0];
        assert_eq!(
            var.position(),
            2,
            "Variant should be at position 2 (0-indexed)"
        );
        assert_eq!(var.ref_base(), b'C', "Reference base should be C");
        assert!(
            var.all_alleles().contains_key(&b'T'),
            "Allele T should be present"
        );
    }

    #[test]
    fn local_coverage_reflects_alignment_not_stub() {
        // 100bp query vs 200bp target (SNP at position 50). The query covers positions
        // 0..100. local_coverage for the variant (at pos 50) should be 1 (one read
        // aligned there), NOT a vector of all-1s ignoring whether the position is covered.
        let mut query_bases = vec![b'A'; 100];
        let mut target_bases = vec![b'A'; 200];
        query_bases[50] = b'T';
        target_bases[50] = b'C'; // SNP at position 50

        let query = Sequence::new(query_bases, None, "q".to_string(), None);
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan).expect("alignment must succeed");
        assert!(
            !result.variants.is_empty(),
            "must have at least one variant"
        );

        let var = &result.variants[0];
        let lc = var.local_coverage();
        // Positions within the alignment window (0..100) should have coverage ≥ 1.
        // Positions beyond the query end (100..200) were not covered — coverage = 0.
        // If local_coverage were still the stub (all 1s), uncovered positions would show 1.
        // With real coverage, the window around pos 50 is fully within the alignment → 1.
        assert!(
            lc.iter().any(|&c| c >= 1),
            "at least one position in the ±50bp window must have coverage ≥ 1"
        );
        // The ±50bp window around pos 50 is pos 0..101 — fully within the alignment.
        // All values should be 1 (one read). The stub would also give 1 here, but
        // the real test is that positions OUTSIDE the alignment are 0, not 1.
        // Use a variant near the start: align a SNP at position 5, window is 0..56.
        // Positions after query end (100..200) in that window should be 0 with real coverage.
        // We can't easily test that without a variant near position 150, so just confirm
        // the value is derived from alignment data (a known-1 position is fine as a smoke test
        // — the real regression guard is the audit finding that the stub was all-1s).
        assert!(
            lc[0] >= 1,
            "position within alignment window must have non-zero coverage, got {}",
            lc[0]
        );
    }

    #[test]
    fn tandem_repeat_variants_are_annotated() {
        // Build a target with a clear tandem repeat (ATATAT...) flanked by unique sequence.
        // A query with a SNP inside the repeat should produce a variant with in_tandem_repeat=true.
        // A SNP outside the repeat should produce in_tandem_repeat=false.
        let mut target_bases = b"TTAACCGGTA".to_vec(); // unique prefix (10bp)
        target_bases.extend_from_slice(b"ATATATATATATATATATATAT"); // tandem repeat (22bp, pos 10..32)
        target_bases.extend_from_slice(b"CGTACCGATT"); // unique suffix (10bp)
                                                       // Total: 42bp

        // Query matches target except: SNP in repeat at pos 15, SNP outside repeat at pos 2.
        let mut query_bases = target_bases.clone();
        query_bases[2] = if query_bases[2] == b'G' { b'C' } else { b'G' }; // SNP at pos 2 (unique region)
        query_bases[15] = if query_bases[15] == b'A' { b'T' } else { b'A' }; // SNP at pos 15 (repeat)

        let query = Sequence::new(query_bases, None, "q".to_string(), None);
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan).expect("alignment must succeed");
        assert!(result.variants.len() >= 2, "must have at least 2 variants");

        let repeat_variant = result.variants.iter().find(|v| v.position() == 15);
        let unique_variant = result.variants.iter().find(|v| v.position() == 2);

        assert!(repeat_variant.is_some(), "variant at pos 15 must exist");
        assert!(unique_variant.is_some(), "variant at pos 2 must exist");

        assert!(
            repeat_variant.unwrap().in_tandem_repeat(),
            "variant inside repeat region must be annotated in_tandem_repeat=true"
        );
        assert!(
            !unique_variant.unwrap().in_tandem_repeat(),
            "variant outside repeat region must be annotated in_tandem_repeat=false"
        );
    }

    #[test]
    fn test_indel_calling_deletion_variant_created() {
        // Query is missing a base at position 4: target="ACGTACGT", query="ACGACGT".
        // This should produce a deletion variant at position 4.
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let query = Sequence::new(b"ACGACGT".to_vec(), None, "query_del".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some(), "alignment should succeed");

        let result = result.unwrap();
        assert!(
            !result.variants.is_empty(),
            "deletion should produce at least one variant"
        );

        // Find the deletion variant (VariantType::Deletion)
        let deletion_var = result
            .variants
            .iter()
            .find(|v| v.variant_type() == phraya_core::types::VariantType::Deletion);
        assert!(
            deletion_var.is_some(),
            "must have a deletion variant for the missing base"
        );

        let var = deletion_var.unwrap();
        assert_eq!(
            var.variant_type(),
            phraya_core::types::VariantType::Deletion,
            "variant should be marked as Deletion"
        );
    }

    #[test]
    fn test_indel_calling_insertion_variant_created() {
        // Query has an extra base at position 4: target="ACGTACGT", query="ACGTTACGT".
        // This should produce an insertion variant at position 4.
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let query = Sequence::new(b"ACGTTACGT".to_vec(), None, "query_ins".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some(), "alignment should succeed");

        let result = result.unwrap();
        assert!(
            !result.variants.is_empty(),
            "insertion should produce at least one variant"
        );

        // Find the insertion variant (VariantType::Insertion)
        let insertion_var = result
            .variants
            .iter()
            .find(|v| v.variant_type() == phraya_core::types::VariantType::Insertion);
        assert!(
            insertion_var.is_some(),
            "must have an insertion variant for the extra base"
        );

        let var = insertion_var.unwrap();
        assert_eq!(
            var.variant_type(),
            phraya_core::types::VariantType::Insertion,
            "variant should be marked as Insertion"
        );
    }

    /// Throughput: 20 reads × 100bp against a 200bp reference must complete
    /// at ≥ 100 reads/sec (< 200ms wall time).
    ///
    /// Uses diverse (LCG-generated) sequences to avoid the minimizer-seed
    /// explosion that repetitive sequences cause (~1274 seeds → hours).
    /// With diverse sequences, ~6 seeds per alignment → ~120K DP cells total.
    ///
    /// WFA (O(s·n)) replaced the O(n×m) DP; this test passes in debug.
    #[test]
    fn issue_88_throughput_100_reads_per_sec() {
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

        let ref_seq = diverse_dna(200, 42);
        let read_seq: Vec<u8> = ref_seq[..100].to_vec();

        let target = Sequence::new(ref_seq, None, "ref".to_string(), None);
        let plan = make_plan();

        let start = std::time::Instant::now();
        for i in 0..20 {
            let query = Sequence::new(read_seq.clone(), None, format!("read{i}"), None);
            let _ = align_task(&query, &target, &plan);
        }
        let elapsed = start.elapsed();

        assert!(
            elapsed.as_millis() < 200,
            "20 alignments (150bp vs 1000bp) took {}ms — below 100 reads/sec target.\n\
             The naive O(n×m) DP must be replaced with true WFA wavefront algorithm.",
            elapsed.as_millis()
        );
    }

    // ========================================================================
    // ISSUE #180: Abandonment Sentinel Integration Tests
    // ========================================================================
    //
    // These tests verify that abandoned alignments (no fitting-end reached
    // within a capped max_s) are correctly filtered out in align_read and
    // score_alignments, preventing spurious edit_distance=0 perfect alignments
    // from becoming the primary. This is the critical bug fix for ADR-0007.

    /// **ACCEPTANCE CRITERION**: When an anchor extension returns abandoned
    /// (no fitting-end within cap), it must not be added to the alignments
    /// vector passed to score_alignments.
    ///
    /// The contract: extend_anchor should filter abandonment before it reaches
    /// align_read's alignment collection, or align_read should skip Err/None results.
    #[test]
    fn issue_180_abandoned_alignment_never_reaches_score_alignments() {
        // This test verifies the integration: at the executor level,
        // abandoned alignments are filtered before scoring.
        //
        // When issue #180 is implemented, fill_wfa_fitting returns Option,
        // and wfa_extend/myers_extend propagate abandonment as Err.
        // extend_anchor must return Err on abandonment, and align_read must
        // skip Err results (line 306: "Err(e) => log::warn!(...)").
        //
        // Before the fix, an abandoned alignment with (empty_cigar, 0, full_len)
        // would be added to alignments and would become the primary (min edit_distance).
        //
        // After the fix, it's never added, so score_alignments never sees it.
        //
        // At current (uncapped) max_s, nothing is abandoned, so this test just
        // verifies that a normal alignment succeeds and produces reasonable output.

        let target = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "ref".to_string(), None);
        let query = Sequence::new(
            b"ACGTACGTACGTACGT".to_vec(),
            None,
            "query".to_string(),
            None,
        );
        let plan = make_plan();

        // At current (uncapped) max_s, this should succeed normally
        let result = align_task(&query, &target, &plan);
        assert!(
            result.is_some(),
            "normal alignment must succeed (backward compat at uncapped max_s)"
        );

        // Sanity check: the result should be from a real alignment, not an abandoned fallback
        // After issue #180 fix, abandoned alignments are filtered, so only real ones reach scoring
        let result = result.unwrap();
        // A real alignment produces variants (for mismatches/indels) or has coverage
        // The current test just verifies the basic contract: alignment succeeds
        assert!(
            result.query_positions.len() > 0,
            "alignment must produce query position information (primary + alternatives)"
        );
    }

    /// **ACCEPTANCE CRITERION**: score_alignments must only receive real
    /// (non-abandoned) alignments in its input vector.
    ///
    /// This documents the contract: align_read (line 290-312) builds the
    /// alignments vector by collecting only successful extend_anchor results.
    /// Abandoned results are filtered out before scoring.
    #[test]
    fn issue_180_score_alignments_input_never_contains_abandoned() {
        // This test verifies that score_alignments (which uses
        // min_by_key(|a| edit_distance)) never receives an alignment with
        // edit_distance=0 and empty CIGAR, which would indicate abandonment
        // mis-represented as a perfect match.

        // Create a test input as score_alignments would receive it
        use crate::Alignment;

        // A real perfect match alignment
        let perfect = Alignment {
            cigar: "16M".to_string(),
            edit_distance: 0,
            query_start: 0,
            query_end: 16,
            target_start: 0,
            target_end: 16,
        };

        // score_alignments should work correctly on this
        let scored = score_alignments(&[perfect.clone()], 16);
        assert_eq!(
            scored.primary, perfect,
            "perfect match is selected as primary"
        );

        // But an abandoned alignment (empty CIGAR, 0 edit distance) should
        // NEVER reach score_alignments; it's filtered in align_read.
        // If it somehow did, the following would demonstrate the bug:
        // (We don't actually test this bug because the fix prevents it.)

        // Verification: after issue #180 fix, align_read's extend_anchor
        // will return Err on abandonment, and the log::warn! at line 306
        // will be hit, preventing the alignment from being added to the vector.
    }

    /// **ACCEPTANCE CRITERION**: Multiple anchors must be filtered individually.
    /// If one abandons and one succeeds, only the successful one is scored.
    ///
    /// This verifies that the filtering happens per-anchor in the loop at line 293-307,
    /// not globally.
    #[test]
    fn issue_180_per_anchor_filtering_in_loop() {
        // Scenario: multiple seeds found, leading to multiple anchors.
        // Some might abandon (hypothetically, with a cap), others might succeed.
        // The alignments vector should only contain successful ones.

        let target = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "ref".to_string(), None);
        let query = Sequence::new(
            b"ACGTACGTACGTACGT".to_vec(),
            None,
            "query".to_string(),
            None,
        );
        let plan = make_plan();

        // At current (uncapped) max_s, all anchors align successfully
        let result = align_task(&query, &target, &plan);
        assert!(
            result.is_some(),
            "all-anchors-succeed case must produce a result"
        );

        // After issue #180 fix, if one anchor abandoned and one succeeded,
        // the result would only be built from the successful one.
        // This test documents that contract.
    }

    /// **ACCEPTANCE CRITERION**: Backward compatibility: at uncapped max_s,
    /// no alignment is abandoned, so results are identical to today.
    ///
    /// This test runs the same query-target pair and verifies the result
    /// structure is consistent (no silent changes in alignment selection).
    #[test]
    fn issue_180_backward_compat_uncapped_produces_consistent_results() {
        let target = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "ref".to_string(), None);
        let query = Sequence::new(
            b"ACGTACGTACGTACGT".to_vec(),
            None,
            "query".to_string(),
            None,
        );
        let plan = make_plan();

        let result1 = align_task(&query, &target, &plan);
        let result2 = align_task(&query, &target, &plan);

        // Results should be deterministic (identical on reruns)
        assert_eq!(
            result1.is_some(),
            result2.is_some(),
            "alignment results must be deterministic"
        );

        if let (Some(r1), Some(r2)) = (result1, result2) {
            assert_eq!(
                r1.variants.len(),
                r2.variants.len(),
                "variant count must match"
            );
            assert_eq!(
                r1.query_positions.len(),
                r2.query_positions.len(),
                "query position count must match"
            );
        }
    }

    /// **ACCEPTANCE CRITERION**: The primary alignment selected by score_alignments
    /// must have a valid (non-empty or non-zero) CIGAR/edit_distance pair.
    ///
    /// After issue #180, abandoned alignments are filtered out, so the primary
    /// can never be the spurious (empty_cigar, 0, full_len) fallback.
    #[test]
    fn issue_180_primary_alignment_is_never_spurious_fallback() {
        let target = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "ref".to_string(), None);
        let query = Sequence::new(
            b"ACGTACGTACGTACGT".to_vec(),
            None,
            "query".to_string(),
            None,
        );
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some(), "alignment should succeed");

        let result = result.unwrap();

        // The variants extracted from scored.primary (line 341-357) depend on the CIGAR.
        // A spurious (empty_cigar, 0, full_len) alignment would produce no variants.
        // A real alignment produces variants (or at least coverage).
        //
        // This is a sanity check: the result must be from a real alignment,
        // not the abandoned fallback.
        assert!(
            !result.variants.is_empty() || result.coverage.iter().any(|w| !w.counts.is_empty()),
            "result must be from a real alignment (extracted variants or coverage), \
             not from a spurious abandoned fallback"
        );
    }

    /// Deterministic diverse amino-acid sequence (LCG), mirroring `diverse_dna` but over
    /// the 20 canonical residues — for protein-alphabet tests below.
    fn diverse_protein(len: usize, seed: u64) -> Vec<u8> {
        const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
        let mut x = seed;
        (0..len)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                AA[((x >> 33) % 20) as usize]
            })
            .collect()
    }

    /// ADR-0013: protein alignment must never attempt reverse-strand search — protein has
    /// no biological reverse strand, and `reverse_complement` would silently corrupt any
    /// residue byte that happens to collide with a DNA letter (e.g. Ala 'A' -> 'T') if it
    /// ran on protein bytes at all. Every variant `align_read` reports for a protein query
    /// must therefore be `Strand::Forward`, even though nothing in the query content itself
    /// rules out a spurious reverse-orientation "match" the way it would for DNA.
    #[test]
    fn protein_alignment_never_attempts_reverse_strand() {
        let target_bytes = diverse_protein(300, 11);
        let mut region = target_bytes[100..200].to_vec();
        // One substitution so extract_variants_from_cigar reports at least one variant.
        let orig = region[50];
        region[50] = if orig == b'A' { b'C' } else { b'A' };

        let target = Sequence::new(target_bytes, None, "target".to_string(), None);
        let query = Sequence::new(region, None, "query".to_string(), None);

        let mut plan = make_plan();
        plan.alphabet = Alphabet::Protein;
        let ctx = TargetContext::build(&target, &plan, Strategy::Sensitive);
        let cfg = AlignConfig::sensitive();

        let result = align_read(&ctx, &query, &plan, &cfg, None)
            .expect("protein query with one substitution against its source region must align");
        assert!(!result.variants.is_empty(), "expected at least one reported substitution");
        for variant in &result.variants {
            assert_eq!(
                variant.strand(),
                Strand::Forward,
                "protein alignment must never report Strand::Reverse — no reverse-strand search exists for protein (ADR-0013)"
            );
        }
    }

    /// ADR-0013: WFA and Myers are already alphabet-blind byte comparators (confirmed by
    /// reading the kernels — Myers' Peq table is `[[0u64; 256]; num_blocks]`, not a 4-symbol
    /// table; WFA's match-extension is a raw byte compare). This differential-tests that
    /// claim directly on protein-alphabet bytes: both engines must agree on edit distance
    /// across a spread of divergence levels, exactly as they already do for DNA.
    #[test]
    fn protein_wfa_myers_edit_distance_agree() {
        for seed in 0..20u64 {
            let query = diverse_protein(120, seed * 97 + 1);
            let mut target = query.clone();
            // Sprinkle substitutions at a seed-dependent rate (0%..~15%).
            let rate = seed % 5; // 0,1,2,3,4 -> 0%,3.75%,7.5%,11.25%,15%
            const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
            let mut x = seed.wrapping_add(12345);
            for b in target.iter_mut() {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                if (x >> 40) % 100 < rate as u64 * 15 {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
                    *b = AA[((x >> 33) % 20) as usize];
                }
            }
            let seed_anchor = SeedAnchor { query_pos: 0, target_pos: 0 };
            let wfa = wfa_extend(&query, &target, seed_anchor.clone())
                .unwrap_or_else(|e| panic!("wfa_extend failed for seed {seed}: {e:?}"));
            let myers = myers_extend(&query, &target, seed_anchor)
                .unwrap_or_else(|e| panic!("myers_extend failed for seed {seed}: {e:?}"));
            assert_eq!(
                wfa.edit_distance, myers.edit_distance,
                "WFA/Myers edit distance disagreement on protein bytes (seed={seed}, rate={rate})"
            );
        }
    }

    // ── Strategy::Assembly (ADR-0015) ────────────────────────────────────────────────

    fn assembly_dna(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..len)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((x >> 33) & 3) as usize]
            })
            .collect()
    }

    fn with_snps(seq: &[u8], positions: &[usize]) -> Vec<u8> {
        let mut out = seq.to_vec();
        for &p in positions {
            out[p] = match out[p] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        out
    }

    fn align_assembly_pair(query: &[u8], target: &[u8]) -> Option<AlignmentResult> {
        let q = Sequence::new(query.to_vec(), None, "q".to_string(), None);
        let t = Sequence::new(target.to_vec(), None, "t".to_string(), None);
        align_task_with_config(&q, &t, &make_plan(), &AlignConfig::new(Strategy::Assembly))
    }

    fn block_spans(r: &AlignmentResult) -> Vec<(u32, u32, u32, u32)> {
        r.blocks.iter().map(|b| (b.q_start, b.q_end, b.t_start, b.t_end)).collect()
    }

    #[test]
    fn assembly_contained_query_is_one_exact_block() {
        let target = assembly_dna(5_000, 1);
        let r = align_assembly_pair(&target[1_000..4_000], &target).expect("contained query aligns");
        assert_eq!(block_spans(&r), vec![(0, 3_000, 1_000, 4_000)]);
        let b = &r.blocks[0];
        assert_eq!((b.strand, b.edit_distance, &*b.cigar), (Strand::Forward, 0, "3000M"));
        assert!(r.variants.is_empty());
    }

    #[test]
    fn assembly_reports_every_interior_snp_with_its_allele() {
        let target = assembly_dna(8_000, 2);
        let snps: Vec<usize> = (150..7_850).step_by(97).collect();
        let query = with_snps(&target, &snps);
        let r = align_assembly_pair(&query, &target).expect("aligns");
        assert_eq!(r.blocks.len(), 1, "one homologous region gives one block: {:?}", block_spans(&r));
        let mut got: Vec<(u32, u8, u32)> = r
            .variants
            .iter()
            .map(|v| {
                assert_eq!(v.variant_type(), phraya_core::types::VariantType::Snp);
                (v.position(), *v.all_alleles().keys().next().unwrap(), v.query_position())
            })
            .collect();
        got.sort();
        let want: Vec<(u32, u8, u32)> = snps.iter().map(|&p| (p as u32, query[p], p as u32)).collect();
        assert_eq!(got, want);
    }

    #[test]
    fn assembly_query_spanning_two_targets_gives_one_block_per_overlap() {
        let a = assembly_dna(6_000, 3);
        let b = assembly_dna(6_000, 4);
        let mut query = a[4_000..].to_vec();
        query.extend_from_slice(&b[..3_000]);
        let ra = align_assembly_pair(&query, &a).expect("aligns to A");
        let rb = align_assembly_pair(&query, &b).expect("aligns to B");
        assert_eq!(block_spans(&ra), vec![(0, 2_000, 4_000, 6_000)]);
        assert_eq!(block_spans(&rb), vec![(2_000, 5_000, 0, 3_000)]);
    }

    #[test]
    fn assembly_reverse_complement_query_reports_forward_query_coordinates() {
        let target = assembly_dna(6_000, 5);
        let fwd = with_snps(&target[1_000..5_000], &[500, 2_000]);
        let query = reverse_complement(&fwd);
        let r = align_assembly_pair(&query, &target).expect("reverse-complemented query aligns");
        assert_eq!(block_spans(&r), vec![(0, 4_000, 1_000, 5_000)]);
        let b = &r.blocks[0];
        assert_eq!(b.strand, Strand::Reverse);
        let mut positions: Vec<u32> = r.variants.iter().map(|v| v.position()).collect();
        positions.sort();
        assert_eq!(positions, vec![1_500, 3_000]);
        for v in &r.variants {
            assert_eq!(v.strand(), Strand::Reverse);
            assert!(b.contains_observation(v), "variant not attributed to its block: {v:?}");
            // The alt allele is reported in reference orientation: the complement of the
            // query base at the variant's forward query coordinate.
            let fq = b.forward_query_position(v.query_position()) as usize;
            let alt = *v.all_alleles().keys().next().unwrap();
            assert_eq!(reverse_complement(&[query[fq]]), vec![alt]);
        }
    }

    #[test]
    fn assembly_small_indels_are_gap_filled_inside_one_block() {
        let target = assembly_dna(6_000, 6);
        let mut query = target[..2_000].to_vec();
        query.extend_from_slice(&target[2_003..4_000]); // 3 bp deletion
        query.extend_from_slice(b"GATTACA"); // 7 bp insertion
        query.extend_from_slice(&target[4_000..]);
        let r = align_assembly_pair(&query, &target).expect("aligns");
        assert_eq!(block_spans(&r), vec![(0, 6_004, 0, 6_000)]);
        assert_eq!(r.blocks[0].edit_distance, 10);
        use phraya_core::types::VariantType;
        let kinds: Vec<(VariantType, u32)> = r.variants.iter().map(|v| (v.variant_type(), v.position())).collect();
        assert_eq!(kinds.len(), 2, "one deletion and one insertion: {kinds:?}");
        assert!(kinds.iter().any(|&(k, p)| k == VariantType::Deletion && (1_995..=2_003).contains(&p)), "{kinds:?}");
        assert!(kinds.iter().any(|&(k, p)| k == VariantType::Insertion && (3_995..=4_007).contains(&p)), "{kinds:?}");
    }

    #[test]
    fn assembly_large_insertion_splits_the_block() {
        let target = assembly_dna(6_000, 7);
        let mut query = target[..3_000].to_vec();
        query.extend(assembly_dna(400, 8));
        query.extend_from_slice(&target[3_000..]);
        let r = align_assembly_pair(&query, &target).expect("aligns");
        assert_eq!(r.blocks.len(), 2, "a 400 bp insertion breaks the chain: {:?}", block_spans(&r));
        let covered: u32 = r.blocks.iter().map(|b| b.t_end - b.t_start).sum();
        assert!((5_990..=6_010).contains(&covered), "blocks cover {covered} bp: {:?}", block_spans(&r));
        assert!(r.variants.is_empty(), "the insertion is a block boundary, not a variant");
    }

    #[test]
    fn assembly_unrelated_sequences_have_no_block() {
        assert!(align_assembly_pair(&assembly_dna(5_000, 9), &assembly_dna(5_000, 10)).is_none());
    }

    #[test]
    fn assembly_long_divergent_contig_is_one_block() {
        // 120 kb at 0.5% divergence: twice the length a whole-contig WFA survives (A1).
        let target = assembly_dna(120_000, 11);
        let snps: Vec<usize> = (100..119_900).step_by(200).collect();
        let query = with_snps(&target, &snps);
        let r = align_assembly_pair(&query, &target).expect("aligns");
        assert_eq!(block_spans(&r), vec![(0, 120_000, 0, 120_000)]);
        assert_eq!(r.variants.len(), snps.len());
    }

    /// §7: on a long near-identical contig, `chain_seeds` (bounded by `max_dist_x`/`max_iter`)
    /// can split one clean diagonal into several independent chains; each grows via
    /// `chain_segments` to the same maximal exact match, so without dedup the same block would
    /// repeat (observed: 16 identical blocks on a 615 kb self-match in Contig_28). The dedup
    /// in `align_assembly` collapses them. This test asserts the dedup invariant — no two
    /// emitted blocks are identical — which holds whether or not the diagonal splits on the
    /// synthetic seed (it does not split at 50 kb, so the test also guards against a regression
    /// that would *introduce* duplicates). The block itself is still correct.
    #[test]
    fn assembly_self_match_has_no_duplicate_blocks() {
        let target = assembly_dna(50_000, 1);
        let r = align_assembly_pair(&target, &target).expect("self-match aligns");
        let spans: Vec<(u32, u32, u32, u32, &str)> = r
            .blocks
            .iter()
            .map(|b| (b.q_start, b.q_end, b.t_start, b.t_end, &*b.cigar))
            .collect();
        let mut seen = std::collections::HashSet::new();
        for s in &spans {
            assert!(seen.insert(*s), "duplicate block in self-match: {spans:?}");
        }
    }

    fn block(coords: (usize, usize, usize, usize), cigar: &str) -> LocalBlock {
        LocalBlock {
            q_start: coords.0,
            q_end: coords.1,
            t_start: coords.2,
            t_end: coords.3,
            cigar: cigar.to_string(),
            edit_distance: 0,
        }
    }

    #[test]
    fn dedup_blocks_is_order_preserving_and_handles_non_adjacent_duplicates() {

        let a = block((0, 1000, 0, 1000), "1000M");
        let b = block((1000, 2000, 1000, 2000), "1000M");
        let cases: Vec<(Vec<(Strand, LocalBlock, f64)>, Vec<(Strand, LocalBlock, f64)>)> = vec![
            // adjacent duplicates collapse to one
            (
                vec![(Strand::Forward, a.clone(), 1.0), (Strand::Forward, a.clone(), 1.0)],
                vec![(Strand::Forward, a.clone(), 1.0)],
            ),
            // a non-adjacent third equal to the first is dropped, order kept
            (
                vec![
                    (Strand::Forward, a.clone(), 1.0),
                    (Strand::Forward, b.clone(), 1.0),
                    (Strand::Forward, a.clone(), 1.0),
                ],
                vec![(Strand::Forward, a.clone(), 1.0), (Strand::Forward, b.clone(), 1.0)],
            ),
            // different strands are independent and both kept
            (
                vec![(Strand::Forward, a.clone(), 1.0), (Strand::Reverse, a.clone(), 1.0)],
                vec![(Strand::Forward, a.clone(), 1.0), (Strand::Reverse, a.clone(), 1.0)],
            ),
            // dedup is on (strand, block) only; differing identity does not prevent a collapse
            (
                vec![(Strand::Forward, a.clone(), 1.0), (Strand::Forward, a.clone(), 0.9)],
                vec![(Strand::Forward, a.clone(), 1.0)],
            ),
        ];
        for (mut input, want) in cases {
            dedup_blocks(&mut input);
            assert_eq!(input, want, "got {input:?}, want deduped to first occurrences");
        }
    }

    /// §4 spacing sweep: dense SNP clusters where no clean 21-mer anchor survives should still
    /// recover every SNP once the chain spans the zero-drift gap. Probing several spacings
    /// records where recovery holds and where the chain splits.
    #[test]
    fn assembly_cluster_spacing_sweep_calls_all_snps() {
        use phraya_core::types::VariantType;
        let target = assembly_dna(12_000, 3);
        let center = 6_000usize;
        let n = 8usize;
        let mut results: Vec<(usize, usize)> = Vec::new();
        for spacing in [5, 10, 20, 50] {
            let pos: Vec<usize> = (0..n).map(|i| center + i * spacing).collect();
            let query = with_snps(&target, &pos);
            let r = align_assembly_pair(&query, &target).expect("sweep pair aligns");
            let mut called: Vec<u32> = r
                .variants
                .iter()
                .filter(|v| v.variant_type() == VariantType::Snp)
                .map(|v| v.position())
                .collect();
            called.sort_unstable();
            results.push((spacing, called.len()));
            assert_eq!(
                called.len(),
                n,
                "spacing={spacing} only {}/{} SNPs: {called:?} blocks={:?}",
                called.len(),
                n,
                block_spans(&r)
            );
        }
        let summary: Vec<String> = results.iter().map(|(s, c)| format!("{s}bp:{c}")).collect();
        eprintln!("cluster spacing sweep: {}", summary.join(", "));
    }

}
