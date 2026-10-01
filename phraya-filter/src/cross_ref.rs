//! Cross-reference filter operations (ADR-0011, issue #202).
//!
//! These filters operate on the **superposition** stored in the cross-space
//! queries sidecar (`.phraya.queries`), enacting decisions that `align` deliberately
//! deferred: read depletion against a designated space, and best-space typing.
//!
//! Per the design philosophy in AGENTS.md, align surfaces possibilities and filter
//! decides. These operations read the absolute normalized identity per (read, space)
//! straight off the sidecar — exactly as the cross-reference superposition was designed
//! for.

use phraya_io::queries::{CrossSpacePlacement, CrossSpaceQueryIndex};

/// Decision outcome for cross-reference read depletion.
///
/// Depletion conservatively preserves ambiguous reads (near-ties across spaces)
/// and retains reads that don't confidently hit the target space.
#[derive(Debug, Clone, PartialEq)]
pub enum DepletionDecision {
    /// Read should be kept (not in the depletion target space, or ambiguous).
    Keep {
        /// Why the read was kept.
        reason: KeepReason,
    },
    /// Read should be depleted (dropped).
    Drop {
        /// Why the read was dropped.
        reason: DropReason,
    },
}

/// Why a read was kept during depletion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KeepReason {
    /// No placement in the depletion target space — can't deplete what isn't there.
    NotInTargetSpace,
    /// Read's best hit is in another space, so it's kept by default.
    BetterElsewhere,
    /// Margin between top two placements is below the ambiguity threshold —
    /// ambiguous, so we conservatively preserve it.
    Ambiguous { margin: f64 },
}

/// Why a read was dropped during depletion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    /// Read's best placement is confidently in the target space (identity high,
    /// margin to other spaces exceeds the ambiguity threshold).
    ConfidentHit,
    /// Read only has placements in the target space (no competing placement).
    SoleHit,
}

/// Decision outcome for best-space typing.
///
/// Typing assigns each read to its best-identity reference space and reports
/// the cross-space margin. Ambiguous and unassigned reads are surfaced, not
/// silently binned.
#[derive(Debug, Clone, PartialEq)]
pub enum TypingDecision {
    /// Read assigned to its best-identity space.
    ///
    /// `margin` is the identity gap between the best and second-best placement.
    /// A larger margin indicates higher confidence in the assignment.
    Assign {
        best_space: String,
        best_identity: f64,
        margin: f64,
    },
    /// Read has near-equally-good placements across spaces — cannot confidently assign.
    ///
    /// Both the best and second-best placements are reported so the caller can decide.
    Ambiguous {
        best_space: String,
        best_identity: f64,
        second_space: String,
        second_identity: f64,
        margin: f64,
    },
    /// Read has no placements passing the minimum identity threshold in any space.
    Unassigned,
}

/// Cross-reference filter for depletion and typing decisions.
///
/// Operates on the cross-space queries sidecar, using absolute normalized identity
/// per (read, space) to make decisions that `align` deliberately deferred.
///
/// # Design
///
/// The `ambiguity_margin` (default 0.02) is the identity gap below which a
/// near-tie between spaces is treated as ambiguous. This mirrors the existing
/// 0.95 score-ratio threshold's philosophy — hard-coded opinion with a knob
/// for override, matching the existing pattern used by `ThresholdFilter`.
///
/// For depletion: reads are only dropped when they confidently hit the target
/// space (identity ≥ minimum AND margin to other placements ≥ ambiguity_margin).
/// Ambiguous (near-tie) reads are always preserved — depletion must not discard
/// reads that could belong elsewhere.
#[derive(Debug, Clone)]
pub struct CrossRefFilter {
    /// Identity gap below which two placements are considered ambiguous (near-tie).
    /// Reads within this margin of each other are treated as ambiguous.
    pub ambiguity_margin: f64,
    /// Minimum identity (0–1) for a placement to be considered "real."
    /// Placements below this threshold are ignored in all decisions.
    pub min_identity: f64,
    /// The reference space to deplete against, if depletion mode is active.
    pub deplete_space: Option<String>,
}

impl Default for CrossRefFilter {
    fn default() -> Self {
        CrossRefFilter {
            ambiguity_margin: 0.02,
            min_identity: 0.95,
            deplete_space: None,
        }
    }
}

impl CrossRefFilter {
    /// Create a new cross-reference filter with default thresholds.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the ambiguity margin (identity gap below which placements are near-ties).
    pub fn with_ambiguity_margin(mut self, margin: f64) -> Self {
        self.ambiguity_margin = margin;
        self
    }

    /// Set the minimum identity threshold for considering a placement real.
    pub fn with_min_identity(mut self, identity: f64) -> Self {
        self.min_identity = identity;
        self
    }

    /// Set the reference space to deplete against.
    pub fn with_deplete_space(mut self, space: &str) -> Self {
        self.deplete_space = Some(space.to_string());
        self
    }

    /// Get placements passing the minimum identity threshold for a read.
    ///
    /// Returns placements sorted by descending identity (best first).
    fn filtered_placements<'a>(
        &self,
        placements: &'a [CrossSpacePlacement],
    ) -> Vec<&'a CrossSpacePlacement> {
        let mut kept: Vec<&CrossSpacePlacement> = placements
            .iter()
            .filter(|p| p.identity >= self.min_identity)
            .collect();
        kept.sort_by(|a, b| {
            b.identity
                .partial_cmp(&a.identity)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        kept
    }

    /// Classify a read for depletion: keep or drop based on confidence it belongs
    /// to the designated target space.
    ///
    /// Conservative: ambiguous (near-tie) reads are always preserved.
    ///
    /// # Arguments
    /// * `placements` - all placements for this read across all spaces
    ///
    /// # Returns
    /// `DepletionDecision::Drop` only when the read's best placement is in the
    /// target space with high confidence (margin to second-best exceeds ambiguity).
    /// Otherwise returns `DepletionDecision::Keep` with the reason.
    pub fn classify_depletion(
        &self,
        placements: &[CrossSpacePlacement],
    ) -> DepletionDecision {
        let target = match &self.deplete_space {
            Some(s) => s,
            None => {
                // No deplete target configured — keep everything.
                return DepletionDecision::Keep {
                    reason: KeepReason::NotInTargetSpace,
                };
            }
        };

        let real = self.filtered_placements(placements);

        if real.is_empty() {
            // No placement passes the identity threshold — there's nothing to deplete.
            return DepletionDecision::Keep {
                reason: KeepReason::NotInTargetSpace,
            };
        }

        let best = &real[0];

        // Check if the best placement is in the target space.
        if best.space != *target {
            // Best hit is elsewhere — keep the read.
            return DepletionDecision::Keep {
                reason: KeepReason::BetterElsewhere,
            };
        }

        // Best placement IS in the target space — check if it's unambiguous.
        match real.len() {
            1 => {
                // Only one placement, and it's in the target space — sole hit.
                // No competing placement to cause ambiguity, so deplete.
                DepletionDecision::Drop {
                    reason: DropReason::SoleHit,
                }
            }
            _ => {
                let second = &real[1];
                let margin = best.identity - second.identity;

                if margin >= self.ambiguity_margin {
                    // Confident that this read belongs to the target space.
                    DepletionDecision::Drop {
                        reason: DropReason::ConfidentHit,
                    }
                } else {
                    // Near-tie between target and another space — ambiguous, keep.
                    DepletionDecision::Keep {
                        reason: KeepReason::Ambiguous { margin },
                    }
                }
            }
        }
    }

    /// Classify a read for typing: assign to best-identity space, or flag as ambiguous/unassigned.
    ///
    /// Unlike depletion, typing always produces an assignment (or ambiguity surface)
    /// using the best hit, reporting the margin so the caller can apply their own
    /// cutoff if desired.
    ///
    /// # Arguments
    /// * `placements` - all placements for this read across all spaces
    ///
    /// # Returns
    /// - `Assign` if a single space has a confident best hit (margin ≥ ambiguity)
    /// - `Ambiguous` if the top two placements are within ambiguity_margin
    /// - `Unassigned` if no placement passes min_identity
    pub fn classify_typing(
        &self,
        placements: &[CrossSpacePlacement],
    ) -> TypingDecision {
        let real = self.filtered_placements(placements);

        if real.is_empty() {
            return TypingDecision::Unassigned;
        }

        let best = &real[0];

        match real.len() {
            1 => {
                // Only one placement passes the threshold — assign with infinite margin.
                TypingDecision::Assign {
                    best_space: best.space.clone(),
                    best_identity: best.identity,
                    margin: f64::INFINITY,
                }
            }
            _ => {
                let second = &real[1];
                let margin = best.identity - second.identity;

                if margin >= self.ambiguity_margin {
                    TypingDecision::Assign {
                        best_space: best.space.clone(),
                        best_identity: best.identity,
                        margin,
                    }
                } else {
                    TypingDecision::Ambiguous {
                        best_space: best.space.clone(),
                        best_identity: best.identity,
                        second_space: second.space.clone(),
                        second_identity: second.identity,
                        margin,
                    }
                }
            }
        }
    }

    /// Classify all reads in a cross-space query index for depletion.
    ///
    /// Returns a hashmap of read_id → DepletionDecision.
    pub fn classify_all_depletion(
        &self,
        index: &CrossSpaceQueryIndex,
    ) -> std::collections::HashMap<String, DepletionDecision> {
        index
            .iter()
            .map(|(read_id, placements)| {
                (read_id.clone(), self.classify_depletion(placements))
            })
            .collect()
    }

    /// Classify all reads in a cross-space query index for typing.
    ///
    /// Returns a hashmap of read_id → TypingDecision.
    pub fn classify_all_typing(
        &self,
        index: &CrossSpaceQueryIndex,
    ) -> std::collections::HashMap<String, TypingDecision> {
        index
            .iter()
            .map(|(read_id, placements)| {
                (read_id.clone(), self.classify_typing(placements))
            })
            .collect()
    }
}

impl Default for &CrossRefFilter {
    fn default() -> Self {
        &CrossRefFilter {
            ambiguity_margin: 0.02,
            min_identity: 0.95,
            deplete_space: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use phraya_io::queries::CrossSpacePlacement;

    fn placement(space: &str, identity: f64) -> CrossSpacePlacement {
        CrossSpacePlacement {
            space: space.to_string(),
            pos: 100,
            identity,
        }
    }

    fn placements(specs: &[(&str, f64)]) -> Vec<CrossSpacePlacement> {
        specs
            .iter()
            .map(|(s, i)| placement(s, *i))
            .collect()
    }

    // ---- Depletion classification tests ----

    #[test]
    fn test_depletion_confident_target_hit() {
        let filter = CrossRefFilter::new().with_deplete_space("host");
        // Read hits host at 0.99, target at 0.95 — both pass threshold
        let placements = placements(&[("host", 0.99), ("target", 0.95)]);
        let decision = filter.classify_depletion(&placements);

        assert_eq!(
            decision,
            DepletionDecision::Drop {
                reason: DropReason::ConfidentHit
            }
        );
    }

    #[test]
    fn test_depletion_ambiguous_near_tie() {
        let filter = CrossRefFilter::new().with_deplete_space("host");
        // Read hits host at 0.96, target at 0.95 — margin < 0.02
        let placements = placements(&[("host", 0.96), ("target", 0.95)]);
        let decision = filter.classify_depletion(&placements);

        match &decision {
            DepletionDecision::Keep {
                reason: KeepReason::Ambiguous { margin },
            } => assert!(*margin < 0.02),
            other => panic!("expected ambiguous keep, got {:?}", other),
        }
    }

    #[test]
    fn test_depletion_better_elsewhere() {
        let filter = CrossRefFilter::new().with_deplete_space("host");
        // Read hits target at 0.99, host at 0.80 — below min_identity
        let placements = placements(&[("target", 0.99), ("host", 0.80)]);
        let decision = filter.classify_depletion(&placements);

        assert_eq!(
            decision,
            DepletionDecision::Keep {
                reason: KeepReason::BetterElsewhere
            }
        );
    }

    #[test]
    fn test_depletion_sole_hit() {
        let filter = CrossRefFilter::new().with_deplete_space("host");
        // Only one placement, in the target space
        let placements = placements(&[("host", 0.97)]);
        let decision = filter.classify_depletion(&placements);

        assert_eq!(
            decision,
            DepletionDecision::Drop {
                reason: DropReason::SoleHit
            }
        );
    }

    #[test]
    fn test_depletion_no_target_space_match() {
        let filter = CrossRefFilter::new().with_deplete_space("host");
        // Read only has placements outside the target space
        let placements = placements(&[("target", 0.99), ("vector", 0.97)]);
        let decision = filter.classify_depletion(&placements);

        assert_eq!(
            decision,
            DepletionDecision::Keep {
                reason: KeepReason::BetterElsewhere
            }
        );
    }

    // ---- Typing classification tests ----

    #[test]
    fn test_typing_clear_winner() {
        let filter = CrossRefFilter::new();
        let placements = placements(&[("host", 0.99), ("target", 0.80)]);
        let decision = filter.classify_typing(&placements);

        match decision {
            TypingDecision::Assign {
                best_space,
                margin,
                ..
            } => {
                assert_eq!(best_space, "host");
                assert!(margin >= 0.02);
            }
            other => panic!("expected assign, got {:?}", other),
        }
    }

    #[test]
    fn test_typing_ambiguous() {
        let filter = CrossRefFilter::new();
        let placements = placements(&[("host", 0.96), ("target", 0.95)]);
        let decision = filter.classify_typing(&placements);

        match &decision {
            TypingDecision::Ambiguous {
                best_space,
                second_space,
                margin,
                ..
            } => {
                assert_eq!(best_space, "host");
                assert_eq!(second_space, "target");
                assert!(*margin < 0.02);
            }
            other => panic!("expected ambiguous, got {:?}", other),
        }
    }

    #[test]
    fn test_typing_unassigned() {
        let filter = CrossRefFilter::new();
        // All placements below min_identity (0.95)
        let placements = placements(&[("host", 0.85), ("target", 0.80)]);
        let decision = filter.classify_typing(&placements);

        assert_eq!(decision, TypingDecision::Unassigned);
    }

    #[test]
    fn test_default_ambiguity_margin() {
        let filter = CrossRefFilter::new();
        assert_eq!(filter.ambiguity_margin, 0.02);
        assert_eq!(filter.min_identity, 0.95);
        assert!(filter.deplete_space.is_none());
    }

    #[test]
    fn test_custom_ambiguity_margin() {
        let filter = CrossRefFilter::new()
            .with_deplete_space("host")
            .with_ambiguity_margin(0.05);
        assert_eq!(filter.ambiguity_margin, 0.05);

        // Now 0.99 vs 0.96 has margin 0.03, which is < 0.05 → ambiguous
        let placements = placements(&[("host", 0.99), ("target", 0.96)]);
        let decision = filter.classify_depletion(&placements);
        match &decision {
            DepletionDecision::Keep {
                reason: KeepReason::Ambiguous { .. },
            } => {}
            other => panic!("expected ambiguous keep with margin 0.05, got {:?}", other),
        }
    }
}
