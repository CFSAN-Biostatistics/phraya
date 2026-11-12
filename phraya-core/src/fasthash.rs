//! A fast, non-cryptographic hasher for maps keyed by minimizer values.
//!
//! Minimizer values are already the output of a hash function, and the maps built on them
//! (the per-target seed index, the read-minimizer prefilter, k-mer uniqueness counting) are
//! probed once per minimizer on the alignment hot path. `std`'s default SipHash gives
//! DoS-resistance that sequence data does not need, at several times the per-probe cost.
//!
//! The algorithm is the one `rustc-hash` 2.x uses: add the key, multiply by an odd constant,
//! and rotate in `finish` so the low bits `hashbrown` indexes with come from the
//! well-mixed high bits of the product. Iteration order of these maps differs from the
//! default hasher's, so nothing may depend on it.

use std::hash::{BuildHasherDefault, Hasher};

const K: u64 = 0xf135_7aea_2e62_a9c5;

#[derive(Default, Clone, Copy)]
pub struct MinimizerHasher(u64);

impl MinimizerHasher {
    #[inline]
    fn add(&mut self, v: u64) {
        self.0 = self.0.wrapping_add(v).wrapping_mul(K);
    }
}

impl Hasher for MinimizerHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.add(u64::from(b));
        }
    }
    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(u64::from(i));
    }
    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(u64::from(i));
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.0.rotate_left(26)
    }
}

/// `BuildHasher` for `HashMap`/`HashSet` keyed by minimizer values (or tuples of them).
pub type MinimizerBuildHasher = BuildHasherDefault<MinimizerHasher>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashMap, HashSet};
    use std::hash::{BuildHasher, Hash};

    fn hash_of<T: Hash>(v: &T) -> u64 {
        MinimizerBuildHasher::default().hash_one(v)
    }

    /// Keys that differ only in their high bits must still spread across the low bits
    /// `hashbrown` uses for bucket selection; an identity-style hash would put all 4096 in
    /// one bucket.
    #[test]
    fn high_bit_only_keys_spread_across_low_bits() {
        let buckets: HashSet<u64> = (0..4096u64).map(|i| hash_of(&(i << 32)) & 0xFFF).collect();
        assert!(buckets.len() > 2300, "only {} of 4096 low-12-bit buckets used", buckets.len());
    }

    #[test]
    fn tuple_keys_distinguish_both_components() {
        assert_ne!(hash_of(&(7u64, 1u32)), hash_of(&(7u64, 2u32)));
        assert_ne!(hash_of(&(7u64, 1u32)), hash_of(&(8u64, 1u32)));
    }

    #[test]
    fn works_as_a_map_hasher() {
        let mut m: HashMap<u64, u32, MinimizerBuildHasher> = HashMap::default();
        for i in 0..50_000u64 {
            m.insert(i.wrapping_mul(0x9E37_79B9_7F4A_7C15), i as u32);
        }
        assert_eq!(m.len(), 50_000);
        for i in 0..50_000u64 {
            assert_eq!(m[&i.wrapping_mul(0x9E37_79B9_7F4A_7C15)], i as u32);
        }
    }
}
