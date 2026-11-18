//! k-fold scoring of the training sites.
//!
//! Scoring the training windows with a profile built from those same windows
//! inflates the PSSM term. Here every window is scored with a profile built
//! only from the other folds, so its score is a fair estimate of how a
//! comparable unseen site would score.

use crate::model::{score, Profile, Score, Window};
use anyhow::{bail, Result};

pub struct CvRow {
    /// 0-based index of the window in the training input.
    pub index: usize,
    pub fold: usize,
    /// Score using a profile that excludes this window's fold.
    pub held_out: Score,
    /// Score using the profile built from all windows (the circular estimate).
    pub in_sample_total: f64,
}

/// SplitMix64: small, dependency-free, deterministic PRNG for fold assignment.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Assigns each of `n` items to one of `k` folds of near-equal size.
pub fn assign_folds(n: usize, k: usize, seed: u64) -> Vec<usize> {
    let mut order: Vec<usize> = (0..n).collect();
    let mut state = seed;
    for i in (1..n).rev() {
        let j = (splitmix64(&mut state) % (i as u64 + 1)) as usize;
        order.swap(i, j);
    }
    let mut folds = vec![0; n];
    for (rank, &item) in order.iter().enumerate() {
        folds[item] = rank % k;
    }
    folds
}

pub fn run(windows: &[Window], k: usize, seed: u64) -> Result<Vec<CvRow>> {
    if k < 2 {
        bail!("need at least 2 folds, got {k}");
    }
    if windows.len() < k {
        bail!(
            "{} training windows cannot be split into {k} folds",
            windows.len()
        );
    }
    let folds = assign_folds(windows.len(), k, seed);
    let full = Profile::from_windows(windows)?;
    let mut rows = Vec::with_capacity(windows.len());
    for fold in 0..k {
        let train = windows
            .iter()
            .zip(&folds)
            .filter(|(_, &f)| f != fold)
            .map(|(w, _)| w);
        let profile = Profile::from_windows(train)?;
        for (index, (w, _)) in windows
            .iter()
            .zip(&folds)
            .enumerate()
            .filter(|(_, (_, &f))| f == fold)
        {
            rows.push(CvRow {
                index,
                fold,
                held_out: score(w, &profile),
                in_sample_total: score(w, &full).total,
            });
        }
    }
    rows.sort_by_key(|r| r.index);
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn folds_are_balanced_and_deterministic() {
        let a = assign_folds(103, 5, 7);
        assert_eq!(a, assign_folds(103, 5, 7));
        let mut sizes = [0; 5];
        a.iter().for_each(|&f| sizes[f] += 1);
        assert!(sizes.iter().all(|&s| s == 20 || s == 21), "{sizes:?}");
    }

    #[test]
    fn rejects_bad_k() {
        let w = [Window::from_bytes(b"AAAAAAAAAAKAAAAAAAAAA").unwrap()];
        assert!(run(&w, 1, 0).is_err());
        assert!(run(&w, 2, 0).is_err());
    }
}
