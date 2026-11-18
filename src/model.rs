//! Position profile, window extraction and the integrated score.

use crate::amino;
use crate::tables::{DIMER, DIMER_COEFFICIENT, PROPERTY_SCALES, PSSM_COEFFICIENT};
use anyhow::{bail, Context, Result};
use std::io::BufRead;

/// Residues on each side of the central lysine.
pub const FLANK: usize = 10;
/// Total window length (10 + K + 10).
pub const WINDOW_LEN: usize = 2 * FLANK + 1;

/// A validated 21-residue window, stored as alphabet indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Window([usize; WINDOW_LEN]);

impl Window {
    /// Builds a window from exactly 21 standard residues.
    /// On failure returns the offending byte (or `None` for a wrong length).
    pub fn from_bytes(bytes: &[u8]) -> std::result::Result<Self, Option<u8>> {
        if bytes.len() != WINDOW_LEN {
            return Err(None);
        }
        let mut idx = [0usize; WINDOW_LEN];
        for (slot, &b) in idx.iter_mut().zip(bytes) {
            *slot = amino::index(b).ok_or(Some(b))?;
        }
        Ok(Self(idx))
    }

    pub fn indices(&self) -> &[usize; WINDOW_LEN] {
        &self.0
    }

    pub fn to_string_lossy(&self) -> String {
        self.0.iter().map(|&i| amino::ALPHABET[i] as char).collect()
    }
}

/// Position-specific residue frequencies estimated from known sites.
#[derive(Debug, Clone)]
pub struct Profile {
    /// `freq[position][amino_acid]`.
    pub freq: [[f64; 20]; WINDOW_LEN],
    pub n_sites: usize,
}

impl Profile {
    pub fn from_windows<'a>(windows: impl IntoIterator<Item = &'a Window>) -> Result<Self> {
        let mut counts = [[0usize; 20]; WINDOW_LEN];
        let mut n = 0usize;
        for w in windows {
            n += 1;
            for (pos, &aa) in w.indices().iter().enumerate() {
                counts[pos][aa] += 1;
            }
        }
        if n == 0 {
            bail!("training set is empty");
        }
        let mut freq = [[0.0; 20]; WINDOW_LEN];
        for (f_row, c_row) in freq.iter_mut().zip(&counts) {
            for (f, &c) in f_row.iter_mut().zip(c_row) {
                *f = c as f64 / n as f64;
            }
        }
        Ok(Self { freq, n_sites: n })
    }
}

/// Reads and validates a training file: one 21-residue window per line with
/// K at the centre. Blank lines and `>` comment lines are ignored.
pub fn read_training(reader: impl BufRead, source: &str) -> Result<Vec<Window>> {
    let mut windows = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line_no = i + 1;
        let line = line.with_context(|| format!("{source}: reading line {line_no}"))?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('>') {
            continue;
        }
        let w = match Window::from_bytes(s.as_bytes()) {
            Ok(w) => w,
            Err(None) => bail!(
                "{source}:{line_no}: expected {WINDOW_LEN} residues, found {}",
                s.len()
            ),
            Err(Some(b)) => bail!(
                "{source}:{line_no}: non-standard residue {:?} in {s:?}",
                b as char
            ),
        };
        if w.indices()[FLANK] != amino::index(b'K').unwrap() {
            bail!(
                "{source}:{line_no}: centre residue (position {}) is not K in {s:?}",
                FLANK + 1
            );
        }
        windows.push(w);
    }
    if windows.is_empty() {
        bail!("{source}: no training windows found");
    }
    Ok(windows)
}

/// Feature breakdown for one window.
#[derive(Debug, Clone, PartialEq)]
pub struct Score {
    /// Weighted property contributions, in [`PROPERTY_SCALES`] order.
    pub properties: [f64; 7],
    pub pssm: f64,
    pub dimer: f64,
    /// Dipeptides in the window that have no value in the original table.
    pub missing_dimers: usize,
    pub total: f64,
}

pub fn score(window: &Window, profile: &Profile) -> Score {
    let idx = window.indices();
    let mut properties = [0.0; 7];
    for (out, scale) in properties.iter_mut().zip(PROPERTY_SCALES) {
        *out = idx.iter().map(|&a| scale.values[a]).sum::<f64>() * scale.coefficient;
    }
    let pssm = idx
        .iter()
        .enumerate()
        .map(|(pos, &a)| profile.freq[pos][a])
        .sum::<f64>()
        * PSSM_COEFFICIENT;
    let mut dimer_sum = 0.0;
    let mut missing_dimers = 0;
    for pair in idx.windows(2) {
        match DIMER[pair[0]][pair[1]] {
            Some(v) => dimer_sum += v,
            None => missing_dimers += 1,
        }
    }
    let dimer = dimer_sum * DIMER_COEFFICIENT;
    let total = properties.iter().sum::<f64>() + pssm + dimer;
    Score {
        properties,
        pssm,
        dimer,
        missing_dimers,
        total,
    }
}

/// Why a lysine in a query sequence could not be scored.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SkipReason {
    /// Fewer than 10 residues on one side.
    NearTerminus,
    /// The window contains a residue outside the 20 standard amino acids.
    NonStandardResidue,
}

impl SkipReason {
    pub fn as_str(self) -> &'static str {
        match self {
            SkipReason::NearTerminus => "near_terminus",
            SkipReason::NonStandardResidue => "non_standard_residue",
        }
    }
}

/// Result for one lysine in a query sequence.
#[derive(Debug, Clone)]
// Sites are yielded one at a time and consumed immediately; boxing would only
// add an allocation per lysine.
#[allow(clippy::large_enum_variant)]
pub enum Site {
    Scored {
        position: usize,
        window: Window,
        score: Score,
    },
    Skipped {
        position: usize,
        reason: SkipReason,
    },
}

/// Finds every K in `sequence` (expected uppercase) and scores those with a
/// full window. `position` is 1-based.
pub fn scan<'a>(sequence: &'a [u8], profile: &'a Profile) -> impl Iterator<Item = Site> + 'a {
    sequence
        .iter()
        .enumerate()
        .filter(|(_, &b)| b == b'K')
        .map(move |(i, _)| {
            let position = i + 1;
            if i < FLANK || i + FLANK >= sequence.len() {
                return Site::Skipped {
                    position,
                    reason: SkipReason::NearTerminus,
                };
            }
            match Window::from_bytes(&sequence[i - FLANK..=i + FLANK]) {
                Ok(window) => Site::Scored {
                    position,
                    window,
                    score: score(&window, profile),
                },
                Err(_) => Site::Skipped {
                    position,
                    reason: SkipReason::NonStandardResidue,
                },
            }
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform_profile() -> Profile {
        Profile {
            freq: [[0.05; 20]; WINDOW_LEN],
            n_sites: 1,
        }
    }

    #[test]
    fn edges_are_inclusive() {
        // K at 0-based 10 and at len-11: both have exactly 10 flanking residues.
        let mut seq = vec![b'A'; 41];
        seq[10] = b'K';
        seq[30] = b'K';
        let p = uniform_profile();
        let sites: Vec<_> = scan(&seq, &p).collect();
        assert_eq!(sites.len(), 2);
        assert!(sites.iter().all(|s| matches!(s, Site::Scored { .. })));
        assert!(matches!(sites[0], Site::Scored { position: 11, .. }));
        assert!(matches!(sites[1], Site::Scored { position: 31, .. }));
    }

    #[test]
    fn near_terminus_and_non_standard_are_reported() {
        let p = uniform_profile();
        let seq = b"AKAAAAAAAAAAXAAAAAAAAAKAAAAAAAAAAAK";
        let reasons: Vec<_> = scan(seq, &p)
            .map(|s| match s {
                Site::Skipped { reason, .. } => Some(reason),
                Site::Scored { .. } => None,
            })
            .collect();
        assert_eq!(
            reasons,
            vec![
                Some(SkipReason::NearTerminus),
                Some(SkipReason::NonStandardResidue),
                Some(SkipReason::NearTerminus)
            ]
        );
    }

    #[test]
    fn missing_dimers_are_counted() {
        let mut w = *b"AAAAAAAAAAKAAAAAAAAAA";
        w[0] = b'C';
        w[1] = b'W'; // CW
        w[2] = b'M'; // WM
        w[3] = b'W'; // MW
        let s = score(&Window::from_bytes(&w).unwrap(), &uniform_profile());
        assert_eq!(s.missing_dimers, 3);
    }

    #[test]
    fn total_is_sum_of_parts() {
        let w = Window::from_bytes(b"RLNACFNLQQKKGTREKGERD").unwrap();
        let s = score(&w, &uniform_profile());
        let parts: f64 = s.properties.iter().sum::<f64>() + s.pssm + s.dimer;
        assert!((parts - s.total).abs() < 1e-9);
    }

    #[test]
    fn training_validation_errors_name_the_line() {
        let bad_len = read_training("AAAAAAAAAAKAAAAAAAAAA\nAAAK\n".as_bytes(), "t").unwrap_err();
        assert!(bad_len.to_string().contains("t:2"), "{bad_len}");
        let bad_centre = read_training("AAAAAAAAAAAAAAAAAAAAA\n".as_bytes(), "t").unwrap_err();
        assert!(bad_centre.to_string().contains("not K"), "{bad_centre}");
        let bad_res = read_training("AAAAAAAAAAKAAAAAAAAAX\n".as_bytes(), "t").unwrap_err();
        assert!(bad_res.to_string().contains("'X'"), "{bad_res}");
    }

    #[test]
    fn training_tolerates_crlf_blank_and_comment_lines() {
        let w = read_training(">c\r\n\r\nAAAAAAAAAAKAAAAAAAAAA\r\n".as_bytes(), "t").unwrap();
        assert_eq!(w.len(), 1);
    }
}
