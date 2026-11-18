//! Lysine acetylation site scoring.
//!
//! Rust port of a 2011 Java tool. Each lysine with 10 residues on either side
//! is scored as a linear combination of nine features: seven amino-acid
//! property scales, a position-specific frequency profile built from known
//! sites, and a dipeptide table.

pub mod amino;
pub mod crossval;
pub mod fasta;
pub mod model;
pub mod tables;

/// Training windows bundled with the binary (2,000 human positive sites).
pub const BUNDLED_TRAINING: &str = include_str!("../data/human_positive_sites.txt");
