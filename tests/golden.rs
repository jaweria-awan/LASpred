//! Checks the port against output captured from the original Java program.
//!
//! Expected differences, both intentional bug fixes:
//! * The Java code rounded each feature and profile frequency to 4 decimals;
//!   this port does not, so totals agree to within a small tolerance.
//! * The Java code skipped the first and last lysine that had a full window
//!   (off-by-one); this port scores them, so it reports strictly more sites.

use acetyl_scan::{fasta, model, BUNDLED_TRAINING};
use std::collections::{HashMap, HashSet};

// Worst case of the Java rounding: 21 frequencies * 5e-5 * 17.66 + 9 features * 5e-5.
const TOLERANCE: f64 = 0.02;

#[test]
fn matches_original_java_output() {
    let windows = model::read_training(BUNDLED_TRAINING.as_bytes(), "bundled").unwrap();
    let profile = model::Profile::from_windows(&windows).unwrap();

    let fasta_text = include_str!("fixtures/reference.fasta");
    let records: Vec<_> = fasta::Reader::new(fasta_text.as_bytes())
        .map(|r| r.unwrap())
        .collect();

    // (protein, 1-based position) -> (window, total, sequence length)
    let mut rust: HashMap<(String, usize), (String, f64, usize)> = HashMap::new();
    for rec in &records {
        for site in model::scan(&rec.sequence, &profile) {
            if let model::Site::Scored {
                position,
                window,
                score,
            } = site
            {
                rust.insert(
                    (rec.id().to_string(), position),
                    (window.to_string_lossy(), score.total, rec.sequence.len()),
                );
            }
        }
    }

    let mut java_keys = HashSet::new();
    let mut max_diff: f64 = 0.0;
    for line in include_str!("fixtures/java_reference.tsv").lines().skip(3) {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let key = (f[0].to_string(), f[1].parse::<usize>().unwrap() + 1);
        let java_total: f64 = f[3].parse().unwrap();
        let (window, total, _) = rust
            .get(&key)
            .unwrap_or_else(|| panic!("Rust did not score Java site {key:?}"));
        assert_eq!(window, f[2], "window mismatch at {key:?}");
        let d = (total - java_total).abs();
        assert!(d < TOLERANCE, "{key:?}: rust {total} vs java {java_total}");
        max_diff = max_diff.max(d);
        java_keys.insert(key);
    }
    assert!(java_keys.len() > 200, "fixture looks truncated");
    eprintln!(
        "max |rust - java| = {max_diff:.6} over {} sites",
        java_keys.len()
    );

    // Every extra Rust site must be an edge lysine the Java code missed.
    let extras: Vec<_> = rust
        .iter()
        .filter(|(k, _)| !java_keys.contains(*k))
        .collect();
    assert!(
        !extras.is_empty(),
        "fixture contains edge cases; expected extra sites"
    );
    for ((id, pos), (_, _, len)) in extras {
        let i = pos - 1;
        assert!(
            i == model::FLANK || i + model::FLANK + 1 == *len,
            "{id}:{pos} scored by Rust but not Java and is not an edge position"
        );
    }
}
