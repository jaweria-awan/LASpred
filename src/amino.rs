//! The 20 standard amino acids and their table indices.

/// Canonical alphabet; every table in [`crate::tables`] is indexed in this order.
pub const ALPHABET: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";

/// Index of an amino acid in [`ALPHABET`], case-insensitive.
/// Returns `None` for anything that is not one of the 20 standard residues
/// (e.g. X, U, B, Z, O, gaps).
#[inline]
pub fn index(residue: u8) -> Option<usize> {
    match residue.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'D' => Some(2),
        b'E' => Some(3),
        b'F' => Some(4),
        b'G' => Some(5),
        b'H' => Some(6),
        b'I' => Some(7),
        b'K' => Some(8),
        b'L' => Some(9),
        b'M' => Some(10),
        b'N' => Some(11),
        b'P' => Some(12),
        b'Q' => Some(13),
        b'R' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'V' => Some(17),
        b'W' => Some(18),
        b'Y' => Some(19),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_round_trips_alphabet() {
        for (i, &aa) in ALPHABET.iter().enumerate() {
            assert_eq!(index(aa), Some(i));
            assert_eq!(index(aa.to_ascii_lowercase()), Some(i));
        }
    }

    #[test]
    fn non_standard_residues_have_no_index() {
        for c in b"XUBZOJ*-. 1" {
            assert_eq!(index(*c), None, "{}", *c as char);
        }
    }
}
