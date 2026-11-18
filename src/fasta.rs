//! Minimal, strict-but-tolerant streaming FASTA reader.
//!
//! Tolerates CRLF line endings, blank lines, lowercase residues, internal
//! whitespace and a trailing `*` stop symbol. Sequence data appearing before
//! the first `>` header is an error (the original Java looped forever on it).

use anyhow::{bail, Context, Result};
use std::io::BufRead;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// Full header line without the leading `>`.
    pub header: String,
    /// Uppercased residues with whitespace and a terminal `*` removed.
    pub sequence: Vec<u8>,
}

impl Record {
    /// First whitespace-delimited token of the header.
    pub fn id(&self) -> &str {
        self.header.split_whitespace().next().unwrap_or("")
    }
}

/// Iterator over FASTA records read from `reader`.
pub struct Reader<R> {
    lines: std::io::Lines<R>,
    pending_header: Option<String>,
    line_no: usize,
    started: bool,
}

impl<R: BufRead> Reader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            lines: reader.lines(),
            pending_header: None,
            line_no: 0,
            started: false,
        }
    }

    fn next_record(&mut self) -> Result<Option<Record>> {
        let mut header = self.pending_header.take();
        let mut sequence = Vec::new();
        while let Some(line) = self.lines.next() {
            self.line_no += 1;
            let line = line.with_context(|| format!("reading line {}", self.line_no))?;
            let line = line.trim_end_matches(['\r', '\n']);
            if let Some(h) = line.strip_prefix('>') {
                let h = h.trim().to_string();
                if !self.started {
                    self.started = true;
                    header = Some(h);
                    continue;
                }
                self.pending_header = Some(h);
                break;
            }
            if line.trim().is_empty() {
                continue;
            }
            if !self.started {
                bail!(
                    "line {}: sequence data before the first '>' header; input is not FASTA",
                    self.line_no
                );
            }
            sequence.extend(
                line.bytes()
                    .filter(|b| !b.is_ascii_whitespace())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
        match header {
            None => Ok(None),
            Some(header) => {
                if sequence.last() == Some(&b'*') {
                    sequence.pop();
                }
                Ok(Some(Record { header, sequence }))
            }
        }
    }
}

impl<R: BufRead> Iterator for Reader<R> {
    type Item = Result<Record>;
    fn next(&mut self) -> Option<Self::Item> {
        self.next_record().transpose()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> Result<Vec<Record>> {
        Reader::new(s.as_bytes()).collect()
    }

    #[test]
    fn parses_multiline_crlf_lowercase_and_stop() {
        let r = parse(">p1 desc\r\nmkv\r\n\r\nAAk*\r\n>p2\nKK\n").unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].id(), "p1");
        assert_eq!(r[0].header, "p1 desc");
        assert_eq!(r[0].sequence, b"MKVAAK");
        assert_eq!(r[1].sequence, b"KK");
    }

    #[test]
    fn empty_record_and_trailing_header() {
        let r = parse(">a\n>b\nAC\n>c\n").unwrap();
        assert_eq!(
            r.iter().map(|r| r.sequence.len()).collect::<Vec<_>>(),
            vec![0, 2, 0]
        );
    }

    #[test]
    fn leading_blank_lines_are_ok() {
        assert_eq!(parse("\n\n>a\nK\n").unwrap().len(), 1);
    }

    #[test]
    fn sequence_before_header_is_an_error() {
        let err = parse("MKV\n>a\nK\n").unwrap_err().to_string();
        assert!(err.contains("line 1"), "{err}");
    }

    #[test]
    fn empty_input_yields_nothing() {
        assert!(parse("").unwrap().is_empty());
    }
}
