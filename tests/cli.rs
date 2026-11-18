use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run(args: &[&str], stdin: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_acetyl-scan"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(stdin.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

fn stdout(o: &Output) -> String {
    String::from_utf8(o.stdout.clone()).unwrap()
}

const FASTA: &str = ">p1 test\nAAAAAAAAAAKAAAAAAAAAAAXAAAAAAAAAKAAAAAAAAAA\n>p2\nKAAA\n";

#[test]
fn scan_reports_scored_and_skipped_with_consistent_columns() {
    let o = run(&["scan", "-", "--include-skipped"], FASTA);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let out = stdout(&o);
    let lines: Vec<&str> = out.lines().collect();
    let ncol = lines[0].split('\t').count();
    assert_eq!(ncol, 16);
    for l in &lines {
        assert_eq!(l.split('\t').count(), ncol, "ragged row: {l:?}");
    }
    let status: Vec<&str> = lines[1..]
        .iter()
        .map(|l| l.rsplit('\t').next().unwrap())
        .collect();
    assert_eq!(status, ["scored", "non_standard_residue", "near_terminus"]);
    assert!(lines[1].starts_with("p1\t11\tAAAAAAAAAAKAAAAAAAAAA\t"));
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("proteins=2 lysines=3 scored=1"), "{err}");
}

#[test]
fn threshold_and_only_predicted_filter_rows() {
    let o = run(
        &["scan", "-", "--threshold", "1e9", "--only-predicted"],
        FASTA,
    );
    assert!(o.status.success());
    assert_eq!(
        stdout(&o).lines().count(),
        1,
        "only the header should remain"
    );
}

#[test]
fn non_fasta_input_fails_cleanly() {
    let o = run(&["scan", "-"], "MKVLAAGK\n");
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("before the first '>' header"));
}

#[test]
fn bad_training_file_is_rejected_with_location() {
    let dir = std::env::temp_dir().join(format!("acetyl-scan-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("train.txt");
    std::fs::write(&path, "AAAAAAAAAAKAAAAAAAAAA\nAAAAAAAAAAAAAAAAAAAAA\n").unwrap();
    let o = run(&["scan", "-", "--training", path.to_str().unwrap()], FASTA);
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("train.txt:2"));
}

#[test]
fn profile_and_crossval_run_on_bundled_data() {
    let p = run(&["profile"], "");
    assert!(p.status.success());
    let out = stdout(&p);
    assert_eq!(out.lines().count(), 21); // header + 20 residues
    let k_row = out.lines().find(|l| l.starts_with("K\t")).unwrap();
    assert_eq!(k_row.split('\t').nth(11), Some("1.0000")); // centre column is always K

    let c = run(&["crossval", "-k", "5"], "");
    assert!(c.status.success());
    assert_eq!(stdout(&c).lines().count(), 2001);
    assert!(String::from_utf8_lossy(&c.stderr).contains("sites=2000 folds=5"));
}
