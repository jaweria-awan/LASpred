use acetyl_scan::{
    amino, crossval, fasta,
    model::{self, Profile, Site, SkipReason, Window},
    tables::{LEGACY_THRESHOLD, PROPERTY_SCALES},
    BUNDLED_TRAINING,
};
use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(version, about = "Score lysine (K) sites for potential acetylation", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Score every lysine in a protein FASTA file (TSV output).
    Scan(ScanArgs),
    /// Print the position-specific residue frequency profile (TSV).
    Profile(ProfileArgs),
    /// Score the training sites with k-fold held-out profiles (TSV).
    Crossval(CrossvalArgs),
}

#[derive(Args)]
struct TrainingArg {
    /// Training windows (21 residues, K at centre, one per line).
    /// Defaults to the bundled 2,000 human positive sites.
    #[arg(long, value_name = "FILE")]
    training: Option<PathBuf>,
}

#[derive(Args)]
struct ScanArgs {
    /// Protein FASTA file, or '-' for stdin.
    input: PathBuf,
    /// Output file (default: stdout).
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,
    #[command(flatten)]
    training: TrainingArg,
    /// Score cutoff for the `predicted` column. The default is the original
    /// program's value, which nearly every lysine exceeds.
    #[arg(long, default_value_t = LEGACY_THRESHOLD)]
    threshold: f64,
    /// Write only sites with total > threshold.
    #[arg(long)]
    only_predicted: bool,
    /// Also write rows for lysines that could not be scored.
    #[arg(long)]
    include_skipped: bool,
    /// Decimal places in score columns.
    #[arg(long, default_value_t = 4)]
    precision: usize,
}

#[derive(Args)]
struct ProfileArgs {
    #[command(flatten)]
    training: TrainingArg,
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,
}

#[derive(Args)]
struct CrossvalArgs {
    #[command(flatten)]
    training: TrainingArg,
    /// Number of folds.
    #[arg(short = 'k', long, default_value_t = 5)]
    folds: usize,
    /// Seed for fold assignment.
    #[arg(long, default_value_t = 42)]
    seed: u64,
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,
    #[arg(long, default_value_t = 4)]
    precision: usize,
}

/// Number of blank fields between `position` and `status` in a skipped row.
const SCAN_EMPTY_FIELDS: usize = 1 + 7 + 5;

fn open_output(path: &Option<PathBuf>) -> Result<Box<dyn Write>> {
    Ok(match path {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).with_context(|| format!("creating {}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    })
}

fn load_training(arg: &TrainingArg) -> Result<Vec<Window>> {
    match &arg.training {
        Some(p) => {
            let f = File::open(p).with_context(|| format!("opening {}", p.display()))?;
            model::read_training(BufReader::new(f), &p.display().to_string())
        }
        None => model::read_training(BUNDLED_TRAINING.as_bytes(), "bundled training set"),
    }
}

fn score_header() -> String {
    let props: Vec<&str> = PROPERTY_SCALES.iter().map(|s| s.name).collect();
    format!("{}\tpssm\tdimer\tmissing_dimers\ttotal", props.join("\t"))
}

fn write_score(out: &mut dyn Write, s: &model::Score, p: usize) -> io::Result<()> {
    for v in s.properties {
        write!(out, "{v:.p$}\t")?;
    }
    write!(
        out,
        "{:.p$}\t{:.p$}\t{}\t{:.p$}",
        s.pssm, s.dimer, s.missing_dimers, s.total
    )
}

fn cmd_scan(a: ScanArgs) -> Result<()> {
    let windows = load_training(&a.training)?;
    let profile = Profile::from_windows(&windows)?;
    let input: Box<dyn BufRead> = if a.input.as_os_str() == "-" {
        Box::new(BufReader::new(io::stdin().lock()))
    } else {
        let f = File::open(&a.input).with_context(|| format!("opening {}", a.input.display()))?;
        Box::new(BufReader::new(f))
    };
    let mut out = open_output(&a.output)?;
    let p = a.precision;
    writeln!(
        out,
        "protein_id\tposition\twindow\t{}\tpredicted\tstatus",
        score_header()
    )?;

    let (mut proteins, mut lysines, mut scored, mut predicted) = (0usize, 0usize, 0usize, 0usize);
    let mut skipped: BTreeMap<SkipReason, usize> = BTreeMap::new();
    let mut ids_seen: BTreeMap<String, usize> = BTreeMap::new();
    let mut empty_ids = 0usize;
    for rec in fasta::Reader::new(input) {
        let rec = rec.with_context(|| format!("parsing {}", a.input.display()))?;
        proteins += 1;
        let id = match rec.id() {
            "" => {
                empty_ids += 1;
                format!("record_{proteins}")
            }
            s => s.to_string(),
        };
        *ids_seen.entry(id.clone()).or_default() += 1;
        for site in model::scan(&rec.sequence, &profile) {
            lysines += 1;
            match site {
                Site::Scored {
                    position,
                    window,
                    score,
                } => {
                    scored += 1;
                    let hit = score.total > a.threshold;
                    predicted += hit as usize;
                    if a.only_predicted && !hit {
                        continue;
                    }
                    write!(out, "{id}\t{position}\t{}\t", window.to_string_lossy())?;
                    write_score(&mut *out, &score, p)?;
                    writeln!(out, "\t{}\tscored", hit as u8)?;
                }
                Site::Skipped { position, reason } => {
                    *skipped.entry(reason).or_default() += 1;
                    if a.include_skipped && !a.only_predicted {
                        // window + 7 properties + pssm, dimer, missing_dimers, total, predicted
                        let blanks = "\t".repeat(SCAN_EMPTY_FIELDS);
                        writeln!(out, "{id}\t{position}{blanks}\t{}", reason.as_str())?;
                    }
                }
            }
        }
    }
    out.flush()?;

    let skipped_str: Vec<String> = skipped
        .iter()
        .map(|(r, n)| format!("{}={n}", r.as_str()))
        .collect();
    eprintln!(
        "proteins={proteins} lysines={lysines} scored={scored} predicted={predicted} skipped: {}",
        if skipped_str.is_empty() {
            "none".into()
        } else {
            skipped_str.join(" ")
        }
    );
    let dupes = ids_seen.values().filter(|&&n| n > 1).count();
    if dupes > 0 {
        eprintln!("warning: {dupes} protein id(s) occur more than once; rows are ambiguous");
    }
    if empty_ids > 0 {
        eprintln!("warning: {empty_ids} record(s) had an empty header; named record_<n>");
    }
    if a.threshold == LEGACY_THRESHOLD && scored > 0 {
        eprintln!(
            "note: threshold {LEGACY_THRESHOLD} is the original program's value and was not \
             calibrated for 21-residue windows; it marks nearly every lysine. Use the `total` \
             column or set --threshold."
        );
    }
    Ok(())
}

fn cmd_profile(a: ProfileArgs) -> Result<()> {
    let windows = load_training(&a.training)?;
    let profile = Profile::from_windows(&windows)?;
    let mut out = open_output(&a.output)?;
    write!(out, "residue")?;
    for pos in 0..model::WINDOW_LEN {
        write!(out, "\t{}", pos as i64 - model::FLANK as i64)?;
    }
    writeln!(out)?;
    for (aa_i, &aa) in amino::ALPHABET.iter().enumerate() {
        write!(out, "{}", aa as char)?;
        for pos in 0..model::WINDOW_LEN {
            write!(out, "\t{:.4}", profile.freq[pos][aa_i])?;
        }
        writeln!(out)?;
    }
    out.flush()?;
    eprintln!("sites={}", profile.n_sites);
    Ok(())
}

fn mean_sd(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0).max(1.0);
    (m, var.sqrt())
}

fn cmd_crossval(a: CrossvalArgs) -> Result<()> {
    let windows = load_training(&a.training)?;
    let rows = crossval::run(&windows, a.folds, a.seed)?;
    let mut out = open_output(&a.output)?;
    let p = a.precision;
    writeln!(
        out,
        "index\tfold\twindow\t{}\tin_sample_total",
        score_header()
    )?;
    for r in &rows {
        write!(
            out,
            "{}\t{}\t{}\t",
            r.index + 1,
            r.fold,
            windows[r.index].to_string_lossy()
        )?;
        write_score(&mut *out, &r.held_out, p)?;
        writeln!(out, "\t{:.p$}", r.in_sample_total)?;
    }
    out.flush()?;
    let held: Vec<f64> = rows.iter().map(|r| r.held_out.total).collect();
    let ins: Vec<f64> = rows.iter().map(|r| r.in_sample_total).collect();
    let (hm, hs) = mean_sd(&held);
    let (im, is) = mean_sd(&ins);
    eprintln!(
        "sites={} folds={} held_out_total mean={hm:.3} sd={hs:.3} | in_sample_total mean={im:.3} sd={is:.3}",
        rows.len(),
        a.folds
    );
    Ok(())
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Scan(a) => cmd_scan(a),
        Command::Profile(a) => cmd_profile(a),
        Command::Crossval(a) => cmd_crossval(a),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            // Exiting quietly on a closed pipe (e.g. `| head`) is expected behaviour.
            if e.chain().any(|c| {
                c.downcast_ref::<io::Error>()
                    .is_some_and(|io| io.kind() == io::ErrorKind::BrokenPipe)
            }) {
                return ExitCode::SUCCESS;
            }
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}
