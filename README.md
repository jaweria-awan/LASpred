# acetyl-scan

Scores lysine (K) residues in protein sequences for potential acetylation.

This is a Rust port of a Java tool written in 2011 (`FinalCode`). It uses the
same features and coefficients, and it fixes the original's bugs (listed
below). It checks the port against output captured from the original program.

> **Legacy reference.** This repository is a new implementation, not the
> original code. The 2011 Java version is kept unchanged, with steps to run it,
> at **https://github.com/jaweria-awan/LASpred_legacy**. Use it to see where the
> model came from or to reproduce the original results. Use this repository
> for actual work. The differences are listed under
> [Fixes relative to the original Java](#fixes-relative-to-the-original-java).

## Quick start

Only Docker is needed.

```sh
git clone git@github.com:jaweria-awan/LASpred.git && cd LASpred
docker build -t acetyl-scan .              # build the tool
docker build --target test .               # optional: run all tests

cd examples
docker run --rm -v "$PWD":/data acetyl-scan scan test_proteins.fasta > scores.tsv
```

`examples/test_proteins.fasta` contains histone H3.1, histone H4, ubiquitin and
one deliberately messy record (lowercase, an `X`, a blank line, a trailing `*`).
The terminal prints a summary:

```
proteins=4 lysines=42 scored=30 predicted=30 skipped: near_terminus=9 non_standard_residue=3
```

Try next:

```sh
# Also list lysines that could not be scored, with the reason
docker run --rm -v "$PWD":/data acetyl-scan scan test_proteins.fasta --include-skipped

# Only sites scoring above your own cutoff
docker run --rm -v "$PWD":/data acetyl-scan scan test_proteins.fasta --threshold 285 --only-predicted

# Your own proteins: put the file in the current folder
docker run --rm -v "$PWD":/data acetyl-scan scan my_proteins.fasta > my_scores.tsv

# All options
docker run --rm acetyl-scan scan --help
```

### Reading the output

`scores.tsv` is a tab-separated table (opens in Excel), one row per lysine:

| Column | Meaning |
|---|---|
| `protein_id` | First word of the FASTA header |
| `position` | Position of the K in the protein, counting from 1 (including the initial Met) |
| `window` | The 21 residues around it, with the K in the middle |
| 9 feature columns | Each feature's contribution to the score (see [Model](#model)) |
| `missing_dimers` | Pairs in the window with no table value (CW, MW, WM); usually 0 |
| `total` | **The score. Higher means more similar to known acetylation sites** |
| `predicted` | 1 if `total` > threshold, else 0 |
| `status` | `scored`, `near_terminus` (fewer than 10 residues on one side) or `non_standard_residue` |

**Sort by `total` to rank sites.** Don't rely on `predicted` with the default
threshold; it marks nearly every lysine (see below). In the example file the
top hits are H4 positions 13 and 17 and H3 position 15, which are the known
acetylation sites H4K12, H4K16 and H3K14 (histone names skip the initial Met).

The example sequences are for trying the tool. Download them from UniProt
(P68431, P62805, P0CG48) before using them in research.

## Model

For every lysine with 10 residues on each side, the 21-residue window gets
an integrated score made of nine weighted parts:

| Column | Feature | Weight |
|---|---|---|
| `absolute_entropy` | Absolute entropy scale, summed over window | −0.026 |
| `average_volume` | Average volume of residues | −0.009 |
| `nonbonded_energy` | Average non-bonded energy per residue | −0.695 |
| `heat_capacity` | Heat capacity | 0.019 |
| `isoelectric_point` | Isoelectric point | 0.131 |
| `sequence_frequency` | Sequence frequency | 0.003 |
| `chou_fasman_coil` | Chou–Fasman coil-conformation parameter | 1.781 |
| `pssm` | Sum of position-specific residue frequencies from known sites | 17.66 |
| `dimer` | Sum of dipeptide scores over adjacent pairs | 5.553 |

The position profile is built from the training windows. By default these are
the 2,000 human positive sites in `data/human_positive_sites.txt`, which are
built into the binary.

### About the threshold

The original program called a site positive when `total > 212.368`. That value
was never calibrated for 21-residue windows. Known sites score about 284 ± 5,
and random K-centred windows score about 276 ± 6, so **nearly every lysine
passes**. The value is kept as the default only for compatibility. Every row
includes `total`, so you can rank sites or apply your own cutoff with
`--threshold`. A proper calibration needs a set of negative examples
(non-acetylated lysines), which this repository does not include.

## Usage

### Docker

```sh
docker build -t acetyl-scan .

# Score a FASTA file in the current directory
docker run --rm -v "$PWD":/data acetyl-scan scan proteins.fasta > scores.tsv

# Or stream through stdin
docker run --rm -i acetyl-scan scan - < proteins.fasta > scores.tsv

# Write to a file inside the mount (on Linux, run as your own user so the file is yours)
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/data acetyl-scan \
  scan proteins.fasta -o scores.tsv --threshold 285 --only-predicted
```

### Native

```sh
cargo build --release
./target/release/acetyl-scan scan proteins.fasta > scores.tsv
```

### Commands

| Command | Purpose |
|---|---|
| `scan <FASTA\|->` | One TSV row per lysine with all nine feature contributions, `total` and `predicted` |
| `profile` | The 20 × 21 position frequency table (columns −10 … +10) |
| `crossval [-k 5] [--seed 42]` | Scores each training site with a profile built *without* its fold, next to the in-sample score |

Common options:
- `--training FILE`: use other training windows (21 residues per line, K at the centre).
- `-o FILE`: write the TSV to a file instead of stdout.
- `--precision N`: number of decimal places in the score columns.

`scan` options:
- `--threshold X`: score cutoff for the `predicted` column.
- `--only-predicted`: write only sites above the threshold.
- `--include-skipped`: also write rows for lysines that could not be scored.

A summary goes to stderr, for example
`proteins=28 lysines=321 scored=267 predicted=44 skipped: near_terminus=54`.

### `scan` output columns

`protein_id` (first word of the header), `position` (1-based), `window`, the
nine feature columns above, `missing_dimers`, `total`, `predicted` (0/1),
`status` (`scored`, `near_terminus`, `non_standard_residue`).

## Fixes relative to the original Java

| Original behaviour | Now |
|---|---|
| Window bounds off by one: the first and last lysine with a full window were never scored | Every lysine with 10 residues on each side is scored |
| Infinite loop if the file didn't start with a `>` header | Clear error with line number |
| Rows without a newline ran together; positions were 0-based | One TSV row per lysine, 1-based positions |
| Lowercase or non-standard residues (X, U, B, Z) silently scored as 0 | Input is uppercased; windows with non-standard residues are skipped and counted |
| Dipeptides CW, MW, WM missing from the table, silently 0 | Still 0 (the original values are unknown), but counted per row in `missing_dimers` |
| Rounding to 4 decimals in the middle, applied inconsistently | Full precision throughout; rounding only when printing |
| Hard-coded `E:/…` paths | Command-line arguments, stdin/stdout |
| Training file not validated | Length, centre K and residues checked, with `file:line` errors |
| Training sites scored with their own profile (circular) | `crossval` gives held-out scores |
| ~2,000 lines of if/else chains | Lookup tables generated mechanically from the Java source |

## Verification

`tests/golden.rs` runs the port on `tests/fixtures/reference.fasta` and compares
it with `tests/fixtures/java_reference.tsv`. That reference was produced by the
original Java program, unchanged except for its file paths. All 262 sites the
Java program scored match, with a largest difference of 0.00019 (from removing
the intermediate rounding). The only extra sites are the edge lysines the Java
program skipped.

Run all checks without a local Rust toolchain:

```sh
docker build --target test .
```
