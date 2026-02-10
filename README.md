# qPCR Oligo Inclusivity Tool

A desktop application for evaluating how well a set of qPCR primers and probes match against a collection of reference sequences. Given a FASTA file of target sequences and user-defined forward primers, reverse primers, and probes, the tool aligns each oligo against every sequence and reports match statistics, mismatch patterns, and amplicon information.

Built with Rust using eframe/egui for the GUI and rust-bio for pairwise alignment.

## Building

Requires the Rust toolchain (install from https://rustup.rs).

```
cargo build --release
```

The compiled binary will be in `target/release/`. Release builds use LTO and single codegen unit, so compilation takes several minutes.

## Usage

1. **Load sequences** -- Click "Browse" to select a FASTA file containing the reference sequences to analyze.

2. **Enter oligos** -- Use the category tabs (Forward Primers, Reverse Primers, Probes) to enter oligo sequences in FASTA format. Forward and reverse primers are required. Probes are optional. Oligo sets can be saved/loaded as JSON files.

3. **Configure settings** -- Adjust match thresholds (minimum number of matched oligos per category), coverage requirements, maximum mismatches per oligo, and amplicon size constraints. The "Advanced Settings" dialog exposes alignment scoring parameters (match/mismatch/gap scores) and alignment mode (local or global).

4. **Run analysis** -- Click "Run Analysis". Progress is displayed in a progress bar. Results appear in a window with the full text report.

5. **Export** -- Results can be saved as plain text or exported to Excel (.xlsx).

## How the Analysis Works

For each reference sequence, the tool performs the following steps:

### Oligo Alignment

Each oligo (primer or probe) is aligned against the reference sequence in both the sense and antisense orientations using pairwise alignment from rust-bio. The orientation with the higher alignment score is selected. An oligo is considered matched if:

- The alignment coverage (fraction of the oligo length covered by match/substitution operations) meets the minimum coverage threshold (default 0.8).
- The number of mismatches does not exceed the per-oligo maximum (default 7).

IUPAC ambiguity codes are supported. Ambiguous bases that are compatible with the aligned target base are counted as matches. N is always counted as a mismatch.

### Amplicon Detection

After aligning forward and reverse primers independently, the tool looks for convergent primer pairs -- a forward primer whose alignment start position is upstream of a reverse primer's alignment end position. If amplicon size constraints are enabled, only pairs producing an amplicon within the specified size range are accepted. When multiple valid pairs exist, the largest amplicon is selected.

If no valid forward+reverse pair is found for a sequence, all oligo matches for that sequence are discarded. When a valid amplicon is found, any oligo match (including other forward/reverse primers and probes) that falls outside the amplicon boundaries is also discarded. Probes are only aligned if at least one forward and one reverse primer initially match.

### Pattern Aggregation

Each sequence that meets the per-category minimum match thresholds produces a signature string. The signature encodes, for each oligo, either the mismatch pattern (positions where the reference differs from the oligo) or `NO_MATCH`. Signatures are formatted as `fwd1 | fwd2 || probe1 || rev1 | rev2`, with `|` separating oligos within a category and `||` separating categories.

Sequences sharing identical signatures are grouped into patterns. The output reports each pattern with its count, percentage, mismatch totals, and example sequence IDs.

### Output Statistics

The results include:

- Per-oligo match counts and sense/antisense breakdown.
- Per-category mismatch distributions (0 mismatches, 1 mismatch, >1 mismatches, no match), based on the best-matching oligo in each category.
- Overall mismatch quality (all categories perfect, all within 1 mismatch, or 2+ mismatches in any category).
- Amplicon statistics (valid amplicon count, mean/min/max amplicon size).
- Pattern table sorted by frequency.

All sequence processing is parallelized using rayon.
