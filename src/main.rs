//! qPCR Oligo Inclusivity Tool v2.0
//! A fast Rust implementation with eframe GUI
//!
//! This tool evaluates primer/probe inclusivity in large sets of reference sequences
//! using pairwise alignment from rust-bio. Users designate forward primers, reverse
//! primers, and probes separately. The tool searches for amplificates between forward
//! and reverse primers and matches probes within the amplicon region.

use bio::alignment::pairwise::{Aligner, Scoring};
use bio::alignment::AlignmentOperation;
use bio::io::fasta;
use eframe::egui;
use rayon::prelude::*;
use rfd::FileDialog;
use rust_xlsxwriter::{Format, Workbook};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

// ============================================================================
// Data Structures
// ============================================================================

/// A FASTA record with id and sequence
#[derive(Clone, Debug)]
struct FastaRecord {
    id: String,
    seq: String,
}

/// Oligo/Primer sequence with id
#[derive(Clone, Debug)]
struct Oligo {
    id: String,
    seq: String,
}

/// Serializable oligo for JSON save/load
#[derive(Clone, Debug, Serialize, Deserialize)]
struct OligoSerde {
    id: String,
    seq: String,
}

/// JSON-serializable set of all oligo categories
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
struct OligoSet {
    forward_primers: Vec<OligoSerde>,
    reverse_primers: Vec<OligoSerde>,
    probes: Vec<OligoSerde>,
}

/// Which oligo category is currently active in the UI
#[derive(Clone, Debug, Copy, PartialEq)]
enum OligoCategory {
    ForwardPrimer,
    ReversePrimer,
    Probe,
}

/// Result of aligning an oligo to a sequence
#[derive(Clone, Debug)]
struct OligoResult {
    matched: bool,
    orientation: Option<Orientation>,
    signature: String,
    mismatches: usize,
    score: i32,
    coverage: f64,
    start_pos: Option<usize>,
    end_pos: Option<usize>,
}

impl OligoResult {
    fn no_match() -> Self {
        Self {
            matched: false,
            orientation: None,
            signature: String::new(),
            mismatches: 0,
            score: 0,
            coverage: 0.0,
            start_pos: None,
            end_pos: None,
        }
    }
}

/// Orientation of the oligo match
#[derive(Clone, Debug, Copy, PartialEq)]
enum Orientation {
    Sense,
    Antisense,
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Orientation::Sense => write!(f, "fwd"),
            Orientation::Antisense => write!(f, "rev"),
        }
    }
}

/// Information about the selected amplicon
#[derive(Clone, Debug, Default)]
struct AmpliconInfo {
    found: bool,
    forward_oligo_id: Option<String>,
    reverse_oligo_id: Option<String>,
    start: usize,
    end: usize,
    size: usize,
}

/// Aggregated pattern data
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
struct PatternData {
    count: usize,
    total_mismatches: usize,
    matched_fwd: usize,
    matched_rev: usize,
    matched_probe: usize,
    examples: Vec<String>,
    amplicon_lengths: Vec<usize>,
}

/// Oligo match statistics
#[derive(Clone, Debug, Default)]
struct OligoStats {
    total_matches: usize,
    sense_matches: usize,
    antisense_matches: usize,
}

/// Alignment settings
#[derive(Clone, Debug, Serialize, Deserialize)]
struct AlignmentSettings {
    mode: AlignmentMode,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    min_fwd_matched: usize,
    min_rev_matched: usize,
    min_probe_matched: usize,
    min_coverage: f64,
    max_mismatches_per_oligo: usize,
    ambiguity_display: AmbiguityDisplayMode,
    min_amplicon_size: Option<usize>,
    max_amplicon_size: Option<usize>,
}

impl Default for AlignmentSettings {
    fn default() -> Self {
        Self {
            mode: AlignmentMode::Local,
            match_score: 2,
            mismatch_score: -1,
            gap_open_score: -2,
            gap_extend_score: -1,
            min_fwd_matched: 1,
            min_rev_matched: 1,
            min_probe_matched: 0,
            min_coverage: 0.8,
            max_mismatches_per_oligo: 7,
            ambiguity_display: AmbiguityDisplayMode::ShowDots,
            min_amplicon_size: None,
            max_amplicon_size: None,
        }
    }
}

#[derive(Clone, Debug, Copy, PartialEq, Serialize, Deserialize)]
enum AlignmentMode {
    Local,
    Global,
}

#[derive(Clone, Debug, Copy, PartialEq, Serialize, Deserialize)]
enum AmbiguityDisplayMode {
    ShowBases,
    ShowDots,
}

/// Per-sequence analysis result
#[derive(Clone, Debug)]
struct SequenceResult {
    fwd_results: HashMap<String, OligoResult>,
    rev_results: HashMap<String, OligoResult>,
    probe_results: HashMap<String, OligoResult>,
    fwd_matched: usize,
    rev_matched: usize,
    probe_matched: usize,
    total_mismatches: usize,
    amplicon_info: Option<AmpliconInfo>,
    /// Best (minimum) mismatch count among matched fwd primers, None if no match
    best_fwd_mm: Option<usize>,
    /// Best (minimum) mismatch count among matched rev primers, None if no match
    best_rev_mm: Option<usize>,
    /// Best (minimum) mismatch count among matched probes, None if no match
    best_probe_mm: Option<usize>,
}

/// Per-category mismatch distribution counts
#[derive(Clone, Debug, Default)]
struct MismatchDistribution {
    zero_mm: usize,   // best oligo in category has 0 mismatches
    one_mm: usize,    // best oligo in category has exactly 1 mismatch
    more_mm: usize,   // best oligo in category has >1 mismatches
    no_match: usize,  // no oligo matched in this category
}

/// Analysis results
#[derive(Clone, Debug, Default)]
struct AnalysisResults {
    alignment_dict: HashMap<String, PatternData>,
    oligo_stats: HashMap<String, OligoStats>,
    total_sequences: usize,
    sequences_with_min_matches: usize,
    sequences_with_valid_amplicon: usize,
    sequences_failed_amplicon: usize,
    fwd_mm_dist: MismatchDistribution,
    rev_mm_dist: MismatchDistribution,
    probe_mm_dist: MismatchDistribution,
    /// Overall pattern: all active categories have best oligo with 0 mismatches
    overall_all_perfect: usize,
    /// Overall pattern: worst category has exactly 1 mismatch (all ≤1mm but not all 0)
    overall_max_one_mm: usize,
    /// Overall pattern: at least one category has ≥2 mismatches
    overall_two_plus_mm: usize,
    /// Overall pattern: at least one category has no match at all
    overall_no_match: usize,
    output_text: String,
}

/// Progress tracking for analysis
#[derive(Clone)]
struct ProgressTracker {
    current: Arc<AtomicUsize>,
    total: Arc<AtomicUsize>,
    status: Arc<Mutex<String>>,
    running: Arc<AtomicBool>,
}

impl Default for ProgressTracker {
    fn default() -> Self {
        Self {
            current: Arc::new(AtomicUsize::new(0)),
            total: Arc::new(AtomicUsize::new(0)),
            status: Arc::new(Mutex::new("Ready".to_string())),
            running: Arc::new(AtomicBool::new(false)),
        }
    }
}

// ============================================================================
// Bioinformatics Functions
// ============================================================================

/// Check if two nucleotides match considering IUPAC ambiguity codes
/// Returns true if they match, false otherwise
/// N always returns false (mismatch) regardless of the other base
fn iupac_match(a: u8, b: u8) -> bool {
    let a_upper = a.to_ascii_uppercase();
    let b_upper = b.to_ascii_uppercase();

    // N always mismatches
    if a_upper == b'N' || b_upper == b'N' {
        return false;
    }

    // Exact match
    if a_upper == b_upper {
        return true;
    }

    // Define IUPAC ambiguity codes and their possible bases
    let get_bases = |code: u8| -> &'static [u8] {
        match code.to_ascii_uppercase() {
            b'A' => &[b'A'],
            b'T' => &[b'T'],
            b'G' => &[b'G'],
            b'C' => &[b'C'],
            b'U' => &[b'U'],
            b'R' => &[b'A', b'G'],
            b'Y' => &[b'C', b'T'],
            b'S' => &[b'G', b'C'],
            b'W' => &[b'A', b'T'],
            b'K' => &[b'G', b'T'],
            b'M' => &[b'A', b'C'],
            b'B' => &[b'C', b'G', b'T'],
            b'D' => &[b'A', b'G', b'T'],
            b'H' => &[b'A', b'C', b'T'],
            b'V' => &[b'A', b'C', b'G'],
            _ => &[],
        }
    };

    let a_bases = get_bases(a_upper);
    let b_bases = get_bases(b_upper);

    for &base_a in a_bases {
        for &base_b in b_bases {
            let normalized_a = if base_a == b'U' { b'T' } else { base_a };
            let normalized_b = if base_b == b'U' { b'T' } else { base_b };
            if normalized_a == normalized_b {
                return true;
            }
        }
    }

    false
}

/// Get reverse complement of a DNA sequence
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'U' => 'A',
            'R' => 'Y',
            'Y' => 'R',
            'S' => 'S',
            'W' => 'W',
            'K' => 'M',
            'M' => 'K',
            'B' => 'V',
            'D' => 'H',
            'H' => 'D',
            'V' => 'B',
            'N' => 'N',
            _ => 'N',
        })
        .collect()
}

/// Parse a FASTA file and return records
fn parse_fasta(path: &PathBuf) -> Result<Vec<FastaRecord>, String> {
    let file = File::open(path).map_err(|e| format!("Failed to open file: {}", e))?;
    let reader = fasta::Reader::new(BufReader::new(file));

    let mut records = Vec::new();
    for result in reader.records() {
        match result {
            Ok(record) => {
                let id = record.id().to_string();
                let seq = std::str::from_utf8(record.seq())
                    .map_err(|e| format!("Invalid UTF-8 in sequence: {}", e))?
                    .to_uppercase()
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .collect();
                records.push(FastaRecord { id, seq });
            }
            Err(e) => return Err(format!("Error reading FASTA record: {}", e)),
        }
    }
    Ok(records)
}

/// Parse FASTA format from a string and return oligos
fn parse_fasta_string(text: &str) -> Result<Vec<Oligo>, String> {
    let cursor = std::io::Cursor::new(text);
    let reader = fasta::Reader::new(cursor);

    let mut oligos = Vec::new();
    for result in reader.records() {
        match result {
            Ok(record) => {
                let id = record.id().to_string();
                let seq = std::str::from_utf8(record.seq())
                    .map_err(|e| format!("Invalid UTF-8 in sequence: {}", e))?
                    .to_uppercase()
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .collect();
                oligos.push(Oligo { id, seq });
            }
            Err(e) => return Err(format!("Error parsing FASTA: {}", e)),
        }
    }
    Ok(oligos)
}

/// Convert oligos to FASTA format string
fn oligos_to_fasta_string(oligos: &[Oligo]) -> String {
    oligos
        .iter()
        .map(|o| format!(">{}\n{}", o.id, o.seq))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Perform alignment and return best alignment with orientation
fn get_best_alignment(
    target_seq: &str,
    oligo_seq: &str,
    settings: &AlignmentSettings,
) -> (
    Option<bio::alignment::Alignment>,
    Option<Orientation>,
    String,
) {
    let scoring = Scoring::new(
        settings.gap_open_score,
        settings.gap_extend_score,
        |a: u8, b: u8| {
            if iupac_match(a, b) {
                settings.match_score
            } else {
                settings.mismatch_score
            }
        },
    );

    let target_bytes = target_seq.as_bytes();
    let oligo_bytes = oligo_seq.as_bytes();
    let oligo_rc = reverse_complement(oligo_seq);
    let oligo_rc_bytes = oligo_rc.as_bytes();

    let mut aligner = Aligner::with_capacity_and_scoring(
        target_bytes.len(),
        oligo_bytes.len(),
        scoring.clone(),
    );

    let alignment_sense = match settings.mode {
        AlignmentMode::Local => aligner.local(target_bytes, oligo_bytes),
        AlignmentMode::Global => aligner.global(target_bytes, oligo_bytes),
    };

    let alignment_antisense = match settings.mode {
        AlignmentMode::Local => aligner.local(target_bytes, oligo_rc_bytes),
        AlignmentMode::Global => aligner.global(target_bytes, oligo_rc_bytes),
    };

    if alignment_sense.score >= alignment_antisense.score {
        (
            Some(alignment_sense),
            Some(Orientation::Sense),
            oligo_seq.to_string(),
        )
    } else {
        (
            Some(alignment_antisense),
            Some(Orientation::Antisense),
            oligo_rc,
        )
    }
}

/// Check if alignment is valid based on coverage criteria
fn is_valid_alignment(
    alignment: &bio::alignment::Alignment,
    oligo_len: usize,
    min_coverage: f64,
) -> bool {
    let aligned_length = alignment
        .operations
        .iter()
        .filter(|op| matches!(op, AlignmentOperation::Match | AlignmentOperation::Subst))
        .count();

    let coverage = aligned_length as f64 / oligo_len as f64;
    coverage >= min_coverage
}

/// Calculate alignment coverage
fn get_alignment_coverage(alignment: &bio::alignment::Alignment, oligo_len: usize) -> f64 {
    let aligned_length = alignment
        .operations
        .iter()
        .filter(|op| matches!(op, AlignmentOperation::Match | AlignmentOperation::Subst))
        .count();

    aligned_length as f64 / oligo_len as f64
}

/// Check if two nucleotides are an exact match (not just IUPAC ambiguity match)
fn is_exact_match(a: u8, b: u8) -> bool {
    let a_upper = a.to_ascii_uppercase();
    let b_upper = b.to_ascii_uppercase();

    let normalized_a = if a_upper == b'U' { b'T' } else { a_upper };
    let normalized_b = if b_upper == b'U' { b'T' } else { b_upper };

    normalized_a == normalized_b
}

/// Generate signature from alignment
/// For reverse orientation hits, the pattern is reversed to match the original oligo orientation
fn generate_signature(
    alignment: &bio::alignment::Alignment,
    target_seq: &str,
    oligo_seq: &str,
    orientation: Orientation,
    ambiguity_display: AmbiguityDisplayMode,
) -> (String, usize) {
    let oligo_len = oligo_seq.len();
    let mut sig: Vec<char> = vec!['-'; oligo_len];
    let mut mismatches = 0;

    let target_bytes = target_seq.as_bytes();
    let oligo_bytes = oligo_seq.as_bytes();

    let mut t_pos = alignment.xstart;
    let mut q_pos = alignment.ystart;

    for op in &alignment.operations {
        match op {
            AlignmentOperation::Match => {
                if q_pos < oligo_len && t_pos < target_bytes.len() {
                    let target_base = target_bytes[t_pos];
                    let oligo_base = oligo_bytes[q_pos];

                    if is_exact_match(target_base, oligo_base) {
                        sig[q_pos] = '.';
                    } else {
                        match ambiguity_display {
                            AmbiguityDisplayMode::ShowDots => {
                                sig[q_pos] = '.';
                            }
                            AmbiguityDisplayMode::ShowBases => {
                                sig[q_pos] = target_bytes[t_pos] as char;
                                mismatches += 1;
                            }
                        }
                    }
                }
                t_pos += 1;
                q_pos += 1;
            }
            AlignmentOperation::Subst => {
                if q_pos < oligo_len && t_pos < target_bytes.len() {
                    let target_base = target_bytes[t_pos];
                    let oligo_base = oligo_bytes[q_pos];

                    if iupac_match(target_base, oligo_base) {
                        match ambiguity_display {
                            AmbiguityDisplayMode::ShowDots => {
                                sig[q_pos] = '.';
                            }
                            AmbiguityDisplayMode::ShowBases => {
                                sig[q_pos] = target_base as char;
                                mismatches += 1;
                            }
                        }
                    } else {
                        sig[q_pos] = target_base as char;
                        mismatches += 1;
                    }
                }
                t_pos += 1;
                q_pos += 1;
            }
            AlignmentOperation::Del => {
                t_pos += 1;
            }
            AlignmentOperation::Ins => {
                if q_pos < oligo_len {
                    sig[q_pos] = '-';
                    mismatches += 1;
                }
                q_pos += 1;
            }
            AlignmentOperation::Xclip(_) => {}
            AlignmentOperation::Yclip(_) => {}
        }
    }

    mismatches += sig.iter().filter(|&&c| c == '-').count();

    let mut signature: String = sig.into_iter().collect();

    if orientation == Orientation::Antisense {
        signature = signature
            .chars()
            .map(|c| match c.to_ascii_uppercase() {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                'U' => 'A',
                'R' => 'Y',
                'Y' => 'R',
                'S' => 'S',
                'W' => 'W',
                'K' => 'M',
                'M' => 'K',
                'B' => 'V',
                'D' => 'H',
                'H' => 'D',
                'V' => 'B',
                _ => c,
            })
            .collect();
        signature = signature.chars().rev().collect();
    }

    (signature, mismatches)
}

/// Align a set of oligos against a target sequence, returning results map
/// Oligos exceeding max_mismatches_per_oligo are discarded as NO_MATCH
fn align_oligos(
    sequence: &FastaRecord,
    oligos: &[Oligo],
    settings: &AlignmentSettings,
) -> HashMap<String, OligoResult> {
    let mut results = HashMap::new();

    for oligo in oligos {
        let (alignment_opt, orientation_opt, actual_oligo) =
            get_best_alignment(&sequence.seq, &oligo.seq, settings);

        let result =
            if let (Some(alignment), Some(orientation)) = (alignment_opt, orientation_opt) {
                if is_valid_alignment(&alignment, actual_oligo.len(), settings.min_coverage) {
                    let (signature, mismatches) = generate_signature(
                        &alignment,
                        &sequence.seq,
                        &actual_oligo,
                        orientation,
                        settings.ambiguity_display,
                    );

                    // Discard if mismatches exceed threshold
                    if mismatches > settings.max_mismatches_per_oligo {
                        OligoResult::no_match()
                    } else {
                        let coverage = get_alignment_coverage(&alignment, actual_oligo.len());

                        OligoResult {
                            matched: true,
                            orientation: Some(orientation),
                            signature,
                            mismatches,
                            score: alignment.score,
                            coverage,
                            start_pos: Some(alignment.xstart),
                            end_pos: Some(alignment.xend),
                        }
                    }
                } else {
                    OligoResult::no_match()
                }
            } else {
                OligoResult::no_match()
            };

        results.insert(oligo.id.clone(), result);
    }

    results
}

/// Find the best amplicon from forward and reverse primer results.
/// Selects the largest convergent fwd+rev pair within size constraints.
/// Filters ALL oligo categories: only matches within the amplicon bounds are kept.
fn find_best_amplicon(
    fwd_results: &mut HashMap<String, OligoResult>,
    rev_results: &mut HashMap<String, OligoResult>,
    probe_results: &mut HashMap<String, OligoResult>,
    min_size: Option<usize>,
    max_size: Option<usize>,
) -> AmpliconInfo {
    // Collect matched forward primers with positions
    let forward_hits: Vec<(String, usize, usize)> = fwd_results
        .iter()
        .filter(|(_, r)| r.matched)
        .filter_map(|(id, r)| match (r.start_pos, r.end_pos) {
            (Some(s), Some(e)) => Some((id.clone(), s, e)),
            _ => None,
        })
        .collect();

    // Collect matched reverse primers with positions
    let reverse_hits: Vec<(String, usize, usize)> = rev_results
        .iter()
        .filter(|(_, r)| r.matched)
        .filter_map(|(id, r)| match (r.start_pos, r.end_pos) {
            (Some(s), Some(e)) => Some((id.clone(), s, e)),
            _ => None,
        })
        .collect();

    // Find all convergent pairs within size constraints
    // Convergent: forward start < reverse end (primers face each other)
    let mut valid_amplicons: Vec<(String, String, usize, usize, usize)> = Vec::new();

    for (fwd_id, fwd_start, _fwd_end) in &forward_hits {
        for (rev_id, _rev_start, rev_end) in &reverse_hits {
            if fwd_start < rev_end {
                let amplicon_size = rev_end - fwd_start + 1;
                let size_ok = match (min_size, max_size) {
                    (Some(min), Some(max)) => amplicon_size >= min && amplicon_size <= max,
                    (Some(min), None) => amplicon_size >= min,
                    (None, Some(max)) => amplicon_size <= max,
                    (None, None) => true,
                };
                if size_ok {
                    valid_amplicons.push((
                        fwd_id.clone(),
                        rev_id.clone(),
                        *fwd_start,
                        *rev_end,
                        amplicon_size,
                    ));
                }
            }
        }
    }

    // Helper: mark all results in a map as no_match
    let mark_all_no_match = |results: &mut HashMap<String, OligoResult>| {
        for result in results.values_mut() {
            *result = OligoResult::no_match();
        }
    };

    if valid_amplicons.is_empty() {
        // No valid amplicon: ALL categories become no_match
        mark_all_no_match(fwd_results);
        mark_all_no_match(rev_results);
        mark_all_no_match(probe_results);
        return AmpliconInfo::default();
    }

    // Select the LARGEST valid amplicon
    let best = valid_amplicons
        .iter()
        .max_by_key(|a| a.4)
        .unwrap()
        .clone();

    let (best_fwd_id, best_rev_id, amp_start, amp_end, amp_size) = best;

    // Filter all categories: only matches entirely within [amp_start, amp_end] are kept
    let filter_by_bounds = |results: &mut HashMap<String, OligoResult>| {
        for result in results.values_mut() {
            if !result.matched {
                continue;
            }
            if let (Some(start), Some(end)) = (result.start_pos, result.end_pos) {
                if start < amp_start || end > amp_end {
                    *result = OligoResult::no_match();
                }
            } else {
                *result = OligoResult::no_match();
            }
        }
    };

    filter_by_bounds(fwd_results);
    filter_by_bounds(rev_results);
    filter_by_bounds(probe_results);

    AmpliconInfo {
        found: true,
        forward_oligo_id: Some(best_fwd_id),
        reverse_oligo_id: Some(best_rev_id),
        start: amp_start,
        end: amp_end,
        size: amp_size,
    }
}

/// Analyze a sequence against categorized oligos (forward primers, reverse primers, probes)
fn analyze_sequence(
    sequence: &FastaRecord,
    fwd_primers: &[Oligo],
    rev_primers: &[Oligo],
    probes: &[Oligo],
    settings: &AlignmentSettings,
) -> SequenceResult {
    // Align forward and reverse primers
    let mut fwd_results = align_oligos(sequence, fwd_primers, settings);
    let mut rev_results = align_oligos(sequence, rev_primers, settings);

    // Probes are only aligned if we have at least one fwd and rev match
    let mut probe_results: HashMap<String, OligoResult> = HashMap::new();
    let amplicon_info;

    let has_fwd_match = fwd_results.values().any(|r| r.matched);
    let has_rev_match = rev_results.values().any(|r| r.matched);

    if has_fwd_match && has_rev_match {
        // Align probes, then find amplicon and filter ALL categories by its bounds
        probe_results = align_oligos(sequence, probes, settings);

        amplicon_info = Some(find_best_amplicon(
            &mut fwd_results,
            &mut rev_results,
            &mut probe_results,
            settings.min_amplicon_size,
            settings.max_amplicon_size,
        ));
    } else {
        // No valid amplicon possible - ALL categories become NO_MATCH
        for result in fwd_results.values_mut() {
            *result = OligoResult::no_match();
        }
        for result in rev_results.values_mut() {
            *result = OligoResult::no_match();
        }
        for probe in probes {
            probe_results.insert(probe.id.clone(), OligoResult::no_match());
        }
        amplicon_info = Some(AmpliconInfo::default());
    }

    // Count matches per category
    let fwd_matched = fwd_results.values().filter(|r| r.matched).count();
    let rev_matched = rev_results.values().filter(|r| r.matched).count();
    let probe_matched = probe_results.values().filter(|r| r.matched).count();

    // Best (minimum) mismatch count per category
    let best_fwd_mm = fwd_results.values().filter(|r| r.matched).map(|r| r.mismatches).min();
    let best_rev_mm = rev_results.values().filter(|r| r.matched).map(|r| r.mismatches).min();
    let best_probe_mm = probe_results.values().filter(|r| r.matched).map(|r| r.mismatches).min();

    // Total mismatches: sum of per-category minimums (best oligo per category)
    let total_mismatches: usize = best_fwd_mm.unwrap_or(0)
        + best_rev_mm.unwrap_or(0)
        + best_probe_mm.unwrap_or(0);

    SequenceResult {
        fwd_results,
        rev_results,
        probe_results,
        fwd_matched,
        rev_matched,
        probe_matched,
        total_mismatches,
        amplicon_info,
        best_fwd_mm,
        best_rev_mm,
        best_probe_mm,
    }
}

/// Build a signature part for a set of oligo results in the given oligo order
fn build_signature_parts(
    oligos: &[Oligo],
    results: &HashMap<String, OligoResult>,
) -> Vec<String> {
    oligos
        .iter()
        .map(|oligo| {
            if let Some(result) = results.get(&oligo.id) {
                if result.matched {
                    let orientation_symbol = match result.orientation {
                        Some(Orientation::Sense) => "(fwd)",
                        Some(Orientation::Antisense) => "(rev)",
                        None => "",
                    };
                    format!("{}{}", result.signature, orientation_symbol)
                } else {
                    "NO_MATCH".to_string()
                }
            } else {
                "NO_MATCH".to_string()
            }
        })
        .collect()
}

/// Run the full analysis with parallel processing
fn run_analysis(
    sequences: &[FastaRecord],
    fwd_primers: &[Oligo],
    rev_primers: &[Oligo],
    probes: &[Oligo],
    settings: &AlignmentSettings,
    progress: &ProgressTracker,
) -> AnalysisResults {
    let total_sequences = sequences.len();
    progress.total.store(total_sequences, Ordering::SeqCst);
    progress.current.store(0, Ordering::SeqCst);

    // Initialize oligo stats for all categories
    let mut initial_stats: HashMap<String, OligoStats> = HashMap::new();
    for oligo in fwd_primers.iter().chain(rev_primers.iter()).chain(probes.iter()) {
        initial_stats.insert(oligo.id.clone(), OligoStats::default());
    }
    let oligo_stats = Arc::new(Mutex::new(initial_stats));

    let alignment_dict: Arc<Mutex<HashMap<String, PatternData>>> =
        Arc::new(Mutex::new(HashMap::new()));
    let sequences_with_min_matches = Arc::new(AtomicUsize::new(0));
    let sequences_with_valid_amplicon = Arc::new(AtomicUsize::new(0));
    let sequences_failed_amplicon = Arc::new(AtomicUsize::new(0));

    // Mismatch distribution counters per category
    let fwd_mm_zero = Arc::new(AtomicUsize::new(0));
    let fwd_mm_one = Arc::new(AtomicUsize::new(0));
    let fwd_mm_more = Arc::new(AtomicUsize::new(0));
    let fwd_mm_none = Arc::new(AtomicUsize::new(0));
    let rev_mm_zero = Arc::new(AtomicUsize::new(0));
    let rev_mm_one = Arc::new(AtomicUsize::new(0));
    let rev_mm_more = Arc::new(AtomicUsize::new(0));
    let rev_mm_none = Arc::new(AtomicUsize::new(0));
    let probe_mm_zero = Arc::new(AtomicUsize::new(0));
    let probe_mm_one = Arc::new(AtomicUsize::new(0));
    let probe_mm_more = Arc::new(AtomicUsize::new(0));
    let probe_mm_none = Arc::new(AtomicUsize::new(0));

    // Overall pattern counters
    let overall_all_perfect = Arc::new(AtomicUsize::new(0));
    let overall_max_one_mm = Arc::new(AtomicUsize::new(0));
    let overall_two_plus_mm = Arc::new(AtomicUsize::new(0));
    let overall_no_match = Arc::new(AtomicUsize::new(0));

    let processed = Arc::new(AtomicUsize::new(0));

    sequences.par_iter().for_each(|record| {
        if !progress.running.load(Ordering::SeqCst) {
            return;
        }

        let seq_result = analyze_sequence(record, fwd_primers, rev_primers, probes, settings);

        // Track amplicon statistics
        if let Some(ref info) = seq_result.amplicon_info {
            if info.found {
                sequences_with_valid_amplicon.fetch_add(1, Ordering::SeqCst);
            } else {
                sequences_failed_amplicon.fetch_add(1, Ordering::SeqCst);
            }
        }

        // Track per-category mismatch distributions (best oligo per category)
        match seq_result.best_fwd_mm {
            Some(0) => { fwd_mm_zero.fetch_add(1, Ordering::SeqCst); }
            Some(1) => { fwd_mm_one.fetch_add(1, Ordering::SeqCst); }
            Some(_) => { fwd_mm_more.fetch_add(1, Ordering::SeqCst); }
            None => { fwd_mm_none.fetch_add(1, Ordering::SeqCst); }
        }
        match seq_result.best_rev_mm {
            Some(0) => { rev_mm_zero.fetch_add(1, Ordering::SeqCst); }
            Some(1) => { rev_mm_one.fetch_add(1, Ordering::SeqCst); }
            Some(_) => { rev_mm_more.fetch_add(1, Ordering::SeqCst); }
            None => { rev_mm_none.fetch_add(1, Ordering::SeqCst); }
        }
        if !probes.is_empty() {
            match seq_result.best_probe_mm {
                Some(0) => { probe_mm_zero.fetch_add(1, Ordering::SeqCst); }
                Some(1) => { probe_mm_one.fetch_add(1, Ordering::SeqCst); }
                Some(_) => { probe_mm_more.fetch_add(1, Ordering::SeqCst); }
                None => { probe_mm_none.fetch_add(1, Ordering::SeqCst); }
            }
        }

        // Overall pattern: worst best-per-category across all active categories
        {
            let mut category_bests: Vec<Option<usize>> = vec![
                seq_result.best_fwd_mm,
                seq_result.best_rev_mm,
            ];
            if !probes.is_empty() {
                category_bests.push(seq_result.best_probe_mm);
            }

            if category_bests.iter().any(|b| b.is_none()) {
                // At least one category has no match at all
                overall_no_match.fetch_add(1, Ordering::SeqCst);
            } else {
                let worst_best = category_bests.iter().filter_map(|b| *b).max().unwrap_or(0);
                match worst_best {
                    0 => { overall_all_perfect.fetch_add(1, Ordering::SeqCst); }
                    1 => { overall_max_one_mm.fetch_add(1, Ordering::SeqCst); }
                    _ => { overall_two_plus_mm.fetch_add(1, Ordering::SeqCst); }
                }
            }
        }

        // Update oligo stats
        {
            let mut stats = oligo_stats.lock().unwrap();
            let all_results = seq_result
                .fwd_results
                .iter()
                .chain(seq_result.rev_results.iter())
                .chain(seq_result.probe_results.iter());
            for (oligo_id, result) in all_results {
                if result.matched {
                    if let Some(s) = stats.get_mut(oligo_id) {
                        s.total_matches += 1;
                        match result.orientation {
                            Some(Orientation::Sense) => s.sense_matches += 1,
                            Some(Orientation::Antisense) => s.antisense_matches += 1,
                            None => {}
                        }
                    }
                }
            }
        }

        // Check per-category minimum match thresholds
        let meets_fwd = seq_result.fwd_matched >= settings.min_fwd_matched;
        let meets_rev = seq_result.rev_matched >= settings.min_rev_matched;
        let meets_probe = seq_result.probe_matched >= settings.min_probe_matched;

        if meets_fwd && meets_rev && meets_probe {
            sequences_with_min_matches.fetch_add(1, Ordering::SeqCst);

            // Build combined signature: FWD || PROBES || REV
            let fwd_parts = build_signature_parts(fwd_primers, &seq_result.fwd_results);
            let probe_parts = build_signature_parts(probes, &seq_result.probe_results);
            let rev_parts = build_signature_parts(rev_primers, &seq_result.rev_results);

            let mut all_parts = Vec::new();
            all_parts.push(fwd_parts.join(" | "));
            if !probes.is_empty() {
                all_parts.push(probe_parts.join(" | "));
            }
            all_parts.push(rev_parts.join(" | "));

            let combined_signature = all_parts.join(" || ");

            // Update alignment dict
            {
                let mut dict = alignment_dict.lock().unwrap();
                let entry = dict
                    .entry(combined_signature)
                    .or_insert_with(|| PatternData {
                        count: 0,
                        total_mismatches: seq_result.total_mismatches,
                        matched_fwd: seq_result.fwd_matched,
                        matched_rev: seq_result.rev_matched,
                        matched_probe: seq_result.probe_matched,
                        examples: Vec::new(),
                        amplicon_lengths: Vec::new(),
                    });
                entry.count += 1;
                if entry.examples.len() < 10 {
                    entry.examples.push(record.id.clone());
                }
                // Always record amplicon length when a valid fwd+rev pair was found
                if let Some(ref info) = seq_result.amplicon_info {
                    if info.found && info.size > 0 {
                        entry.amplicon_lengths.push(info.size);
                    }
                }
            }
        }

        // Update progress
        let current = processed.fetch_add(1, Ordering::SeqCst) + 1;
        progress.current.store(current, Ordering::SeqCst);

        if current % 100 == 0 || current == total_sequences {
            if let Ok(mut status) = progress.status.lock() {
                *status = format!("Processing sequence {}/{}", current, total_sequences);
            }
        }
    });

    let final_oligo_stats = oligo_stats.lock().unwrap().clone();
    let final_alignment_dict = alignment_dict.lock().unwrap().clone();
    let final_sequences_with_min_matches = sequences_with_min_matches.load(Ordering::SeqCst);
    let final_sequences_with_valid_amplicon =
        sequences_with_valid_amplicon.load(Ordering::SeqCst);
    let final_sequences_failed_amplicon = sequences_failed_amplicon.load(Ordering::SeqCst);

    let final_fwd_mm_dist = MismatchDistribution {
        zero_mm: fwd_mm_zero.load(Ordering::SeqCst),
        one_mm: fwd_mm_one.load(Ordering::SeqCst),
        more_mm: fwd_mm_more.load(Ordering::SeqCst),
        no_match: fwd_mm_none.load(Ordering::SeqCst),
    };
    let final_rev_mm_dist = MismatchDistribution {
        zero_mm: rev_mm_zero.load(Ordering::SeqCst),
        one_mm: rev_mm_one.load(Ordering::SeqCst),
        more_mm: rev_mm_more.load(Ordering::SeqCst),
        no_match: rev_mm_none.load(Ordering::SeqCst),
    };
    let final_probe_mm_dist = MismatchDistribution {
        zero_mm: probe_mm_zero.load(Ordering::SeqCst),
        one_mm: probe_mm_one.load(Ordering::SeqCst),
        more_mm: probe_mm_more.load(Ordering::SeqCst),
        no_match: probe_mm_none.load(Ordering::SeqCst),
    };

    let final_overall_all_perfect = overall_all_perfect.load(Ordering::SeqCst);
    let final_overall_max_one_mm = overall_max_one_mm.load(Ordering::SeqCst);
    let final_overall_two_plus_mm = overall_two_plus_mm.load(Ordering::SeqCst);
    let final_overall_no_match = overall_no_match.load(Ordering::SeqCst);

    let output_text = generate_output_text(
        fwd_primers,
        rev_primers,
        probes,
        &final_oligo_stats,
        &final_alignment_dict,
        total_sequences,
        final_sequences_with_min_matches,
        final_sequences_with_valid_amplicon,
        final_sequences_failed_amplicon,
        &final_fwd_mm_dist,
        &final_rev_mm_dist,
        &final_probe_mm_dist,
        final_overall_all_perfect,
        final_overall_max_one_mm,
        final_overall_two_plus_mm,
        final_overall_no_match,
        settings,
    );

    AnalysisResults {
        alignment_dict: final_alignment_dict,
        oligo_stats: final_oligo_stats,
        total_sequences,
        sequences_with_min_matches: final_sequences_with_min_matches,
        sequences_with_valid_amplicon: final_sequences_with_valid_amplicon,
        sequences_failed_amplicon: final_sequences_failed_amplicon,
        fwd_mm_dist: final_fwd_mm_dist,
        rev_mm_dist: final_rev_mm_dist,
        probe_mm_dist: final_probe_mm_dist,
        overall_all_perfect: final_overall_all_perfect,
        overall_max_one_mm: final_overall_max_one_mm,
        overall_two_plus_mm: final_overall_two_plus_mm,
        overall_no_match: final_overall_no_match,
        output_text,
    }
}

/// Generate formatted output text
fn generate_output_text(
    fwd_primers: &[Oligo],
    rev_primers: &[Oligo],
    probes: &[Oligo],
    oligo_stats: &HashMap<String, OligoStats>,
    alignment_dict: &HashMap<String, PatternData>,
    total_sequences: usize,
    sequences_with_min_matches: usize,
    sequences_with_valid_amplicon: usize,
    sequences_failed_amplicon: usize,
    fwd_mm_dist: &MismatchDistribution,
    rev_mm_dist: &MismatchDistribution,
    probe_mm_dist: &MismatchDistribution,
    overall_all_perfect: usize,
    overall_max_one_mm: usize,
    overall_two_plus_mm: usize,
    overall_no_match: usize,
    settings: &AlignmentSettings,
) -> String {
    let mut out = Vec::new();

    out.push("=".repeat(80));
    out.push("qPCR OLIGO INCLUSIVITY ANALYSIS RESULTS".to_string());
    out.push("=".repeat(80));
    out.push(String::new());

    // List oligos by category
    out.push("FORWARD PRIMERS:".to_string());
    for o in fwd_primers {
        out.push(format!("  {} ({})", o.id, o.seq));
    }
    out.push(String::new());

    if !probes.is_empty() {
        out.push("PROBES:".to_string());
        for o in probes {
            out.push(format!("  {} ({})", o.id, o.seq));
        }
        out.push(String::new());
    }

    out.push("REVERSE PRIMERS:".to_string());
    for o in rev_primers {
        out.push(format!("  {} ({})", o.id, o.seq));
    }
    out.push(String::new());

    out.push(format!(
        "Analysis settings: Min fwd matched = {}, Min rev matched = {}, Min probes matched = {}, Min coverage = {}, Max mismatches/oligo = {}",
        settings.min_fwd_matched, settings.min_rev_matched, settings.min_probe_matched, settings.min_coverage, settings.max_mismatches_per_oligo
    ));
    if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
        match (settings.min_amplicon_size, settings.max_amplicon_size) {
            (Some(min), Some(max)) => {
                out.push(format!(
                    "Amplicon size constraint: {} - {} bp",
                    min, max
                ));
            }
            (Some(min), None) => {
                out.push(format!("Amplicon size constraint: >= {} bp", min));
            }
            (None, Some(max)) => {
                out.push(format!("Amplicon size constraint: <= {} bp", max));
            }
            (None, None) => {}
        }
    }
    out.push(String::new());

    out.push("SEQUENCE SIGNATURE PATTERNS:".to_string());
    out.push(format!(
        "Column order: [Forward Primers] || [Probes] || [Reverse Primers]"
    ));
    out.push("-".repeat(50));

    if !alignment_dict.is_empty() {
        let mut sorted_patterns: Vec<_> = alignment_dict.iter().collect();
        sorted_patterns.sort_by(|a, b| b.1.count.cmp(&a.1.count));

        for (signature, data) in sorted_patterns {
            let mut examples_str = data.examples.iter().take(3).cloned().collect::<Vec<_>>();
            if data.examples.len() > 3 {
                examples_str.push(format!("... (+{} more)", data.examples.len() - 3));
            }

            out.push(format!("Pattern: {}", signature));
            out.push(format!(
                "  Count: {}, Mismatches: {}, Fwd matched: {}, Rev matched: {}, Probes matched: {}",
                data.count, data.total_mismatches, data.matched_fwd, data.matched_rev, data.matched_probe
            ));
            if !data.amplicon_lengths.is_empty() {
                let avg_len: f64 =
                    data.amplicon_lengths.iter().sum::<usize>() as f64 / data.amplicon_lengths.len() as f64;
                out.push(format!("  Amplicon length: ~{:.0} bp", avg_len));
            }
            out.push(format!("  Examples: {}", examples_str.join(", ")));
            out.push(String::new());
        }
    } else {
        out.push("No sequences met the minimum matching criteria.".to_string());
        out.push(String::new());
    }

    // Per-oligo statistics by category
    out.push("PER-OLIGO STATISTICS:".to_string());
    out.push("-".repeat(30));

    out.push("  Forward Primers:".to_string());
    for oligo in fwd_primers {
        if let Some(stats) = oligo_stats.get(&oligo.id) {
            let percentage = if total_sequences > 0 {
                (stats.total_matches as f64 / total_sequences as f64) * 100.0
            } else {
                0.0
            };
            out.push(format!(
                "    {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                oligo.id,
                stats.total_matches,
                total_sequences,
                percentage,
                stats.sense_matches,
                stats.antisense_matches
            ));
        }
    }

    if !probes.is_empty() {
        out.push("  Probes:".to_string());
        for oligo in probes {
            if let Some(stats) = oligo_stats.get(&oligo.id) {
                let percentage = if total_sequences > 0 {
                    (stats.total_matches as f64 / total_sequences as f64) * 100.0
                } else {
                    0.0
                };
                out.push(format!(
                    "    {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                    oligo.id,
                    stats.total_matches,
                    total_sequences,
                    percentage,
                    stats.sense_matches,
                    stats.antisense_matches
                ));
            }
        }
    }

    out.push("  Reverse Primers:".to_string());
    for oligo in rev_primers {
        if let Some(stats) = oligo_stats.get(&oligo.id) {
            let percentage = if total_sequences > 0 {
                (stats.total_matches as f64 / total_sequences as f64) * 100.0
            } else {
                0.0
            };
            out.push(format!(
                "    {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                oligo.id,
                stats.total_matches,
                total_sequences,
                percentage,
                stats.sense_matches,
                stats.antisense_matches
            ));
        }
    }

    out.push(String::new());
    out.push("SUMMARY:".to_string());
    out.push("-".repeat(20));
    out.push(format!("Total sequences analyzed: {}", total_sequences));
    let percentage = if total_sequences > 0 {
        (sequences_with_min_matches as f64 / total_sequences as f64) * 100.0
    } else {
        0.0
    };
    out.push(format!(
        "Sequences meeting all thresholds: {} ({:.1}%)",
        sequences_with_min_matches, percentage
    ));

    out.push(String::new());
    out.push("AMPLICON STATISTICS:".to_string());
    out.push("-".repeat(20));
    let amp_percentage = if total_sequences > 0 {
        (sequences_with_valid_amplicon as f64 / total_sequences as f64) * 100.0
    } else {
        0.0
    };
    out.push(format!(
        "Sequences with valid amplicon (fwd+rev pair): {} ({:.1}%)",
        sequences_with_valid_amplicon, amp_percentage
    ));
    out.push(format!(
        "Sequences without valid amplicon: {}",
        sequences_failed_amplicon
    ));
    if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
        let constraint_desc = match (settings.min_amplicon_size, settings.max_amplicon_size) {
            (Some(min), Some(max)) => format!("between {} and {} bp", min, max),
            (Some(min), None) => format!(">= {} bp", min),
            (None, Some(max)) => format!("<= {} bp", max),
            (None, None) => "".to_string(),
        };
        out.push(format!("(Amplicon size constraint: {})", constraint_desc));
    }

    // Mismatch distribution per category
    out.push(String::new());
    out.push("MISMATCH DISTRIBUTION (best oligo per category per sequence):".to_string());
    out.push("-".repeat(50));

    let fmt_pct = |count: usize| -> String {
        if total_sequences > 0 {
            format!("{} ({:.1}%)", count, count as f64 / total_sequences as f64 * 100.0)
        } else {
            format!("{}", count)
        }
    };

    out.push(format!("  Forward Primers:"));
    out.push(format!("    0 mismatches: {}", fmt_pct(fwd_mm_dist.zero_mm)));
    out.push(format!("    1 mismatch:   {}", fmt_pct(fwd_mm_dist.one_mm)));
    out.push(format!("    >1 mismatches: {}", fmt_pct(fwd_mm_dist.more_mm)));
    out.push(format!("    No match:     {}", fmt_pct(fwd_mm_dist.no_match)));

    if !probes.is_empty() {
        out.push(format!("  Probes:"));
        out.push(format!("    0 mismatches: {}", fmt_pct(probe_mm_dist.zero_mm)));
        out.push(format!("    1 mismatch:   {}", fmt_pct(probe_mm_dist.one_mm)));
        out.push(format!("    >1 mismatches: {}", fmt_pct(probe_mm_dist.more_mm)));
        out.push(format!("    No match:     {}", fmt_pct(probe_mm_dist.no_match)));
    }

    out.push(format!("  Reverse Primers:"));
    out.push(format!("    0 mismatches: {}", fmt_pct(rev_mm_dist.zero_mm)));
    out.push(format!("    1 mismatch:   {}", fmt_pct(rev_mm_dist.one_mm)));
    out.push(format!("    >1 mismatches: {}", fmt_pct(rev_mm_dist.more_mm)));
    out.push(format!("    No match:     {}", fmt_pct(rev_mm_dist.no_match)));

    out.push(String::new());
    out.push(format!("  Overall Pattern (worst best-match across all categories):"));
    out.push(format!("    All categories 0 mismatches:  {}", fmt_pct(overall_all_perfect)));
    out.push(format!("    All categories ≤1 mismatch:   {}", fmt_pct(overall_max_one_mm)));
    out.push(format!("    ≥2 mismatches in any category: {}", fmt_pct(overall_two_plus_mm)));
    out.push(format!("    No match in any category:     {}", fmt_pct(overall_no_match)));

    out.push(String::new());
    out.push("LEGEND:".to_string());
    out.push("'(fwd)' = sense orientation, '(rev)' = antisense orientation".to_string());
    out.push("'.' = match, letter = mismatch, '-' = gap/unaligned".to_string());
    out.push("'||' separates oligo categories: Forward Primers || Probes || Reverse Primers".to_string());
    out.push("'|' separates individual oligos within a category".to_string());
    out.push("Probes are only matched within the amplicon region defined by the best fwd+rev pair".to_string());
    out.push("=".repeat(80));

    out.join("\n")
}

// ============================================================================
// GUI Application
// ============================================================================

/// Main application state
struct PrimerAlignApp {
    sequence_file: Option<PathBuf>,
    settings: AlignmentSettings,
    sequences: Vec<FastaRecord>,
    // Categorized oligos
    forward_primers: Vec<Oligo>,
    reverse_primers: Vec<Oligo>,
    probes: Vec<Oligo>,
    fwd_text_input: String,
    rev_text_input: String,
    probe_text_input: String,
    active_oligo_category: OligoCategory,
    // Results and progress
    results: Option<AnalysisResults>,
    progress: ProgressTracker,
    // GUI state
    show_advanced_settings: bool,
    show_results_window: bool,
    show_fasta_info: bool,
    show_about: bool,
    fasta_info_text: String,
    status_message: String,
    analysis_thread: Option<thread::JoinHandle<AnalysisResults>>,
    amplicon_enabled: bool,
    amplicon_min_size_input: String,
    amplicon_max_size_input: String,
}

impl Default for PrimerAlignApp {
    fn default() -> Self {
        Self {
            sequence_file: None,
            settings: AlignmentSettings::default(),
            sequences: Vec::new(),
            forward_primers: Vec::new(),
            reverse_primers: Vec::new(),
            probes: Vec::new(),
            fwd_text_input: String::new(),
            rev_text_input: String::new(),
            probe_text_input: String::new(),
            active_oligo_category: OligoCategory::ForwardPrimer,
            results: None,
            progress: ProgressTracker::default(),
            show_advanced_settings: false,
            show_results_window: false,
            show_fasta_info: false,
            show_about: false,
            fasta_info_text: String::new(),
            status_message: "Ready".to_string(),
            analysis_thread: None,
            amplicon_enabled: false,
            amplicon_min_size_input: String::new(),
            amplicon_max_size_input: "1000".to_string(),
        }
    }
}

impl PrimerAlignApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self::default()
    }

    fn select_sequence_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("FASTA files", &["fasta", "fas", "fa", "fna"])
            .add_filter("All files", &["*"])
            .set_title("Select Sequence File")
            .pick_file()
        {
            match parse_fasta(&path) {
                Ok(records) => {
                    self.sequences = records;
                    self.sequence_file = Some(path);
                    self.status_message = format!("Loaded {} sequences", self.sequences.len());
                }
                Err(e) => {
                    self.status_message = format!("Error loading sequences: {}", e);
                }
            }
        }
    }

    /// Save the current text input to the active oligo category
    fn save_active_category(&mut self) {
        let text = self.get_active_text().to_string();
        let oligos = if text.trim().is_empty() {
            Vec::new()
        } else {
            match parse_fasta_string(&text) {
                Ok(o) => o,
                Err(_) => Vec::new(),
            }
        };

        match self.active_oligo_category {
            OligoCategory::ForwardPrimer => {
                self.fwd_text_input = text;
                self.forward_primers = oligos;
            }
            OligoCategory::ReversePrimer => {
                self.rev_text_input = text;
                self.reverse_primers = oligos;
            }
            OligoCategory::Probe => {
                self.probe_text_input = text;
                self.probes = oligos;
            }
        }
    }

    /// Get the text input for the active category
    fn get_active_text(&self) -> &str {
        match self.active_oligo_category {
            OligoCategory::ForwardPrimer => &self.fwd_text_input,
            OligoCategory::ReversePrimer => &self.rev_text_input,
            OligoCategory::Probe => &self.probe_text_input,
        }
    }

    /// Get mutable reference to the active text input
    fn get_active_text_mut(&mut self) -> &mut String {
        match self.active_oligo_category {
            OligoCategory::ForwardPrimer => &mut self.fwd_text_input,
            OligoCategory::ReversePrimer => &mut self.rev_text_input,
            OligoCategory::Probe => &mut self.probe_text_input,
        }
    }

    /// Check if a category has valid FASTA sequences
    fn category_has_valid_oligos(&self, category: OligoCategory) -> bool {
        match category {
            OligoCategory::ForwardPrimer => !self.forward_primers.is_empty(),
            OligoCategory::ReversePrimer => !self.reverse_primers.is_empty(),
            OligoCategory::Probe => !self.probes.is_empty(),
        }
    }

    /// Parse all oligo categories from their text inputs
    fn parse_all_oligos(&mut self) {
        if !self.fwd_text_input.trim().is_empty() {
            if let Ok(oligos) = parse_fasta_string(&self.fwd_text_input) {
                self.forward_primers = oligos;
            }
        } else {
            self.forward_primers.clear();
        }

        if !self.rev_text_input.trim().is_empty() {
            if let Ok(oligos) = parse_fasta_string(&self.rev_text_input) {
                self.reverse_primers = oligos;
            }
        } else {
            self.reverse_primers.clear();
        }

        if !self.probe_text_input.trim().is_empty() {
            if let Ok(oligos) = parse_fasta_string(&self.probe_text_input) {
                self.probes = oligos;
            }
        } else {
            self.probes.clear();
        }
    }

    /// Save all oligo categories to a JSON file
    fn save_oligos_json(&mut self) {
        self.parse_all_oligos();

        let oligo_set = OligoSet {
            forward_primers: self
                .forward_primers
                .iter()
                .map(|o| OligoSerde {
                    id: o.id.clone(),
                    seq: o.seq.clone(),
                })
                .collect(),
            reverse_primers: self
                .reverse_primers
                .iter()
                .map(|o| OligoSerde {
                    id: o.id.clone(),
                    seq: o.seq.clone(),
                })
                .collect(),
            probes: self
                .probes
                .iter()
                .map(|o| OligoSerde {
                    id: o.id.clone(),
                    seq: o.seq.clone(),
                })
                .collect(),
        };

        if let Some(path) = FileDialog::new()
            .add_filter("JSON files", &["json"])
            .set_title("Save Oligo Set")
            .save_file()
        {
            match serde_json::to_string_pretty(&oligo_set) {
                Ok(json) => match std::fs::write(&path, json) {
                    Ok(_) => {
                        self.status_message = format!(
                            "Oligo set saved to {}",
                            path.file_name()
                                .map(|s| s.to_string_lossy().to_string())
                                .unwrap_or_default()
                        );
                    }
                    Err(e) => {
                        self.status_message = format!("Error saving oligo set: {}", e);
                    }
                },
                Err(e) => {
                    self.status_message = format!("Error serializing oligo set: {}", e);
                }
            }
        }
    }

    /// Load oligo categories from a JSON file
    fn load_oligos_json(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("JSON files", &["json"])
            .set_title("Load Oligo Set")
            .pick_file()
        {
            match std::fs::read_to_string(&path) {
                Ok(content) => match serde_json::from_str::<OligoSet>(&content) {
                    Ok(oligo_set) => {
                        self.forward_primers = oligo_set
                            .forward_primers
                            .iter()
                            .map(|o| Oligo {
                                id: o.id.clone(),
                                seq: o.seq.clone(),
                            })
                            .collect();
                        self.reverse_primers = oligo_set
                            .reverse_primers
                            .iter()
                            .map(|o| Oligo {
                                id: o.id.clone(),
                                seq: o.seq.clone(),
                            })
                            .collect();
                        self.probes = oligo_set
                            .probes
                            .iter()
                            .map(|o| Oligo {
                                id: o.id.clone(),
                                seq: o.seq.clone(),
                            })
                            .collect();

                        // Update text inputs
                        self.fwd_text_input = oligos_to_fasta_string(&self.forward_primers);
                        self.rev_text_input = oligos_to_fasta_string(&self.reverse_primers);
                        self.probe_text_input = oligos_to_fasta_string(&self.probes);

                        self.status_message = format!(
                            "Loaded {} fwd, {} rev, {} probes from {}",
                            self.forward_primers.len(),
                            self.reverse_primers.len(),
                            self.probes.len(),
                            path.file_name()
                                .map(|s| s.to_string_lossy().to_string())
                                .unwrap_or_default()
                        );
                    }
                    Err(e) => {
                        self.status_message = format!("Error parsing JSON: {}", e);
                    }
                },
                Err(e) => {
                    self.status_message = format!("Error reading file: {}", e);
                }
            }
        }
    }

    fn update_amplicon_setting(&mut self) {
        if self.amplicon_enabled {
            if self.amplicon_min_size_input.trim().is_empty() {
                self.settings.min_amplicon_size = None;
            } else if let Ok(size) = self.amplicon_min_size_input.parse::<usize>() {
                if size > 0 {
                    self.settings.min_amplicon_size = Some(size);
                } else {
                    self.settings.min_amplicon_size = None;
                }
            } else {
                self.settings.min_amplicon_size = None;
            }

            if self.amplicon_max_size_input.trim().is_empty() {
                self.settings.max_amplicon_size = None;
            } else if let Ok(size) = self.amplicon_max_size_input.parse::<usize>() {
                if size > 0 {
                    self.settings.max_amplicon_size = Some(size);
                } else {
                    self.settings.max_amplicon_size = None;
                }
            } else {
                self.settings.max_amplicon_size = None;
            }
        } else {
            self.settings.min_amplicon_size = None;
            self.settings.max_amplicon_size = None;
        }
    }

    fn can_run_analysis(&self) -> bool {
        !self.forward_primers.is_empty() && !self.reverse_primers.is_empty()
    }

    fn start_analysis(&mut self) {
        if self.sequences.is_empty() {
            self.status_message = "Please select a sequence file first".to_string();
            return;
        }

        // Parse all oligo categories
        self.parse_all_oligos();

        if self.forward_primers.is_empty() {
            self.status_message =
                "Please enter at least one forward primer".to_string();
            return;
        }

        if self.reverse_primers.is_empty() {
            self.status_message =
                "Please enter at least one reverse primer".to_string();
            return;
        }

        // Update amplicon setting before starting
        self.update_amplicon_setting();

        let sequences = self.sequences.clone();
        let fwd_primers = self.forward_primers.clone();
        let rev_primers = self.reverse_primers.clone();
        let probes = self.probes.clone();
        let settings = self.settings.clone();
        let progress = self.progress.clone();

        self.progress.running.store(true, Ordering::SeqCst);
        self.progress.current.store(0, Ordering::SeqCst);
        self.progress
            .total
            .store(sequences.len(), Ordering::SeqCst);
        if let Ok(mut status) = self.progress.status.lock() {
            *status = "Starting analysis...".to_string();
        }

        let handle = thread::spawn(move || {
            run_analysis(
                &sequences,
                &fwd_primers,
                &rev_primers,
                &probes,
                &settings,
                &progress,
            )
        });

        self.analysis_thread = Some(handle);
    }

    fn check_analysis_complete(&mut self) {
        if let Some(handle) = self.analysis_thread.take() {
            if handle.is_finished() {
                match handle.join() {
                    Ok(results) => {
                        self.results = Some(results);
                        self.progress.running.store(false, Ordering::SeqCst);
                        self.status_message = "Analysis complete!".to_string();
                        self.show_results_window = true;
                    }
                    Err(_) => {
                        self.progress.running.store(false, Ordering::SeqCst);
                        self.status_message = "Analysis failed".to_string();
                    }
                }
            } else {
                self.analysis_thread = Some(handle);
            }
        }
    }

    fn show_fasta_info(&mut self) {
        if self.sequences.is_empty() {
            self.status_message = "Please select a sequence file first".to_string();
            return;
        }

        let mut info = format!("FASTA File Information:\n{}\n\n", "=".repeat(30));
        info.push_str(&format!("Total sequences: {}\n", self.sequences.len()));

        if let Some(ref path) = self.sequence_file {
            info.push_str(&format!(
                "File: {}\n\n",
                path.file_name()
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_default()
            ));
        }

        if !self.sequences.is_empty() {
            let seq_lengths: Vec<usize> = self.sequences.iter().map(|r| r.seq.len()).collect();
            let min_len = seq_lengths.iter().min().unwrap_or(&0);
            let max_len = seq_lengths.iter().max().unwrap_or(&0);
            let avg_len: f64 =
                seq_lengths.iter().sum::<usize>() as f64 / seq_lengths.len() as f64;

            info.push_str("Sequence length statistics:\n");
            info.push_str(&format!("  Min: {} bp\n", min_len));
            info.push_str(&format!("  Max: {} bp\n", max_len));
            info.push_str(&format!("  Average: {:.1} bp\n\n", avg_len));

            info.push_str("First 5 sequence IDs:\n");
            for (i, record) in self.sequences.iter().take(5).enumerate() {
                info.push_str(&format!("  {}. {}\n", i + 1, record.id));
            }
        }

        self.fasta_info_text = info;
        self.show_fasta_info = true;
    }

    fn export_to_excel(&mut self) {
        let results = match &self.results {
            Some(r) => r,
            None => {
                self.status_message = "No results to export. Run analysis first.".to_string();
                return;
            }
        };

        if let Some(path) = FileDialog::new()
            .add_filter("Excel files", &["xlsx"])
            .set_title("Save Excel Results As")
            .save_file()
        {
            match self.write_excel(&path, results) {
                Ok(_) => {
                    self.status_message = format!(
                        "Excel export complete: {}",
                        path.file_name()
                            .map(|s| s.to_string_lossy().to_string())
                            .unwrap_or_default()
                    );
                }
                Err(e) => {
                    self.status_message = format!("Error exporting Excel: {}", e);
                }
            }
        }
    }

    fn write_excel(&self, path: &PathBuf, results: &AnalysisResults) -> Result<(), String> {
        let mut workbook = Workbook::new();
        let worksheet = workbook.add_worksheet();

        let header_format = Format::new().set_bold();
        let title_format = Format::new().set_bold().set_font_size(14);
        let category_format = Format::new().set_bold().set_font_size(11);

        let mut row: u32 = 0;

        // Title
        worksheet
            .write_string_with_format(
                row,
                0,
                "qPCR Oligo Inclusivity Analysis Results",
                &title_format,
            )
            .map_err(|e| e.to_string())?;
        row += 2;

        // Settings
        worksheet
            .write_string(
                row,
                0,
                &format!(
                    "Settings: Min fwd = {}, Min rev = {}, Min probes = {}, Min coverage = {}, Max mm/oligo = {}",
                    self.settings.min_fwd_matched,
                    self.settings.min_rev_matched,
                    self.settings.min_probe_matched,
                    self.settings.min_coverage,
                    self.settings.max_mismatches_per_oligo
                ),
            )
            .map_err(|e| e.to_string())?;
        row += 1;

        if self.settings.min_amplicon_size.is_some() || self.settings.max_amplicon_size.is_some() {
            let constraint_text = match (
                self.settings.min_amplicon_size,
                self.settings.max_amplicon_size,
            ) {
                (Some(min), Some(max)) => {
                    format!("Amplicon size constraint: {} - {} bp", min, max)
                }
                (Some(min), None) => format!("Amplicon size constraint: >= {} bp", min),
                (None, Some(max)) => format!("Amplicon size constraint: <= {} bp", max),
                (None, None) => "".to_string(),
            };
            worksheet
                .write_string(row, 0, &constraint_text)
                .map_err(|e| e.to_string())?;
            row += 1;
        }
        row += 2;

        // Category labels row
        let mut col: u16 = 1; // skip pattern # column
        if !self.forward_primers.is_empty() {
            worksheet
                .write_string_with_format(row, col, "--- Forward Primers ---", &category_format)
                .map_err(|e| e.to_string())?;
        }
        col += self.forward_primers.len() as u16;
        if !self.probes.is_empty() {
            worksheet
                .write_string_with_format(row, col, "--- Probes ---", &category_format)
                .map_err(|e| e.to_string())?;
        }
        col += self.probes.len() as u16;
        if !self.reverse_primers.is_empty() {
            worksheet
                .write_string_with_format(row, col, "--- Reverse Primers ---", &category_format)
                .map_err(|e| e.to_string())?;
        }
        row += 1;

        // Headers row
        col = 0;
        worksheet
            .write_string_with_format(row, col, "Pattern #", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;

        // Column headers: FWD | PROBES | REV
        let ordered_oligos: Vec<&Oligo> = self
            .forward_primers
            .iter()
            .chain(self.probes.iter())
            .chain(self.reverse_primers.iter())
            .collect();

        for oligo in &ordered_oligos {
            worksheet
                .write_string_with_format(
                    row,
                    col,
                    &format!("{} Pattern", oligo.id),
                    &header_format,
                )
                .map_err(|e| e.to_string())?;
            col += 1;
        }

        worksheet
            .write_string_with_format(row, col, "Count", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;
        worksheet
            .write_string_with_format(row, col, "Percentage", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;
        worksheet
            .write_string_with_format(row, col, "Total Mismatches", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;
        worksheet
            .write_string_with_format(row, col, "Fwd Matched", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;
        worksheet
            .write_string_with_format(row, col, "Rev Matched", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;
        worksheet
            .write_string_with_format(row, col, "Probes Matched", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;
        worksheet
            .write_string_with_format(row, col, "Amplicon Length", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;
        worksheet
            .write_string_with_format(row, col, "Example Sequences", &header_format)
            .map_err(|e| e.to_string())?;
        row += 1;

        // Row with full oligo sequences
        col = 1;
        for oligo in &ordered_oligos {
            worksheet
                .write_string(row, col, &oligo.seq)
                .map_err(|e| e.to_string())?;
            col += 1;
        }
        row += 1;

        // Data rows
        let mut sorted_patterns: Vec<_> = results.alignment_dict.iter().collect();
        sorted_patterns.sort_by(|a, b| b.1.count.cmp(&a.1.count));

        let num_fwd = self.forward_primers.len();
        let num_probes = self.probes.len();
        let num_rev = self.reverse_primers.len();

        let mut pattern_num = 1u32;
        for (signature, data) in sorted_patterns {
            col = 0;
            worksheet
                .write_number(row, col, pattern_num as f64)
                .map_err(|e| e.to_string())?;
            col += 1;

            // Parse the combined signature: "fwd1 | fwd2 || probe1 | probe2 || rev1 | rev2"
            let category_sections: Vec<&str> = signature.split(" || ").collect();

            // Extract individual patterns from each section
            let mut all_patterns: Vec<String> = Vec::new();

            // Forward primers section (index 0)
            let fwd_section = category_sections.first().unwrap_or(&"");
            let fwd_patterns: Vec<&str> = if !fwd_section.is_empty() {
                fwd_section.split(" | ").collect()
            } else {
                Vec::new()
            };
            for i in 0..num_fwd {
                all_patterns.push(
                    fwd_patterns
                        .get(i)
                        .unwrap_or(&"NO_MATCH")
                        .to_string(),
                );
            }

            // Probes section (index 1 if probes exist, otherwise skip)
            if num_probes > 0 {
                let probe_section = if category_sections.len() >= 3 {
                    category_sections.get(1).unwrap_or(&"")
                } else {
                    &""
                };
                let probe_patterns: Vec<&str> = if !probe_section.is_empty() {
                    probe_section.split(" | ").collect()
                } else {
                    Vec::new()
                };
                for i in 0..num_probes {
                    all_patterns.push(
                        probe_patterns
                            .get(i)
                            .unwrap_or(&"NO_MATCH")
                            .to_string(),
                    );
                }
            }

            // Reverse primers section (last section)
            let rev_section_idx = if num_probes > 0 { 2 } else { 1 };
            let rev_section = category_sections.get(rev_section_idx).unwrap_or(&"");
            let rev_patterns: Vec<&str> = if !rev_section.is_empty() {
                rev_section.split(" | ").collect()
            } else {
                Vec::new()
            };
            for i in 0..num_rev {
                all_patterns.push(
                    rev_patterns
                        .get(i)
                        .unwrap_or(&"NO_MATCH")
                        .to_string(),
                );
            }

            // Write patterns
            for pattern in &all_patterns {
                worksheet
                    .write_string(row, col, pattern)
                    .map_err(|e| e.to_string())?;
                col += 1;
            }

            // Count
            worksheet
                .write_number(row, col, data.count as f64)
                .map_err(|e| e.to_string())?;
            col += 1;

            // Percentage
            let percentage = if results.total_sequences > 0 {
                (data.count as f64 / results.total_sequences as f64) * 100.0
            } else {
                0.0
            };
            worksheet
                .write_number(row, col, percentage)
                .map_err(|e| e.to_string())?;
            col += 1;

            // Total Mismatches
            worksheet
                .write_number(row, col, data.total_mismatches as f64)
                .map_err(|e| e.to_string())?;
            col += 1;

            // Fwd Matched
            worksheet
                .write_number(row, col, data.matched_fwd as f64)
                .map_err(|e| e.to_string())?;
            col += 1;

            // Rev Matched
            worksheet
                .write_number(row, col, data.matched_rev as f64)
                .map_err(|e| e.to_string())?;
            col += 1;

            // Probes Matched
            worksheet
                .write_number(row, col, data.matched_probe as f64)
                .map_err(|e| e.to_string())?;
            col += 1;

            // Amplicon Length - show most common length, or empty if none found
            if !data.amplicon_lengths.is_empty() {
                use std::collections::HashMap;
                let mut length_counts: HashMap<usize, usize> = HashMap::new();
                for &len in &data.amplicon_lengths {
                    *length_counts.entry(len).or_insert(0) += 1;
                }
                let amplicon_length = length_counts
                    .iter()
                    .max_by_key(|(_, &count)| count)
                    .map(|(&len, _)| len)
                    .unwrap_or(0);
                worksheet
                    .write_number(row, col, amplicon_length as f64)
                    .map_err(|e| e.to_string())?;
            } else {
                worksheet
                    .write_string(row, col, "")
                    .map_err(|e| e.to_string())?;
            }
            col += 1;

            // Examples
            let examples: Vec<_> = data.examples.iter().take(3).collect();
            let mut examples_str = examples
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(", ");
            if data.examples.len() > 3 {
                examples_str.push_str(&format!(" (+{} more)", data.examples.len() - 3));
            }
            worksheet
                .write_string(row, col, &examples_str)
                .map_err(|e| e.to_string())?;

            row += 1;
            pattern_num += 1;
        }

        // Per-oligo Statistics
        row += 2;
        worksheet
            .write_string_with_format(row, 0, "PER-OLIGO STATISTICS:", &header_format)
            .map_err(|e| e.to_string())?;
        row += 1;

        worksheet
            .write_string_with_format(row, 0, "Forward Primers:", &category_format)
            .map_err(|e| e.to_string())?;
        row += 1;
        for oligo in &self.forward_primers {
            if let Some(stats) = results.oligo_stats.get(&oligo.id) {
                let percentage = if results.total_sequences > 0 {
                    (stats.total_matches as f64 / results.total_sequences as f64) * 100.0
                } else {
                    0.0
                };
                worksheet
                    .write_string(
                        row,
                        0,
                        &format!(
                            "  {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                            oligo.id,
                            stats.total_matches,
                            results.total_sequences,
                            percentage,
                            stats.sense_matches,
                            stats.antisense_matches
                        ),
                    )
                    .map_err(|e| e.to_string())?;
                row += 1;
            }
        }

        if !self.probes.is_empty() {
            worksheet
                .write_string_with_format(row, 0, "Probes:", &category_format)
                .map_err(|e| e.to_string())?;
            row += 1;
            for oligo in &self.probes {
                if let Some(stats) = results.oligo_stats.get(&oligo.id) {
                    let percentage = if results.total_sequences > 0 {
                        (stats.total_matches as f64 / results.total_sequences as f64) * 100.0
                    } else {
                        0.0
                    };
                    worksheet
                        .write_string(
                            row,
                            0,
                            &format!(
                                "  {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                                oligo.id,
                                stats.total_matches,
                                results.total_sequences,
                                percentage,
                                stats.sense_matches,
                                stats.antisense_matches
                            ),
                        )
                        .map_err(|e| e.to_string())?;
                    row += 1;
                }
            }
        }

        worksheet
            .write_string_with_format(row, 0, "Reverse Primers:", &category_format)
            .map_err(|e| e.to_string())?;
        row += 1;
        for oligo in &self.reverse_primers {
            if let Some(stats) = results.oligo_stats.get(&oligo.id) {
                let percentage = if results.total_sequences > 0 {
                    (stats.total_matches as f64 / results.total_sequences as f64) * 100.0
                } else {
                    0.0
                };
                worksheet
                    .write_string(
                        row,
                        0,
                        &format!(
                            "  {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                            oligo.id,
                            stats.total_matches,
                            results.total_sequences,
                            percentage,
                            stats.sense_matches,
                            stats.antisense_matches
                        ),
                    )
                    .map_err(|e| e.to_string())?;
                row += 1;
            }
        }

        // Summary
        row += 1;
        worksheet
            .write_string_with_format(row, 0, "SUMMARY:", &header_format)
            .map_err(|e| e.to_string())?;
        row += 1;
        worksheet
            .write_string(
                row,
                0,
                &format!("Total sequences analyzed: {}", results.total_sequences),
            )
            .map_err(|e| e.to_string())?;
        row += 1;

        let percentage = if results.total_sequences > 0 {
            (results.sequences_with_min_matches as f64 / results.total_sequences as f64) * 100.0
        } else {
            0.0
        };
        worksheet
            .write_string(
                row,
                0,
                &format!(
                    "Sequences meeting all thresholds: {} ({:.1}%)",
                    results.sequences_with_min_matches, percentage
                ),
            )
            .map_err(|e| e.to_string())?;
        row += 1;

        // Amplicon statistics
        row += 1;
        worksheet
            .write_string_with_format(row, 0, "AMPLICON STATISTICS:", &header_format)
            .map_err(|e| e.to_string())?;
        row += 1;

        let amp_percentage = if results.total_sequences > 0 {
            (results.sequences_with_valid_amplicon as f64 / results.total_sequences as f64) * 100.0
        } else {
            0.0
        };
        worksheet
            .write_string(
                row,
                0,
                &format!(
                    "Sequences with valid amplicon: {} ({:.1}%)",
                    results.sequences_with_valid_amplicon, amp_percentage
                ),
            )
            .map_err(|e| e.to_string())?;
        row += 1;

        worksheet
            .write_string(
                row,
                0,
                &format!(
                    "Sequences without valid amplicon: {}",
                    results.sequences_failed_amplicon
                ),
            )
            .map_err(|e| e.to_string())?;

        // Mismatch distribution per category
        row += 2;
        worksheet
            .write_string_with_format(
                row,
                0,
                "MISMATCH DISTRIBUTION (best oligo per category per sequence):",
                &header_format,
            )
            .map_err(|e| e.to_string())?;
        row += 1;

        let fmt_pct_xl = |count: usize| -> String {
            if results.total_sequences > 0 {
                format!(
                    "{} ({:.1}%)",
                    count,
                    count as f64 / results.total_sequences as f64 * 100.0
                )
            } else {
                format!("{}", count)
            }
        };

        worksheet
            .write_string_with_format(row, 0, "Forward Primers:", &category_format)
            .map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  0 mismatches: {}", fmt_pct_xl(results.fwd_mm_dist.zero_mm))).map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  1 mismatch: {}", fmt_pct_xl(results.fwd_mm_dist.one_mm))).map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  >1 mismatches: {}", fmt_pct_xl(results.fwd_mm_dist.more_mm))).map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  No match: {}", fmt_pct_xl(results.fwd_mm_dist.no_match))).map_err(|e| e.to_string())?;
        row += 1;

        if !self.probes.is_empty() {
            worksheet
                .write_string_with_format(row, 0, "Probes:", &category_format)
                .map_err(|e| e.to_string())?;
            row += 1;
            worksheet.write_string(row, 0, &format!("  0 mismatches: {}", fmt_pct_xl(results.probe_mm_dist.zero_mm))).map_err(|e| e.to_string())?;
            row += 1;
            worksheet.write_string(row, 0, &format!("  1 mismatch: {}", fmt_pct_xl(results.probe_mm_dist.one_mm))).map_err(|e| e.to_string())?;
            row += 1;
            worksheet.write_string(row, 0, &format!("  >1 mismatches: {}", fmt_pct_xl(results.probe_mm_dist.more_mm))).map_err(|e| e.to_string())?;
            row += 1;
            worksheet.write_string(row, 0, &format!("  No match: {}", fmt_pct_xl(results.probe_mm_dist.no_match))).map_err(|e| e.to_string())?;
            row += 1;
        }

        worksheet
            .write_string_with_format(row, 0, "Reverse Primers:", &category_format)
            .map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  0 mismatches: {}", fmt_pct_xl(results.rev_mm_dist.zero_mm))).map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  1 mismatch: {}", fmt_pct_xl(results.rev_mm_dist.one_mm))).map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  >1 mismatches: {}", fmt_pct_xl(results.rev_mm_dist.more_mm))).map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  No match: {}", fmt_pct_xl(results.rev_mm_dist.no_match))).map_err(|e| e.to_string())?;
        row += 1;

        // Overall pattern distribution
        row += 1;
        worksheet
            .write_string_with_format(row, 0, "Overall Pattern (worst best-match across all categories):", &category_format)
            .map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  All categories 0 mismatches: {}", fmt_pct_xl(results.overall_all_perfect))).map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  All categories ≤1 mismatch: {}", fmt_pct_xl(results.overall_max_one_mm))).map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  ≥2 mismatches in any category: {}", fmt_pct_xl(results.overall_two_plus_mm))).map_err(|e| e.to_string())?;
        row += 1;
        worksheet.write_string(row, 0, &format!("  No match in any category: {}", fmt_pct_xl(results.overall_no_match))).map_err(|e| e.to_string())?;

        workbook.save(path).map_err(|e| e.to_string())?;

        Ok(())
    }

    fn reset_defaults(&mut self) {
        self.settings = AlignmentSettings::default();
        self.amplicon_enabled = false;
        self.amplicon_min_size_input = String::new();
        self.amplicon_max_size_input = "1000".to_string();
    }
}

impl eframe::App for PrimerAlignApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.check_analysis_complete();

        if self.progress.running.load(Ordering::SeqCst) {
            ctx.request_repaint();
        }

        // Menu bar
        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Select Sequence File...").clicked() {
                        self.select_sequence_file();
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Export to Excel...").clicked() {
                        self.export_to_excel();
                        ui.close_menu();
                    }
                });
                ui.menu_button("Tools", |ui| {
                    if ui.button("FASTA Info").clicked() {
                        self.show_fasta_info();
                        ui.close_menu();
                    }
                    if ui.button("Advanced Settings...").clicked() {
                        self.show_advanced_settings = true;
                        ui.close_menu();
                    }
                });
                ui.menu_button("Help", |ui| {
                    if ui.button("About").clicked() {
                        self.show_about = true;
                        ui.close_menu();
                    }
                });
            });
        });

        // Status bar
        egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                let status = if self.progress.running.load(Ordering::SeqCst) {
                    self.progress
                        .status
                        .lock()
                        .map(|s| s.clone())
                        .unwrap_or_else(|_| "Processing...".to_string())
                } else {
                    self.status_message.clone()
                };
                ui.label(&status);
            });
        });

        // Main content
        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.heading("qPCR Oligo Inclusivity Tool");
                ui.label("Categorized oligo analysis with amplicon detection");
                ui.add_space(10.0);

                // File selection group
                ui.group(|ui| {
                    ui.label(egui::RichText::new("Sequence File").strong());
                    ui.add_space(5.0);

                    ui.horizontal(|ui| {
                        ui.label("Sequence File:");
                        if ui.button("Browse...").clicked() {
                            self.select_sequence_file();
                        }
                        if let Some(ref path) = self.sequence_file {
                            ui.colored_label(
                                egui::Color32::DARK_GREEN,
                                format!(
                                    "{} ({} sequences)",
                                    path.file_name()
                                        .map(|s| s.to_string_lossy().to_string())
                                        .unwrap_or_default(),
                                    self.sequences.len()
                                ),
                            );
                        } else {
                            ui.label("No file selected");
                        }
                    });
                });

                ui.add_space(10.0);

                // Oligo sequences input group with category tabs
                ui.group(|ui| {
                    ui.label(egui::RichText::new("Oligo Sequences").strong());
                    ui.add_space(5.0);

                    // Category buttons and Save/Load
                    ui.horizontal(|ui| {
                        // Forward Primers button
                        let fwd_color = if self.category_has_valid_oligos(OligoCategory::ForwardPrimer) {
                            egui::Color32::DARK_GREEN
                        } else {
                            egui::Color32::from_rgb(180, 40, 40)
                        };
                        let fwd_selected = self.active_oligo_category == OligoCategory::ForwardPrimer;
                        let fwd_text = egui::RichText::new("Forward Primers")
                            .color(if fwd_selected { egui::Color32::WHITE } else { fwd_color });
                        let fwd_btn = if fwd_selected {
                            ui.add(egui::Button::new(fwd_text).fill(fwd_color))
                        } else {
                            ui.add(egui::Button::new(fwd_text))
                        };
                        if fwd_btn.clicked() && !fwd_selected {
                            self.save_active_category();
                            self.active_oligo_category = OligoCategory::ForwardPrimer;
                        }

                        // Reverse Primers button
                        let rev_color = if self.category_has_valid_oligos(OligoCategory::ReversePrimer) {
                            egui::Color32::DARK_GREEN
                        } else {
                            egui::Color32::from_rgb(180, 40, 40)
                        };
                        let rev_selected = self.active_oligo_category == OligoCategory::ReversePrimer;
                        let rev_text = egui::RichText::new("Reverse Primers")
                            .color(if rev_selected { egui::Color32::WHITE } else { rev_color });
                        let rev_btn = if rev_selected {
                            ui.add(egui::Button::new(rev_text).fill(rev_color))
                        } else {
                            ui.add(egui::Button::new(rev_text))
                        };
                        if rev_btn.clicked() && !rev_selected {
                            self.save_active_category();
                            self.active_oligo_category = OligoCategory::ReversePrimer;
                        }

                        // Probes button
                        let probe_color = if self.category_has_valid_oligos(OligoCategory::Probe) {
                            egui::Color32::DARK_GREEN
                        } else {
                            egui::Color32::from_rgb(180, 40, 40)
                        };
                        let probe_selected = self.active_oligo_category == OligoCategory::Probe;
                        let probe_text = egui::RichText::new("Probes")
                            .color(if probe_selected { egui::Color32::WHITE } else { probe_color });
                        let probe_btn = if probe_selected {
                            ui.add(egui::Button::new(probe_text).fill(probe_color))
                        } else {
                            ui.add(egui::Button::new(probe_text))
                        };
                        if probe_btn.clicked() && !probe_selected {
                            self.save_active_category();
                            self.active_oligo_category = OligoCategory::Probe;
                        }

                        ui.add_space(20.0);

                        if ui.button("Save JSON").clicked() {
                            self.save_active_category();
                            self.save_oligos_json();
                        }
                        if ui.button("Load JSON").clicked() {
                            self.load_oligos_json();
                        }
                    });

                    ui.add_space(5.0);

                    // Category label
                    let category_label = match self.active_oligo_category {
                        OligoCategory::ForwardPrimer => "Forward Primers (FASTA format)",
                        OligoCategory::ReversePrimer => "Reverse Primers (FASTA format)",
                        OligoCategory::Probe => "Probes (FASTA format)",
                    };
                    ui.label(
                        egui::RichText::new(category_label)
                            .small()
                            .color(egui::Color32::GRAY),
                    );

                    // Text area for active category
                    let text = self.get_active_text_mut();
                    egui::ScrollArea::vertical()
                        .id_salt("oligo_text_scroll")
                        .max_height(120.0)
                        .show(ui, |ui| {
                            ui.add(
                                egui::TextEdit::multiline(text)
                                    .font(egui::TextStyle::Monospace)
                                    .desired_width(f32::INFINITY)
                                    .desired_rows(6)
                                    .hint_text(">Primer1\nATGCGTACGTAGC\n>Primer2\nGCTAGCTAGCTA"),
                            );
                        });

                    ui.add_space(3.0);

                    // Show count for active category
                    let active_text = self.get_active_text();
                    if !active_text.is_empty() {
                        let count = active_text
                            .lines()
                            .filter(|line| line.starts_with('>'))
                            .count();
                        if count > 0 {
                            ui.colored_label(
                                egui::Color32::DARK_GREEN,
                                format!("{} sequence(s) defined", count),
                            );
                        }
                    }
                });

                ui.add_space(10.0);

                // Quick settings group
                ui.group(|ui| {
                    ui.label(egui::RichText::new("Quick Settings").strong());
                    ui.add_space(5.0);

                    ui.horizontal(|ui| {
                        ui.label("Min Fwd Matched:");
                        ui.add(
                            egui::DragValue::new(&mut self.settings.min_fwd_matched)
                                .range(0..=100)
                                .speed(0.1),
                        );

                        ui.add_space(10.0);

                        ui.label("Min Rev Matched:");
                        ui.add(
                            egui::DragValue::new(&mut self.settings.min_rev_matched)
                                .range(0..=100)
                                .speed(0.1),
                        );

                        ui.add_space(10.0);

                        ui.label("Min Probes Matched:");
                        ui.add(
                            egui::DragValue::new(&mut self.settings.min_probe_matched)
                                .range(0..=100)
                                .speed(0.1),
                        );
                    });

                    ui.add_space(5.0);

                    ui.horizontal(|ui| {
                        ui.label("Min Coverage:");
                        ui.add(
                            egui::DragValue::new(&mut self.settings.min_coverage)
                                .range(0.0..=1.0)
                                .speed(0.01)
                                .max_decimals(2),
                        );

                        ui.add_space(20.0);

                        ui.label("Max Mismatches/Oligo:");
                        ui.add(
                            egui::DragValue::new(&mut self.settings.max_mismatches_per_oligo)
                                .range(0..=50)
                                .speed(0.1),
                        );

                        ui.add_space(20.0);

                        ui.label("Ambiguity Match Display:");
                        egui::ComboBox::from_id_salt("ambiguity_display")
                            .selected_text(match self.settings.ambiguity_display {
                                AmbiguityDisplayMode::ShowBases => "Show Base Letters",
                                AmbiguityDisplayMode::ShowDots => "Show Dots",
                            })
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut self.settings.ambiguity_display,
                                    AmbiguityDisplayMode::ShowDots,
                                    "Show Dots",
                                );
                                ui.selectable_value(
                                    &mut self.settings.ambiguity_display,
                                    AmbiguityDisplayMode::ShowBases,
                                    "Show Base Letters",
                                );
                            });

                        ui.add_space(20.0);

                        if ui.button("Advanced Settings...").clicked() {
                            self.show_advanced_settings = true;
                        }
                    });
                });

                ui.add_space(10.0);

                // Amplicon constraint settings
                ui.group(|ui| {
                    ui.label(egui::RichText::new("Amplicon Size Constraints").strong());
                    ui.add_space(5.0);

                    ui.horizontal(|ui| {
                        ui.checkbox(
                            &mut self.amplicon_enabled,
                            "Enable amplicon size limits",
                        );
                    });

                    if self.amplicon_enabled {
                        ui.add_space(5.0);
                        ui.add_enabled_ui(self.amplicon_enabled, |ui| {
                            ui.horizontal(|ui| {
                                ui.label("Min amplicon size (bp):");
                                ui.add(
                                    egui::TextEdit::singleline(&mut self.amplicon_min_size_input)
                                        .desired_width(80.0),
                                );

                                ui.add_space(10.0);

                                ui.label("Max amplicon size (bp):");
                                ui.add(
                                    egui::TextEdit::singleline(&mut self.amplicon_max_size_input)
                                        .desired_width(80.0),
                                );
                            });
                        });

                        ui.add_space(5.0);
                        ui.label(
                            egui::RichText::new(
                                "Limits the amplicon size for valid fwd+rev primer pairs. Without size limits, any convergent pair is accepted."
                            )
                            .small()
                            .color(egui::Color32::GRAY),
                        );
                    } else {
                        ui.add_space(3.0);
                        ui.label(
                            egui::RichText::new(
                                "Amplicon detection is always active (fwd+rev pairing). This only constrains the allowed size range."
                            )
                            .small()
                            .color(egui::Color32::GRAY),
                        );
                    }
                });

                ui.add_space(15.0);

                // Action buttons
                ui.horizontal(|ui| {
                    let is_running = self.progress.running.load(Ordering::SeqCst);

                    // Save active category before checking if we can run
                    // (We check forward_primers/reverse_primers which are updated on tab switch)

                    let can_run = self.can_run_analysis() && !self.sequences.is_empty();

                    ui.add_enabled_ui(!is_running, |ui| {
                        let btn_color = if can_run {
                            egui::Color32::DARK_GREEN
                        } else {
                            egui::Color32::from_rgb(180, 40, 40)
                        };
                        let btn = ui.add(
                            egui::Button::new(
                                egui::RichText::new("Run Analysis").strong().color(egui::Color32::WHITE),
                            )
                            .fill(btn_color),
                        );
                        if btn.clicked() {
                            self.save_active_category();
                            self.start_analysis();
                        }
                    });

                    if ui.button("FASTA Info").clicked() {
                        self.show_fasta_info();
                    }

                    if is_running {
                        ui.add_space(10.0);
                        if ui.button("Stop").clicked() {
                            self.progress.running.store(false, Ordering::SeqCst);
                        }
                        ui.spinner();
                    }
                });

                ui.add_space(15.0);

                // Progress group
                ui.group(|ui| {
                    ui.label(egui::RichText::new("Progress").strong());
                    let current = self.progress.current.load(Ordering::SeqCst);
                    let total = self.progress.total.load(Ordering::SeqCst);
                    let progress = if total > 0 {
                        current as f32 / total as f32
                    } else {
                        0.0
                    };

                    ui.add(egui::ProgressBar::new(progress).show_percentage());

                    if total > 0 {
                        ui.label(format!("{} / {} sequences processed", current, total));
                    }
                });
            });
        });

        // Advanced settings window
        if self.show_advanced_settings {
            egui::Window::new("Advanced Settings")
                .collapsible(false)
                .resizable(false)
                .show(ctx, |ui| {
                    ui.heading("Alignment Parameters");
                    ui.add_space(10.0);

                    egui::Grid::new("settings_grid")
                        .num_columns(2)
                        .show(ui, |ui| {
                            ui.label("Alignment Mode:");
                            egui::ComboBox::from_id_salt("mode")
                                .selected_text(match self.settings.mode {
                                    AlignmentMode::Local => "Local",
                                    AlignmentMode::Global => "Global",
                                })
                                .show_ui(ui, |ui| {
                                    ui.selectable_value(
                                        &mut self.settings.mode,
                                        AlignmentMode::Local,
                                        "Local",
                                    );
                                    ui.selectable_value(
                                        &mut self.settings.mode,
                                        AlignmentMode::Global,
                                        "Global",
                                    );
                                });
                            ui.end_row();

                            ui.label("Match Score:");
                            ui.add(
                                egui::DragValue::new(&mut self.settings.match_score)
                                    .range(-10..=10)
                                    .speed(0.1),
                            );
                            ui.end_row();

                            ui.label("Mismatch Score:");
                            ui.add(
                                egui::DragValue::new(&mut self.settings.mismatch_score)
                                    .range(-10..=10)
                                    .speed(0.1),
                            );
                            ui.end_row();

                            ui.label("Gap Open Score:");
                            ui.add(
                                egui::DragValue::new(&mut self.settings.gap_open_score)
                                    .range(-20..=0)
                                    .speed(0.1),
                            );
                            ui.end_row();

                            ui.label("Gap Extend Score:");
                            ui.add(
                                egui::DragValue::new(&mut self.settings.gap_extend_score)
                                    .range(-10..=0)
                                    .speed(0.1),
                            );
                            ui.end_row();

                            ui.label("Min Coverage:");
                            ui.add(
                                egui::DragValue::new(&mut self.settings.min_coverage)
                                    .range(0.0..=1.0)
                                    .speed(0.01)
                                    .max_decimals(2),
                            );
                            ui.end_row();

                            ui.label("Max Mismatches per Oligo:");
                            ui.add(
                                egui::DragValue::new(&mut self.settings.max_mismatches_per_oligo)
                                    .range(0..=50)
                                    .speed(0.1),
                            );
                            ui.end_row();

                            ui.label("Ambiguity Match Display:");
                            egui::ComboBox::from_id_salt("ambiguity_display_adv")
                                .selected_text(match self.settings.ambiguity_display {
                                    AmbiguityDisplayMode::ShowBases => "Show Base Letters",
                                    AmbiguityDisplayMode::ShowDots => "Show Dots",
                                })
                                .show_ui(ui, |ui| {
                                    ui.selectable_value(
                                        &mut self.settings.ambiguity_display,
                                        AmbiguityDisplayMode::ShowDots,
                                        "Show Dots",
                                    );
                                    ui.selectable_value(
                                        &mut self.settings.ambiguity_display,
                                        AmbiguityDisplayMode::ShowBases,
                                        "Show Base Letters",
                                    );
                                });
                            ui.end_row();
                        });

                    ui.add_space(10.0);
                    ui.separator();

                    ui.heading("Per-Category Match Thresholds");
                    ui.add_space(5.0);

                    egui::Grid::new("threshold_grid")
                        .num_columns(2)
                        .show(ui, |ui| {
                            ui.label("Min Forward Primers Matched:");
                            ui.add(
                                egui::DragValue::new(&mut self.settings.min_fwd_matched)
                                    .range(0..=100)
                                    .speed(0.1),
                            );
                            ui.end_row();

                            ui.label("Min Reverse Primers Matched:");
                            ui.add(
                                egui::DragValue::new(&mut self.settings.min_rev_matched)
                                    .range(0..=100)
                                    .speed(0.1),
                            );
                            ui.end_row();

                            ui.label("Min Probes Matched:");
                            ui.add(
                                egui::DragValue::new(&mut self.settings.min_probe_matched)
                                    .range(0..=100)
                                    .speed(0.1),
                            );
                            ui.end_row();
                        });

                    ui.add_space(10.0);
                    ui.separator();

                    ui.heading("Amplicon Size Constraints");
                    ui.add_space(5.0);

                    ui.checkbox(
                        &mut self.amplicon_enabled,
                        "Enable amplicon size limits",
                    );

                    ui.add_enabled_ui(self.amplicon_enabled, |ui| {
                        ui.horizontal(|ui| {
                            ui.label("Min amplicon size (bp):");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.amplicon_min_size_input)
                                    .desired_width(100.0),
                            );

                            ui.add_space(10.0);

                            ui.label("Max amplicon size (bp):");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.amplicon_max_size_input)
                                    .desired_width(100.0),
                            );
                        });
                    });

                    ui.add_space(10.0);
                    ui.separator();

                    ui.horizontal(|ui| {
                        if ui.button("Reset to Defaults").clicked() {
                            self.reset_defaults();
                        }
                        if ui.button("OK").clicked() {
                            self.show_advanced_settings = false;
                        }
                    });
                });
        }

        // Results window
        if self.show_results_window {
            egui::Window::new("Analysis Results")
                .default_size([800.0, 600.0])
                .show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        if ui.button("Save Text").clicked() {
                            if let Some(ref results) = self.results {
                                if let Some(path) = FileDialog::new()
                                    .add_filter("Text files", &["txt"])
                                    .set_title("Save Results")
                                    .save_file()
                                {
                                    if let Err(e) = std::fs::write(&path, &results.output_text) {
                                        self.status_message = format!("Error saving: {}", e);
                                    } else {
                                        self.status_message = "Results saved".to_string();
                                    }
                                }
                            }
                        }
                        if ui.button("Export Excel").clicked() {
                            self.export_to_excel();
                        }
                        if ui.button("Close").clicked() {
                            self.show_results_window = false;
                        }
                    });

                    ui.separator();

                    egui::ScrollArea::both().show(ui, |ui| {
                        if let Some(ref results) = self.results {
                            ui.add(
                                egui::TextEdit::multiline(&mut results.output_text.as_str())
                                    .font(egui::TextStyle::Monospace)
                                    .desired_width(f32::INFINITY),
                            );
                        }
                    });
                });
        }

        // FASTA info window
        if self.show_fasta_info {
            egui::Window::new("FASTA Information")
                .default_size([400.0, 300.0])
                .show(ctx, |ui| {
                    if ui.button("Close").clicked() {
                        self.show_fasta_info = false;
                    }
                    ui.separator();
                    egui::ScrollArea::both().show(ui, |ui| {
                        ui.add(
                            egui::TextEdit::multiline(&mut self.fasta_info_text.as_str())
                                .font(egui::TextStyle::Monospace)
                                .desired_width(f32::INFINITY),
                        );
                    });
                });
        }

        // About window
        if self.show_about {
            egui::Window::new("About")
                .collapsible(false)
                .resizable(false)
                .show(ctx, |ui| {
                    ui.heading("qPCR Oligo Inclusivity Tool");
                    ui.label("Version 2.0.0");
                    ui.add_space(10.0);
                    ui.label("A high-performance oligo inclusivity analysis tool");
                    ui.label("built with Rust and rust-bio.");
                    ui.add_space(10.0);
                    ui.label("Features:");
                    ui.label("  Categorized oligo input (Forward/Reverse/Probe)");
                    ui.label("  Amplicon-aware analysis with fwd+rev pairing");
                    ui.label("  Probe matching within amplicon regions");
                    ui.label("  Bidirectional primer alignment (sense/antisense)");
                    ui.label("  Per-category match thresholds");
                    ui.label("  Parallel processing with Rayon");
                    ui.label("  Excel export with categorized columns");
                    ui.label("  JSON save/load for oligo sets");
                    ui.add_space(10.0);
                    if ui.button("OK").clicked() {
                        self.show_about = false;
                    }
                });
        }
    }
}

// ============================================================================
// Main Entry Point
// ============================================================================

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([950.0, 800.0])
            .with_min_inner_size([700.0, 500.0]),
        ..Default::default()
    };

    eframe::run_native(
        "qPCR Oligo Inclusivity Tool v2.0",
        options,
        Box::new(|cc| Ok(Box::new(PrimerAlignApp::new(cc)))),
    )
}
