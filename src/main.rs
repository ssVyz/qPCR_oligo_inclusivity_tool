//! qPCR Primer Inclusivity Analysis Tool
//! A fast Rust implementation with eframe GUI
//!
//! This tool evaluates primer inclusivity in large sets of reference sequences
//! using pairwise alignment from rust-bio.

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
    matched_oligos: usize,
    examples: Vec<String>,
    amplicon_lengths: Vec<usize>, // Track amplicon lengths for this pattern
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
    min_oligos_matched: usize,
    min_coverage: f64,
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
            min_oligos_matched: 1,
            min_coverage: 0.8,
            ambiguity_display: AmbiguityDisplayMode::ShowBases,
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

/// Analysis results
#[derive(Clone, Debug, Default)]
struct AnalysisResults {
    alignment_dict: HashMap<String, PatternData>,
    oligo_stats: HashMap<String, OligoStats>,
    total_sequences: usize,
    sequences_with_min_matches: usize,
    sequences_with_valid_amplicon: usize,
    sequences_failed_amplicon: usize,
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
            b'R' => &[b'A', b'G'], // purine
            b'Y' => &[b'C', b'T'], // pyrimidine
            b'S' => &[b'G', b'C'],
            b'W' => &[b'A', b'T'],
            b'K' => &[b'G', b'T'],
            b'M' => &[b'A', b'C'],
            b'B' => &[b'C', b'G', b'T'], // not A
            b'D' => &[b'A', b'G', b'T'], // not C
            b'H' => &[b'A', b'C', b'T'], // not G
            b'V' => &[b'A', b'C', b'G'], // not T
            _ => &[], // unknown character, no match
        }
    };
    
    let a_bases = get_bases(a_upper);
    let b_bases = get_bases(b_upper);
    
    // Check if there's any overlap between the two sets of possible bases
    for &base_a in a_bases {
        for &base_b in b_bases {
            // Normalize U to T for comparison
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
) -> (Option<bio::alignment::Alignment>, Option<Orientation>, String) {
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

    // Try sense orientation
    let alignment_sense = match settings.mode {
        AlignmentMode::Local => aligner.local(target_bytes, oligo_bytes),
        AlignmentMode::Global => aligner.global(target_bytes, oligo_bytes),
    };

    // Try antisense orientation
    let alignment_antisense = match settings.mode {
        AlignmentMode::Local => aligner.local(target_bytes, oligo_rc_bytes),
        AlignmentMode::Global => aligner.global(target_bytes, oligo_rc_bytes),
    };

    // Return the better scoring alignment
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
    
    // Normalize U to T for comparison
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
                    
                    // Check if it's an exact match or ambiguity match
                    if is_exact_match(target_base, oligo_base) {
                        sig[q_pos] = '.';
                    } else {
                        // It's an ambiguity match
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
                    
                    // Check if this is actually an IUPAC ambiguity match
                    // (rust-bio marks it as Subst because characters differ,
                    // but our scoring function treated it as a match)
                    if iupac_match(target_base, oligo_base) {
                        // It's an ambiguity match
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
                        // True mismatch
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
    
    // For reverse orientation, reverse complement mismatch letters and reverse the signature
    // to match original oligo orientation
    if orientation == Orientation::Antisense {
        // Reverse complement any nucleotide letters in the signature (mismatches)
        signature = signature
            .chars()
            .map(|c| {
                match c.to_ascii_uppercase() {
                    'A' => 'T',
                    'T' => 'A',
                    'G' => 'C',
                    'C' => 'G',
                    'U' => 'A',
                    'R' => 'Y', // A|G -> C|T
                    'Y' => 'R', // C|T -> A|G
                    'S' => 'S', // G|C -> C|G (same set)
                    'W' => 'W', // A|T -> T|A (same set)
                    'K' => 'M', // G|T -> A|C
                    'M' => 'K', // A|C -> G|T
                    'B' => 'V', // C|G|T -> A|C|G
                    'D' => 'H', // A|G|T -> A|C|T
                    'H' => 'D', // A|C|T -> A|G|T
                    'V' => 'B', // A|C|G -> C|G|T
                    _ => c, // Keep dots, dashes, and other characters as-is
                }
            })
            .collect();
        // Reverse the entire signature string
        signature = signature.chars().rev().collect();
    }

    (signature, mismatches)
}

/// Validate amplicon constraints and filter oligo results accordingly.
/// Returns (updated results, amplicon info, whether a valid amplicon was found)
fn validate_amplicon_constraints(
    oligo_results: &mut HashMap<String, OligoResult>,
    min_size: Option<usize>,
    max_size: Option<usize>,
) -> AmpliconInfo {
    // Collect all forward (sense) hits with positions
    let forward_hits: Vec<(String, usize, usize)> = oligo_results
        .iter()
        .filter(|(_, r)| r.matched && r.orientation == Some(Orientation::Sense))
        .filter_map(|(id, r)| {
            match (r.start_pos, r.end_pos) {
                (Some(s), Some(e)) => Some((id.clone(), s, e)),
                _ => None,
            }
        })
        .collect();

    // Collect all reverse (antisense) hits with positions
    let reverse_hits: Vec<(String, usize, usize)> = oligo_results
        .iter()
        .filter(|(_, r)| r.matched && r.orientation == Some(Orientation::Antisense))
        .filter_map(|(id, r)| {
            match (r.start_pos, r.end_pos) {
                (Some(s), Some(e)) => Some((id.clone(), s, e)),
                _ => None,
            }
        })
        .collect();

    // Find all valid convergent pairs
    // A convergent pair has: forward start < reverse end (they face each other)
    // Amplicon size = reverse end - forward start + 1
    let mut valid_amplicons: Vec<(String, String, usize, usize, usize)> = Vec::new();

    for (fwd_id, fwd_start, _fwd_end) in &forward_hits {
        for (rev_id, _rev_start, rev_end) in &reverse_hits {
            // Check if convergent (forward before reverse on reference)
            if fwd_start < rev_end {
                let amplicon_size = rev_end - fwd_start + 1;
                // Check both min and max constraints if specified
                let size_valid = match (min_size, max_size) {
                    (Some(min), Some(max)) => amplicon_size >= min && amplicon_size <= max,
                    (Some(min), None) => amplicon_size >= min,
                    (None, Some(max)) => amplicon_size <= max,
                    (None, None) => true,
                };
                if size_valid {
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

    // If no valid amplicon found, mark all as no match
    if valid_amplicons.is_empty() {
        for result in oligo_results.values_mut() {
            *result = OligoResult::no_match();
        }
        return AmpliconInfo::default();
    }

    // Find the largest valid amplicon
    let best_amplicon = valid_amplicons
        .iter()
        .max_by_key(|(_, _, _, _, size)| size)
        .unwrap();

    let (best_fwd_id, best_rev_id, amp_start, amp_end, amp_size) = best_amplicon.clone();

    // Check which oligos are inside the amplicon and update results
    for (oligo_id, result) in oligo_results.iter_mut() {
        if !result.matched {
            continue;
        }

        // The convergent pair primers are always accepted
        if oligo_id == &best_fwd_id || oligo_id == &best_rev_id {
            continue;
        }

        // Check if this oligo is inside the amplicon
        if let (Some(start), Some(end)) = (result.start_pos, result.end_pos) {
            // Oligo must be entirely within the amplicon boundaries
            if start < amp_start || end > amp_end {
                // Outside amplicon - mark as no match
                *result = OligoResult::no_match();
            }
            // If inside, keep the match as-is
        } else {
            // No position info - shouldn't happen for matched oligos, but mark as no match
            *result = OligoResult::no_match();
        }
    }

    AmpliconInfo {
        found: true,
        forward_oligo_id: Some(best_fwd_id),
        reverse_oligo_id: Some(best_rev_id),
        start: amp_start,
        end: amp_end,
        size: amp_size,
    }
}

/// Analyze a sequence against all oligos
fn analyze_sequence(
    sequence: &FastaRecord,
    oligos: &[Oligo],
    settings: &AlignmentSettings,
) -> (HashMap<String, OligoResult>, usize, usize, Option<AmpliconInfo>) {
    let mut oligo_results = HashMap::new();

    // First pass: align all oligos to the sequence
    for oligo in oligos {
        let (alignment_opt, orientation_opt, actual_oligo) =
            get_best_alignment(&sequence.seq, &oligo.seq, settings);

        let result = if let (Some(alignment), Some(orientation)) = (alignment_opt, orientation_opt)
        {
            if is_valid_alignment(&alignment, actual_oligo.len(), settings.min_coverage) {
                let (signature, mismatches) =
                    generate_signature(&alignment, &sequence.seq, &actual_oligo, orientation, settings.ambiguity_display);
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
            } else {
                OligoResult::no_match()
            }
        } else {
            OligoResult::no_match()
        };

        oligo_results.insert(oligo.id.clone(), result);
    }

    // Second pass: apply amplicon constraints if enabled
    let amplicon_info = if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
        Some(validate_amplicon_constraints(&mut oligo_results, settings.min_amplicon_size, settings.max_amplicon_size))
    } else {
        None
    };

    // Calculate final stats after amplicon filtering
    let mut total_matched = 0;
    let mut total_mismatches = 0;
    
    for result in oligo_results.values() {
        if result.matched {
            total_matched += 1;
            total_mismatches += result.mismatches;
        }
    }

    (oligo_results, total_matched, total_mismatches, amplicon_info)
}

/// Run the full analysis with parallel processing
fn run_analysis(
    sequences: &[FastaRecord],
    oligos: &[Oligo],
    settings: &AlignmentSettings,
    progress: &ProgressTracker,
) -> AnalysisResults {
    let total_sequences = sequences.len();
    progress.total.store(total_sequences, Ordering::SeqCst);
    progress.current.store(0, Ordering::SeqCst);

    let oligo_stats: HashMap<String, OligoStats> = oligos
        .iter()
        .map(|o| (o.id.clone(), OligoStats::default()))
        .collect();
    let oligo_stats = Arc::new(Mutex::new(oligo_stats));

    let alignment_dict: Arc<Mutex<HashMap<String, PatternData>>> =
        Arc::new(Mutex::new(HashMap::new()));
    let sequences_with_min_matches = Arc::new(AtomicUsize::new(0));
    let sequences_with_valid_amplicon = Arc::new(AtomicUsize::new(0));
    let sequences_failed_amplicon = Arc::new(AtomicUsize::new(0));

    let processed = Arc::new(AtomicUsize::new(0));

    // Process sequences in parallel using rayon
    sequences.par_iter().for_each(|record| {
        if !progress.running.load(Ordering::SeqCst) {
            return;
        }

        let (oligo_results, matched_count, total_mismatches, amplicon_info) =
            analyze_sequence(record, oligos, settings);

        // Track amplicon statistics
        if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
            if let Some(ref info) = amplicon_info {
                if info.found {
                    sequences_with_valid_amplicon.fetch_add(1, Ordering::SeqCst);
                } else {
                    sequences_failed_amplicon.fetch_add(1, Ordering::SeqCst);
                }
            } else {
                // No amplicon info means no valid amplicon was found
                sequences_failed_amplicon.fetch_add(1, Ordering::SeqCst);
            }
        }

        // Update oligo stats
        {
            let mut stats = oligo_stats.lock().unwrap();
            for (oligo_id, result) in &oligo_results {
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

        if matched_count >= settings.min_oligos_matched {
            sequences_with_min_matches.fetch_add(1, Ordering::SeqCst);

            // Build combined signature
            let mut signature_parts = Vec::new();
            for oligo in oligos {
                if let Some(result) = oligo_results.get(&oligo.id) {
                    if result.matched {
                        let orientation_symbol = match result.orientation {
                            Some(Orientation::Sense) => "(fwd)",
                            Some(Orientation::Antisense) => "(rev)",
                            None => "",
                        };
                        signature_parts.push(format!("{}{}", result.signature, orientation_symbol));
                    } else {
                        signature_parts.push("NO_MATCH".to_string());
                    }
                } else {
                    signature_parts.push("NO_MATCH".to_string());
                }
            }

            let combined_signature = signature_parts.join(" | ");

            // Update alignment dict
            {
                let mut dict = alignment_dict.lock().unwrap();
                let entry = dict.entry(combined_signature).or_insert_with(|| PatternData {
                    count: 0,
                    total_mismatches,
                    matched_oligos: matched_count,
                    examples: Vec::new(),
                    amplicon_lengths: Vec::new(),
                });
                entry.count += 1;
                if entry.examples.len() < 10 {
                    entry.examples.push(record.id.clone());
                }
                // Track amplicon length if available
                if let Some(ref info) = amplicon_info {
                    if info.found {
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
    let final_sequences_with_valid_amplicon = sequences_with_valid_amplicon.load(Ordering::SeqCst);
    let final_sequences_failed_amplicon = sequences_failed_amplicon.load(Ordering::SeqCst);

    let output_text = generate_output_text(
        oligos,
        &final_oligo_stats,
        &final_alignment_dict,
        total_sequences,
        final_sequences_with_min_matches,
        final_sequences_with_valid_amplicon,
        final_sequences_failed_amplicon,
        settings,
    );

    AnalysisResults {
        alignment_dict: final_alignment_dict,
        oligo_stats: final_oligo_stats,
        total_sequences,
        sequences_with_min_matches: final_sequences_with_min_matches,
        sequences_with_valid_amplicon: final_sequences_with_valid_amplicon,
        sequences_failed_amplicon: final_sequences_failed_amplicon,
        output_text,
    }
}

/// Generate formatted output text
fn generate_output_text(
    oligos: &[Oligo],
    oligo_stats: &HashMap<String, OligoStats>,
    alignment_dict: &HashMap<String, PatternData>,
    total_sequences: usize,
    sequences_with_min_matches: usize,
    sequences_with_valid_amplicon: usize,
    sequences_failed_amplicon: usize,
    settings: &AlignmentSettings,
) -> String {
    let mut out = Vec::new();

    out.push("=".repeat(80));
    out.push("ENHANCED qPCR ALIGNMENT ANALYSIS RESULTS".to_string());
    out.push("=".repeat(80));
    out.push(String::new());

    let header_parts: Vec<String> = oligos
        .iter()
        .map(|o| format!("{}({})", o.id, o.seq))
        .collect();
    out.push(format!("Oligos analyzed: {}", header_parts.join(", ")));
    out.push(format!(
        "Analysis settings: Min matches = {}, Min coverage = {}",
        settings.min_oligos_matched, settings.min_coverage
    ));
    if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
        match (settings.min_amplicon_size, settings.max_amplicon_size) {
            (Some(min), Some(max)) => {
                out.push(format!("Amplicon constraint: Min size = {} bp, Max size = {} bp", min, max));
            }
            (Some(min), None) => {
                out.push(format!("Amplicon constraint: Min size = {} bp", min));
            }
            (None, Some(max)) => {
                out.push(format!("Amplicon constraint: Max size = {} bp", max));
            }
            (None, None) => {}
        }
    }
    out.push(String::new());

    out.push("SEQUENCE SIGNATURE PATTERNS:".to_string());
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
                "  Count: {}, Mismatches: {}, Matched oligos: {}",
                data.count, data.total_mismatches, data.matched_oligos
            ));
            out.push(format!("  Examples: {}", examples_str.join(", ")));
            out.push(String::new());
        }
    } else {
        out.push("No sequences met the minimum matching criteria.".to_string());
        out.push(String::new());
    }

    out.push("PER-OLIGO STATISTICS:".to_string());
    out.push("-".repeat(30));
    for oligo in oligos {
        if let Some(stats) = oligo_stats.get(&oligo.id) {
            let percentage = if total_sequences > 0 {
                (stats.total_matches as f64 / total_sequences as f64) * 100.0
            } else {
                0.0
            };
            out.push(format!(
                "{}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
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
        "Sequences with >={} matches: {} ({:.1}%)",
        settings.min_oligos_matched, sequences_with_min_matches, percentage
    ));
    out.push(format!(
        "Sequences with <{} matches: {}",
        settings.min_oligos_matched,
        total_sequences - sequences_with_min_matches
    ));

    // Add amplicon statistics if enabled
    if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
        out.push(String::new());
        out.push("AMPLICON STATISTICS:".to_string());
        out.push("-".repeat(20));
        let amp_percentage = if total_sequences > 0 {
            (sequences_with_valid_amplicon as f64 / total_sequences as f64) * 100.0
        } else {
            0.0
        };
        out.push(format!(
            "Sequences with valid amplicon: {} ({:.1}%)",
            sequences_with_valid_amplicon, amp_percentage
        ));
        out.push(format!(
            "Sequences failed amplicon check: {}",
            sequences_failed_amplicon
        ));
        let constraint_desc = match (settings.min_amplicon_size, settings.max_amplicon_size) {
            (Some(min), Some(max)) => format!("between {} and {} bp", min, max),
            (Some(min), None) => format!(">= {} bp", min),
            (None, Some(max)) => format!("<= {} bp", max),
            (None, None) => "".to_string(),
        };
        out.push(format!(
            "(Note: Amplicon check requires convergent fwd+rev pair {})",
            constraint_desc
        ));
    }

    out.push(String::new());
    out.push("LEGEND:".to_string());
    out.push("'(fwd)' = sense orientation, '(rev)' = antisense orientation".to_string());
    out.push("'.' = match, letter = mismatch, '-' = gap/unaligned".to_string());
    if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
        out.push("Amplicon = largest convergent fwd-rev pair within size constraints".to_string());
        out.push("Oligos outside the amplicon boundaries are marked as NO_MATCH".to_string());
    }
    out.push("=".repeat(80));

    out.join("\n")
}

// ============================================================================
// GUI Application
// ============================================================================

/// Main application state
struct PrimerAlignApp {
    sequence_file: Option<PathBuf>,
    oligo_file: Option<PathBuf>,
    output_file: Option<PathBuf>,
    settings: AlignmentSettings,
    sequences: Vec<FastaRecord>,
    oligos: Vec<Oligo>,
    oligo_text_input: String,
    results: Option<AnalysisResults>,
    progress: ProgressTracker,
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
            oligo_file: None,
            output_file: None,
            settings: AlignmentSettings::default(),
            sequences: Vec::new(),
            oligos: Vec::new(),
            oligo_text_input: String::new(),
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

    fn select_oligo_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("FASTA files", &["fasta", "fas", "fa", "fna"])
            .add_filter("All files", &["*"])
            .set_title("Select Oligo File")
            .pick_file()
        {
            match parse_fasta(&path) {
                Ok(records) => {
                    self.oligos = records
                        .into_iter()
                        .map(|r| Oligo { id: r.id, seq: r.seq })
                        .collect();
                    self.oligo_file = Some(path);
                    self.oligo_text_input = oligos_to_fasta_string(&self.oligos);
                    self.status_message = format!("Loaded {} oligos", self.oligos.len());
                }
                Err(e) => {
                    self.status_message = format!("Error loading oligos: {}", e);
                }
            }
        }
    }

    fn parse_oligos_from_text(&mut self) {
        if self.oligo_text_input.trim().is_empty() {
            self.oligos.clear();
            self.status_message = "No oligo sequences entered".to_string();
            return;
        }

        match parse_fasta_string(&self.oligo_text_input) {
            Ok(oligos) => {
                self.oligos = oligos;
                self.status_message = format!("Parsed {} oligos from text input", self.oligos.len());
            }
            Err(e) => {
                self.status_message = format!("Error parsing oligos: {}", e);
            }
        }
    }

    fn select_output_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("Text files", &["txt"])
            .add_filter("All files", &["*"])
            .set_title("Save Analysis Results As")
            .save_file()
        {
            self.output_file = Some(path);
            self.status_message = "Output file selected".to_string();
        }
    }

    fn update_amplicon_setting(&mut self) {
        if self.amplicon_enabled {
            // Parse minimum amplicon size
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
            
            // Parse maximum amplicon size
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

    fn start_analysis(&mut self) {
        if self.sequences.is_empty() {
            self.status_message = "Please select a sequence file first".to_string();
            return;
        }

        // Parse oligos from text input before analysis
        self.parse_oligos_from_text();

        if self.oligos.is_empty() {
            self.status_message = "Please enter oligo sequences or import an oligo file".to_string();
            return;
        }

        // Update amplicon setting before starting
        self.update_amplicon_setting();

        let sequences = self.sequences.clone();
        let oligos = self.oligos.clone();
        let settings = self.settings.clone();
        let progress = self.progress.clone();

        self.progress.running.store(true, Ordering::SeqCst);
        self.progress.current.store(0, Ordering::SeqCst);
        self.progress.total.store(sequences.len(), Ordering::SeqCst);
        if let Ok(mut status) = self.progress.status.lock() {
            *status = "Starting analysis...".to_string();
        }

        let handle = thread::spawn(move || run_analysis(&sequences, &oligos, &settings, &progress));

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

                        if let Some(ref path) = self.output_file {
                            if let Some(ref results) = self.results {
                                if let Err(e) = std::fs::write(path, &results.output_text) {
                                    self.status_message = format!("Error saving results: {}", e);
                                } else {
                                    self.status_message = format!(
                                        "Results saved to {}",
                                        path.file_name()
                                            .map(|s| s.to_string_lossy().to_string())
                                            .unwrap_or_default()
                                    );
                                }
                            }
                        }
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
            let avg_len: f64 = seq_lengths.iter().sum::<usize>() as f64 / seq_lengths.len() as f64;

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

        let mut row: u32 = 0;

        // Title
        worksheet
            .write_string_with_format(row, 0, "qPCR Alignment Analysis Results", &title_format)
            .map_err(|e| e.to_string())?;
        row += 2;

        // Settings
        let oligo_info: Vec<String> = self
            .oligos
            .iter()
            .map(|o| format!("{}({})", o.id, o.seq))
            .collect();
        worksheet
            .write_string(row, 0, &format!("Oligos analyzed: {}", oligo_info.join(", ")))
            .map_err(|e| e.to_string())?;
        row += 1;
        worksheet
            .write_string(
                row,
                0,
                &format!(
                    "Settings: Min matches = {}, Min coverage = {}",
                    self.settings.min_oligos_matched, self.settings.min_coverage
                ),
            )
            .map_err(|e| e.to_string())?;
        row += 1;
        
        if self.settings.min_amplicon_size.is_some() || self.settings.max_amplicon_size.is_some() {
            let constraint_text = match (self.settings.min_amplicon_size, self.settings.max_amplicon_size) {
                (Some(min), Some(max)) => format!("Amplicon constraint: Min size = {} bp, Max size = {} bp", min, max),
                (Some(min), None) => format!("Amplicon constraint: Min size = {} bp", min),
                (None, Some(max)) => format!("Amplicon constraint: Max size = {} bp", max),
                (None, None) => "".to_string(),
            };
            worksheet
                .write_string(row, 0, &constraint_text)
                .map_err(|e| e.to_string())?;
            row += 1;
        }
        row += 2;

        // Headers
        let mut col: u16 = 0;
        worksheet
            .write_string_with_format(row, col, "Pattern #", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;

        for oligo in &self.oligos {
            worksheet
                .write_string_with_format(row, col, &format!("{} Pattern", oligo.id), &header_format)
                .map_err(|e| e.to_string())?;
            col += 1;
        }

        worksheet
            .write_string_with_format(row, col, "Count", &header_format)
            .map_err(|e| e.to_string())?;
        col += 1;
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
            .write_string_with_format(row, col, "Matched Oligos", &header_format)
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
        col = 0;
        worksheet
            .write_string(row, col, "")
            .map_err(|e| e.to_string())?;
        col += 1;
        for oligo in &self.oligos {
            worksheet
                .write_string(row, col, &oligo.seq)
                .map_err(|e| e.to_string())?;
            col += 1;
        }
        // Empty cells for the remaining columns
        for _ in 0..5 {
            worksheet
                .write_string(row, col, "")
                .map_err(|e| e.to_string())?;
            col += 1;
        }

        row += 1;

        // Data rows
        let mut sorted_patterns: Vec<_> = results.alignment_dict.iter().collect();
        sorted_patterns.sort_by(|a, b| b.1.count.cmp(&a.1.count));

        let mut pattern_num = 1u32;
        for (signature, data) in sorted_patterns {
            col = 0;
            worksheet
                .write_number(row, col, pattern_num as f64)
                .map_err(|e| e.to_string())?;
            col += 1;

            let signature_parts: Vec<&str> = signature.split(" | ").collect();
            for (i, _) in self.oligos.iter().enumerate() {
                let pattern = signature_parts.get(i).unwrap_or(&"NO_MATCH");
                worksheet
                    .write_string(row, col, pattern.to_string())
                    .map_err(|e| e.to_string())?;
                col += 1;
            }

            worksheet
                .write_number(row, col, data.count as f64)
                .map_err(|e| e.to_string())?;
            col += 1;
            // Calculate percentage
            let percentage = if results.total_sequences > 0 {
                (data.count as f64 / results.total_sequences as f64) * 100.0
            } else {
                0.0
            };
            worksheet
                .write_number(row, col, percentage)
                .map_err(|e| e.to_string())?;
            col += 1;
            worksheet
                .write_number(row, col, data.total_mismatches as f64)
                .map_err(|e| e.to_string())?;
            col += 1;
            worksheet
                .write_number(row, col, data.matched_oligos as f64)
                .map_err(|e| e.to_string())?;
            col += 1;
            // Calculate average amplicon length (or most common)
            let amplicon_length = if !data.amplicon_lengths.is_empty() {
                // Use the most common amplicon length, or average if tied
                use std::collections::HashMap;
                let mut length_counts: HashMap<usize, usize> = HashMap::new();
                for &len in &data.amplicon_lengths {
                    *length_counts.entry(len).or_insert(0) += 1;
                }
                // Find the most common length
                length_counts
                    .iter()
                    .max_by_key(|(_, &count)| count)
                    .map(|(&len, _)| len)
                    .unwrap_or(0)
            } else {
                0
            };
            worksheet
                .write_number(row, col, amplicon_length as f64)
                .map_err(|e| e.to_string())?;
            col += 1;

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

        // Statistics
        row += 2;
        worksheet
            .write_string_with_format(row, 0, "PER-OLIGO STATISTICS:", &header_format)
            .map_err(|e| e.to_string())?;
        row += 1;

        for oligo in &self.oligos {
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
                            "{}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
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
                    "Sequences with >={} matches: {} ({:.1}%)",
                    self.settings.min_oligos_matched, results.sequences_with_min_matches, percentage
                ),
            )
            .map_err(|e| e.to_string())?;
        row += 1;

        // Amplicon statistics if enabled
        if self.settings.min_amplicon_size.is_some() || self.settings.max_amplicon_size.is_some() {
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
                        "Sequences failed amplicon check: {}",
                        results.sequences_failed_amplicon
                    ),
                )
                .map_err(|e| e.to_string())?;
        }

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
                    if ui.button("Select Oligo File...").clicked() {
                        self.select_oligo_file();
                        ui.close_menu();
                    }
                    if ui.button("Set Output File...").clicked() {
                        self.select_output_file();
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
            ui.heading("qPCR Primer Alignment Tool");
            ui.label("High-performance primer inclusivity analysis using Rust");
            ui.add_space(10.0);

            // File selection group
            ui.group(|ui| {
                ui.label(egui::RichText::new("Input Files").strong());
                ui.add_space(5.0);

                egui::Grid::new("file_grid").num_columns(3).show(ui, |ui| {
                    ui.label("Sequence File:");
                    if ui.button("Browse...").clicked() {
                        self.select_sequence_file();
                    }
                    if let Some(ref path) = self.sequence_file {
                        ui.colored_label(
                            egui::Color32::DARK_GREEN,
                            format!(
                                "✓ {} ({} sequences)",
                                path.file_name()
                                    .map(|s| s.to_string_lossy().to_string())
                                    .unwrap_or_default(),
                                self.sequences.len()
                            ),
                        );
                    } else {
                        ui.label("No file selected");
                    }
                    ui.end_row();

                    ui.label("Oligo File:");
                    if ui.button("Import...").clicked() {
                        self.select_oligo_file();
                    }
                    if let Some(ref path) = self.oligo_file {
                        ui.colored_label(
                            egui::Color32::DARK_GREEN,
                            format!(
                                "✓ {}",
                                path.file_name()
                                    .map(|s| s.to_string_lossy().to_string())
                                    .unwrap_or_default()
                            ),
                        );
                    } else {
                        ui.label("(optional - use text input below)");
                    }
                    ui.end_row();

                    ui.label("Output File:");
                    if ui.button("Browse...").clicked() {
                        self.select_output_file();
                    }
                    if let Some(ref path) = self.output_file {
                        ui.colored_label(
                            egui::Color32::DARK_GREEN,
                            format!(
                                "✓ {}",
                                path.file_name()
                                    .map(|s| s.to_string_lossy().to_string())
                                    .unwrap_or_default()
                            ),
                        );
                    } else {
                        ui.label("No file selected (optional)");
                    }
                    ui.end_row();
                });
            });

            ui.add_space(10.0);

            // Oligo sequences input group
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new("Oligo Sequences (FASTA format)").strong());
                    ui.add_space(10.0);
                    if !self.oligo_text_input.is_empty() {
                        // Try to count oligos in real-time
                        let oligo_count = self.oligo_text_input.lines()
                            .filter(|line| line.starts_with('>'))
                            .count();
                        ui.colored_label(
                            egui::Color32::DARK_GREEN,
                            format!("{} oligo(s) defined", oligo_count),
                        );
                    }
                });
                ui.add_space(5.0);

                egui::ScrollArea::vertical()
                    .max_height(120.0)
                    .show(ui, |ui| {
                        ui.add(
                            egui::TextEdit::multiline(&mut self.oligo_text_input)
                                .font(egui::TextStyle::Monospace)
                                .desired_width(f32::INFINITY)
                                .desired_rows(6)
                                .hint_text(">Primer1\nATGCGTACGTAGC\n>Primer2\nGCTAGCTAGCTA"),
                        );
                    });

                ui.add_space(5.0);
                ui.label(
                    egui::RichText::new("Enter sequences in FASTA format or import from file above")
                        .small()
                        .color(egui::Color32::GRAY)
                );
            });

            ui.add_space(10.0);

            // Quick settings group
            ui.group(|ui| {
                ui.label(egui::RichText::new("Quick Settings").strong());
                ui.add_space(5.0);

                ui.horizontal(|ui| {
                    ui.label("Min Oligos Matched:");
                    ui.add(
                        egui::DragValue::new(&mut self.settings.min_oligos_matched)
                            .range(1..=100)
                            .speed(0.1),
                    );

                    ui.add_space(20.0);

                    ui.label("Min Coverage:");
                    ui.add(
                        egui::DragValue::new(&mut self.settings.min_coverage)
                            .range(0.0..=1.0)
                            .speed(0.01)
                            .max_decimals(2),
                    );

                    ui.add_space(20.0);

                    if ui.button("Advanced Settings...").clicked() {
                        self.show_advanced_settings = true;
                    }
                });
                
                ui.add_space(5.0);
                
                ui.horizontal(|ui| {
                    ui.label("Ambiguity Match Display:");
                    egui::ComboBox::from_id_salt("ambiguity_display")
                        .selected_text(match self.settings.ambiguity_display {
                            AmbiguityDisplayMode::ShowBases => "Show Base Letters",
                            AmbiguityDisplayMode::ShowDots => "Show Dots",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.settings.ambiguity_display,
                                AmbiguityDisplayMode::ShowBases,
                                "Show Base Letters",
                            );
                            ui.selectable_value(
                                &mut self.settings.ambiguity_display,
                                AmbiguityDisplayMode::ShowDots,
                                "Show Dots",
                            );
                        });
                });
            });

            ui.add_space(10.0);

            // Amplicon constraint settings
            ui.group(|ui| {
                ui.label(egui::RichText::new("Amplicon Constraints").strong());
                ui.add_space(5.0);

                ui.horizontal(|ui| {
                    ui.checkbox(&mut self.amplicon_enabled, "Enable amplicon size check");
                });
                
                if self.amplicon_enabled {
                    ui.add_space(5.0);
                    ui.add_enabled_ui(self.amplicon_enabled, |ui| {
                        ui.horizontal(|ui| {
                            ui.label("Min amplicon size (bp):");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.amplicon_min_size_input)
                                    .desired_width(80.0)
                            );
                            
                            ui.add_space(10.0);
                            
                            ui.label("Max amplicon size (bp):");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.amplicon_max_size_input)
                                    .desired_width(80.0)
                            );
                        });
                    });
                    
                    ui.add_space(5.0);
                    ui.label(
                        egui::RichText::new(
                            "ℹ Requires convergent fwd+rev pair. If no amplicon within bounds, sequence = NO_MATCH"
                        )
                        .small()
                        .color(egui::Color32::GRAY)
                    );
                }
            });

            ui.add_space(15.0);

            // Action buttons
            ui.horizontal(|ui| {
                let is_running = self.progress.running.load(Ordering::SeqCst);

                ui.add_enabled_ui(!is_running, |ui| {
                    if ui
                        .button(egui::RichText::new("▶ Run Analysis").strong())
                        .clicked()
                    {
                        self.start_analysis();
                    }
                });

                ui.add_enabled_ui(!is_running, |ui| {
                    if ui.button("💾 File Only").clicked() {
                        if self.output_file.is_none() {
                            self.status_message = "Please select an output file first".to_string();
                        } else {
                            self.start_analysis();
                        }
                    }
                });

                if ui.button("📊 FASTA Info").clicked() {
                    self.show_fasta_info();
                }

                if ui.button("📈 Export Excel").clicked() {
                    self.export_to_excel();
                }

                if is_running {
                    ui.add_space(10.0);
                    if ui.button("⏹ Stop").clicked() {
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

        // Advanced settings window
        if self.show_advanced_settings {
            egui::Window::new("Advanced Settings")
                .collapsible(false)
                .resizable(false)
                .show(ctx, |ui| {
                    ui.heading("Alignment Parameters");
                    ui.add_space(10.0);

                    egui::Grid::new("settings_grid").num_columns(2).show(ui, |ui| {
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

                        ui.label("Ambiguity Match Display:");
                        egui::ComboBox::from_id_salt("ambiguity_display_adv")
                            .selected_text(match self.settings.ambiguity_display {
                                AmbiguityDisplayMode::ShowBases => "Show Base Letters",
                                AmbiguityDisplayMode::ShowDots => "Show Dots",
                            })
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut self.settings.ambiguity_display,
                                    AmbiguityDisplayMode::ShowBases,
                                    "Show Base Letters",
                                );
                                ui.selectable_value(
                                    &mut self.settings.ambiguity_display,
                                    AmbiguityDisplayMode::ShowDots,
                                    "Show Dots",
                                );
                            });
                        ui.end_row();
                    });

                    ui.add_space(10.0);
                    ui.separator();
                    
                    ui.heading("Amplicon Constraints");
                    ui.add_space(5.0);
                    
                    ui.checkbox(&mut self.amplicon_enabled, "Enable amplicon size check");
                    
                    ui.add_enabled_ui(self.amplicon_enabled, |ui| {
                        ui.horizontal(|ui| {
                            ui.label("Min amplicon size (bp):");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.amplicon_min_size_input)
                                    .desired_width(100.0)
                            );
                            
                            ui.add_space(10.0);
                            
                            ui.label("Max amplicon size (bp):");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.amplicon_max_size_input)
                                    .desired_width(100.0)
                            );
                        });
                        
                        ui.add_space(5.0);
                        ui.label("Algorithm:");
                        ui.label("1. Find all convergent fwd+rev pairs");
                        ui.label("2. Filter pairs within size constraints (min/max)");
                        ui.label("3. Select largest valid amplicon");
                        ui.label("4. Keep oligos inside amplicon");
                        ui.label("5. If no valid amplicon found, mark all as NO_MATCH");
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
                        if ui.button("💾 Save Text").clicked() {
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
                        if ui.button("📈 Export Excel").clicked() {
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
                    ui.heading("qPCR Primer Alignment Tool");
                    ui.label("Version 1.0.0");
                    ui.add_space(10.0);
                    ui.label("A high-performance primer inclusivity analysis tool");
                    ui.label("built with Rust and rust-bio.");
                    ui.add_space(10.0);
                    ui.label("Features:");
                    ui.label("• Bidirectional primer alignment (sense/antisense)");
                    ui.label("• Flexible matching criteria");
                    ui.label("• Amplicon size constraints");
                    ui.label("• Parallel processing with Rayon");
                    ui.label("• Excel export support");
                    ui.label("• Progress tracking");
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
            .with_inner_size([900.0, 750.0])
            .with_min_inner_size([600.0, 400.0]),
        ..Default::default()
    };

    eframe::run_native(
        "qPCR Primer Alignment Tool",
        options,
        Box::new(|cc| Ok(Box::new(PrimerAlignApp::new(cc)))),
    )
}
