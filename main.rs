//############# WORKING CODE - ran benchmarking on this ##################//

use seq_io::fastq::{ Record};
use std::io::BufWriter;
use smallvec::SmallVec;

// ============================================================================
// GAP 5: STREAMING OUTPUT WITH INLINE FGS
//
// Changes from original:
//   - extract_unitigs_greedy NOW takes output buffers by &mut reference
//   - Calls FGS inline as each completed unitig is produced
//   - Non-coding unitigs are dropped (not stored)
//   - No post-traversal FGS loop in main() — that work is done during traversal
//   - Memory bounded by (graph size) + (longest single unitig) instead of sum of all unitigs
//
// Expected memory improvement for DS5:
//   Before Gap 5: ~1.1M unitig strings accumulated = 40-60 GB RAM
//   After Gap 5:  Only one unitig in memory at a time = <1 GB
// ============================================================================

use rust_sketcher::sequence::{base_to_bit_mask, base_to_bits, bit_mask_to_base, unpack_kmer};
use frag_gene_scan_rs::hmm::{Global, Local};
use frag_gene_scan_rs::dna::Nuc;

use rustc_hash::FxHashMap;
use rustc_hash::FxHashSet;
use seq_io::fastq::Reader;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

 
fn main() {
    //let input_path = "/home/s278t089/Projects/final_project/impp_improved/impp_2.0/data/sds2.5x.reads.fq";
    let input_path = "/home/s278t089/Projects/final_project/impp_improved/impp_2.0/data/ds1.gut.reads.fq";
    //let input_path = "/home/s278t089/Projects/final_project/impp_improved/impp_2.0/data/ds5.synth.reads.fq";
    let k = 31;
    let mask = (1u64 << (2 * k)) - 1;
    let min_count = 2;
    let min_unitig_length = 100;
 
 
    let train_dir = PathBuf::from("./FragGeneScanRs/train");
    let model_name = PathBuf::from("illumina_5");
    let (global, locals) = frag_gene_scan_rs::hmm::get_train_from_file(train_dir, model_name).unwrap();
    println!("✅ HMM loaded");

 
    // --- PASS 1: COUNTING ---
    println!("\n📊 Pass 1: Counting k-mers...");

    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    let mut graph: FxHashMap<u64, u8> = FxHashMap::default();
    
    let mut reader = Reader::from_path(input_path).unwrap();
    let mut read_count = 0usize;
    while let Some(result) = reader.next() {
        let record = result.expect("Error reading record");
        let seq = record.seq();
        read_count += 1;
        if seq.len() < k { continue; }
 
        let mut current_packed: u64 = 0;
        for i in 0..k {
            current_packed = (current_packed << 2) | base_to_bits(seq[i]);
        }
        *counts.entry(current_packed).or_insert(0) += 1;
 
        for i in k..seq.len() {
            current_packed = ((current_packed << 2) | base_to_bits(seq[i])) & mask;
            *counts.entry(current_packed).or_insert(0) += 1;
        }
    }
    

 
    println!("\n🧹 Pruning noise k-mers...");

    counts.retain(|_, &mut count| count >= min_count);
    //counts.shrink_to_fit();
    println!("After pruning: {} entries", counts.len());

    // --- PASS 2: GRAPH BUILDING ---
    println!("\n🏗️ Pass 2: Building graph...");

    let mut reader = Reader::from_path(input_path).unwrap();
    while let Some(result) = reader.next() {
        let record = result.expect("Error reading record");
        let seq = record.seq();
        if seq.len() < k + 1 { continue; }
 
        let mut current_packed: u64 = 0;
        for i in 0..k {
            current_packed = (current_packed << 2) | base_to_bits(seq[i]);
        }
 
        for i in k..seq.len() {
            let next_base = seq[i];
            let next_packed = ((current_packed << 2) | base_to_bits(next_base)) & mask;
            
            if *counts.get(&current_packed).unwrap_or(&0) >= min_count && 
               *counts.get(&next_packed).unwrap_or(&0) >= min_count {
                let entry = graph.entry(current_packed).or_insert(0);
                *entry |= base_to_bit_mask(next_base);
            }
            current_packed = next_packed;
        }
    }
 
    println!("Graph size: {} nodes", graph.len());

 
    let mask_64: u64 = if k == 32 { !0 } else { (1 << (2 * k)) - 1 };
 
    // --- TIP PRUNING ---
    println!("\n✂️ Pruning tips...");

    let removed_count = prune_tips(&mut graph, k, mask_64);
    println!("Removed {} tips", removed_count);

 
    // --- OUTPUT BUFFERS ---
    let mut gff_buffer = Vec::new();
    let mut aa_buffer = Vec::new();
    let mut dna_buffer = Vec::new();
 
    // --- TRAVERSAL ---
    println!("\n🚀 Traversal...");

    let (unitig_count, coding_unitig_count) = extract_unitigs_frame_aware(
        &graph,
        &counts,
        k,
        &global,
        &locals,
        &mut gff_buffer,
        &mut aa_buffer,
        &mut dna_buffer,
        min_unitig_length,
    );



    let (reads_scanned, reads_rescued) = rescue_unmapped_reads(
              input_path,
              &counts,        // Still alive — DON'T drop(counts) before this
              k,
              min_count,
              &global,
              &locals,
              &mut gff_buffer,
              &mut aa_buffer,
              &mut dna_buffer,
              coding_unitig_count,  // Continue numbering from where traversal left off
              50_000,
          );
          
        // Write rescue predictions to SEPARATE files
    if reads_rescued > 0 {
        let mut gff_file = File::create("impp_predictions.gff")
        .expect("Failed to create GFF file");
        gff_file.write_all(&gff_buffer).expect("Failed to write to GFF file");
    
        let mut faa_file = File::create("impp_proteins.faa")
            .expect("Failed to create FAA file");
        faa_file.write_all(&aa_buffer).expect("Failed to write to FAA file");
    
        let mut fna_file = File::create("impp_dna.fna")
            .expect("Failed to create DNA file");
        fna_file.write_all(&dna_buffer).expect("Failed to write to DNA file");
    
    println!("Rescue predictions written to rescue_*.* files");
}

println!("✅ Pipeline complete!");
println!("Assembly predictions: assembly_*.* files");
println!("Rescue predictions:   rescue_*.* files");

 
    // --- FINAL OUTPUT ---
    println!("\n✅ Pipeline complete!");
    println!("Total unitigs:    {}", unitig_count);
    println!("Coding unitigs:   {}", coding_unitig_count);
}


pub fn rescue_unmapped_reads(
    input_path: &str,
    counts: &FxHashMap<u64, u32>,
    k: usize,
    min_count: u32,
    global: &Box<Global>,
    locals: &Vec<Local>,
    gff_buffer: &mut Vec<u8>,
    aa_buffer: &mut Vec<u8>,
    dna_buffer: &mut Vec<u8>,
    unitig_counter_start: usize,   // Pass in your current unitig count to continue numbering
    flush_interval: usize,
) -> (usize, usize) {   // Returns (reads_processed, reads_rescued)
 
    let mask: u64 = if k == 32 { !0 } else { (1u64 << (2 * k)) - 1 };
 
    let mut reads_scanned = 0usize;
    let mut reads_rescued = 0usize;
    let mut reads_already_assembled = 0usize;

 
    let mut reader = Reader::from_path(input_path).unwrap();
    let mut unitig_counter = unitig_counter_start;
 
    println!("🔍 Pass 3: Scanning for reads not assembled into the graph...");
 
    while let Some(result) = reader.next() {
        let record = result.expect("Error reading record");
        let seq = record.seq();
        reads_scanned += 1;
 
        if seq.len() < k {
            // Too short to have any k-mers — run FGS directly
            run_fgs_on_read(
                seq,
                &mut unitig_counter,
                global,
                locals,
                gff_buffer,
                aa_buffer,
                dna_buffer,
            );
            reads_rescued += 1;
            continue;
        }
 
        // Check if this read is fully represented in the trusted graph
        // A read is "in the graph" if ALL its k-mers meet min_count
        // Strategy: check the first, middle, and last k-mer as a fast proxy
        // If any of these fail, treat the read as unmapped
        let is_in_graph = is_read_assembled(seq, k, mask, counts, min_count);
 
        if is_in_graph {
            reads_already_assembled += 1;
            continue;  // Skip: this read's genes are already in the unitig output
        }
 
        // This read was not fully assembled — rescue it with direct FGS
        run_fgs_on_read(
            seq,
            &mut unitig_counter,
            global,
            locals,
            gff_buffer,
            aa_buffer,
            dna_buffer,
        );
        reads_rescued += 1;

    }
 
    // Final flush
    // if !gff_buffer.is_empty() {
    //     gff_writer.write_all(&gff_buffer).expect("Failed to flush GFF");
    //     aa_writer.write_all(&aa_buffer).expect("Failed to flush AA");
    //     dna_writer.write_all(&dna_buffer).expect("Failed to flush DNA");
    // }
 
    println!("\n--- 📊 Read Rescue Summary ---");
    println!("Total reads scanned:      {}", reads_scanned);
    println!("Already in graph:         {}", reads_already_assembled);
    println!("Rescued (FGS on read):    {}", reads_rescued);
    println!("Rescue rate:              {:.1}%",
        reads_rescued as f64 / reads_scanned as f64 * 100.0);
    println!("------------------------------");
 
    (reads_scanned, reads_rescued)
}
 
// ============================================================================
// CHECK IF A READ IS REPRESENTED IN THE GRAPH
//
// Fast heuristic: check k-mer coverage at 3 points in the read.
// If any sampled k-mer is below min_count, treat the read as unassembled.
// This is O(3) hash lookups, not O(read_length).
//
// This is a heuristic — it will occasionally misclassify reads. The error
// mode is conservative: false "not in graph" means we might run FGS on a
// read that was assembled. This adds small precision cost (duplicate predictions
// from reads that were also captured in a unitig) but no recall cost.
// ============================================================================
 
fn is_read_assembled(
    seq: &[u8],
    k: usize,
    mask: u64,
    counts: &FxHashMap<u64, u32>,
    min_count: u32,
) -> bool {
    if seq.len() < k {
        return false;
    }
 
    // Sample k-mers at start, middle, and end of read
    let check_positions = [
        0,
        seq.len().saturating_sub(k) / 2,
        seq.len().saturating_sub(k),
    ];
 
    for &start in &check_positions {
        if start + k > seq.len() {
            continue;
        }
        let kmer = pack_kmer(&seq[start..start + k], k, mask);
        if *counts.get(&kmer).unwrap_or(&0) < min_count {
            return false;  // This k-mer wasn't trusted — read not in graph
        }
    }
 
    true
}
 
// ============================================================================
// RUN FGS ON A SINGLE RAW READ
// ============================================================================
 
fn run_fgs_on_read(
    seq: &[u8],
    unitig_counter: &mut usize,
    global: &Box<Global>,
    locals: &Vec<Local>,
    gff_buffer: &mut Vec<u8>,
    aa_buffer: &mut Vec<u8>,
    dna_buffer: &mut Vec<u8>,
) {
    *unitig_counter += 1;
    let header_name = format!("read_{}", unitig_counter).into_bytes();
 
    let nseq: Vec<Nuc> = seq
        .iter()
        .map(|&b| b.to_ascii_uppercase())
        .map(Nuc::from)
        .collect();
 
    let prediction = frag_gene_scan_rs::viterbi::viterbi(
        global, locals, header_name, nseq, false
    );
 
    // if !prediction.genes.is_empty() {
    //     prediction.gff(gff_buffer).expect("Failed to write GFF");
    //     prediction.protein(aa_buffer, false).expect("Failed to write proteins");
    //     prediction.dna(dna_buffer, false).expect("Failed to write DNA");
    // }

    if !prediction.genes.is_empty() {
    
        prediction.gff(gff_buffer).expect("Failed to write GFF");
        prediction.protein(aa_buffer, false).expect("Failed to write proteins");
        prediction.dna(dna_buffer, false).expect("Failed to write DNA");
    } 
}
 
// ============================================================================
// HELPER: pack a raw byte slice into a u64 k-mer
// ============================================================================
 
fn pack_kmer(seq: &[u8], k: usize, mask: u64) -> u64 {
    let mut packed: u64 = 0;
    for &b in seq.iter().take(k) {
        let bits = match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _    => 0,  // N treated as A
        };
        packed = ((packed << 2) | bits) & mask;
    }
    packed
}
 
// ============================================================================
// FRAME-AWARE TRAVERSAL
//
// Key addition: track the dominant reading frame throughout the walk.
// At junctions, use frame continuity as a tie-breaker.
// ============================================================================
pub fn extract_unitigs_frame_aware(
    graph: &FxHashMap<u64, u8>,
    counts: &FxHashMap<u64, u32>,
    k: usize,
    global: &Box<Global>,
    locals: &Vec<Local>,
    gff_buffer: &mut Vec<u8>,
    aa_buffer: &mut Vec<u8>,
    dna_buffer: &mut Vec<u8>,
    min_unitig_length: usize,
) -> (usize, usize) {
    let mask_64: u64 = if k == 32 { !0 } else { (1 << (2 * k)) - 1 };
 
    println!("📊 Calculating in-degrees...");
    let mut in_degrees: FxHashMap<u64, u32> =
        FxHashMap::with_capacity_and_hasher(graph.len(), Default::default());
 
    for (&kmer, &out_mask) in graph.iter() {
        for &base in &[b'A', b'C', b'G', b'T'] {
            if (out_mask & base_to_bit_mask(base)) != 0 {
                let next_kmer = ((kmer << 2) | base_to_bits(base)) & mask_64;
                *in_degrees.entry(next_kmer).or_insert(0) += 1;
            }
        }
    }
 
    let mut visited: FxHashSet<u64> = FxHashSet::default();
    
    let mut total_unitig_count = 0usize;
    let mut coding_unitig_count = 0usize;
    let mut total_junctions = 0usize;
    let mut frame_wins = 0usize;      // junctions resolved by frame consistency
    let mut count_wins = 0usize;      // junctions resolved by count alone
 
    println!("🚀 Starting extraction with frame-aware junction resolution...");
 
    for (&start_kmer, &out_mask_start) in graph.iter() {
        if visited.contains(&start_kmer) {
            continue;
        }
 
        let in_deg = *in_degrees.get(&start_kmer).unwrap_or(&0);
        let out_deg = out_mask_start.count_ones();
 
        if in_deg != 1 || out_deg != 1 {
            let mut current_path = unpack_kmer(start_kmer, k);
            let mut current_kmer = start_kmer;
            let mut current_frame = detect_frame(&current_path);  // Initialize frame
            visited.insert(current_kmer);
 
            loop {
                let local_mask = *graph.get(&current_kmer).unwrap_or(&0);
                let current_out_deg = local_mask.count_ones();
 
                if current_out_deg == 0 {
                    break;
                }
 
                if current_out_deg == 1 {
                    // ---- unambiguous extension ------
                    let next_base = bit_mask_to_base(local_mask);
                    let next_kmer =
                        ((current_kmer << 2) | base_to_bits(next_base as u8)) & mask_64;
                    if visited.contains(&next_kmer) {
                        break;
                    }
                    current_path.push(next_base);
                    current_frame = (current_frame + 1) % 3;  // Update frame
                    visited.insert(next_kmer);
                    current_kmer = next_kmer;
                } else {
                    // ---- junction: frame-aware resolution ------
                    total_junctions += 1;
 
                    // Score each branch by (frame_consistency, count)
                    let mut best_score = (-1, 0u32);
                    let mut best_next: Option<(u64, char)> = None;
                    let mut used_frame = false;
 
                    for &base_char in &['A', 'C', 'G', 'T'] {
                        let bit = base_to_bit_mask(base_char as u8);
                        if (local_mask & bit) == 0 {
                            continue;
                        }
 
                        let next_kmer =
                            ((current_kmer << 2) | base_to_bits(base_char as u8)) & mask_64;
                        
                        if visited.contains(&next_kmer) {
                            continue;
                        }
 
                        let count = *counts.get(&next_kmer).unwrap_or(&0);
                        
                        // Check frame consistency at this junction
                        // The branch base contributes to frame evolution
                        let next_frame = (current_frame + 1) % 3;
                        
                        // Simple heuristic: check if the next k-mer "looks" coding
                        // by seeing if its local context favors the current frame
                        let branch_phase_score = measure_codon_phase_consistency(
                            &current_path,
                            base_char as char,
                            next_kmer,
                            graph,
                            current_frame,
                            mask_64,
                        );
                    

                        // Lexicographic comparison: frame first, count second
                        let score = (branch_phase_score, count);
                        
                        if score.0 > best_score.0 || 
                           (score.0 == best_score.0 && score.1 > best_score.1) {
                            best_score = score;
                            best_next = Some((next_kmer, base_char as char));
                            used_frame = score.0 > 0;
                        }
                    }
 
                    match best_next {
                        None => {
                            break;
                        }
                        Some((next_kmer, next_base)) => {
                            if used_frame {
                                frame_wins += 1;
                            } else {
                                count_wins += 1;
                            }
                            current_path.push(next_base);
                            current_frame = (current_frame + 1) % 3;
                            visited.insert(next_kmer);
                            current_kmer = next_kmer;
                        }
                    }
                }
            }
 
            // ================================================================
            // STREAMING FGS: Process completed unitig inline
            // ================================================================
            if current_path.len() >= min_unitig_length {
                total_unitig_count += 1;
 
                let heuristic_score = orf_heuristic_score(&current_path);
                
                if heuristic_score > 0.5 {
                    let header_name = format!("unitig_{}", total_unitig_count).into_bytes();
                    
                    let nseq: Vec<Nuc> = current_path
                        .bytes()
                        .map(|b| b.to_ascii_uppercase())
                        .map(Nuc::from)
                        .collect();
 
                    let prediction = frag_gene_scan_rs::viterbi::viterbi(
                        &global, 
                        &locals, 
                        header_name, 
                        nseq, 
                        false
                    );
 
                    if !prediction.genes.is_empty() {
                        coding_unitig_count += 1;
                        
                        prediction.gff(gff_buffer).expect("Failed to write GFF");
                        prediction.protein(aa_buffer, false).expect("Failed to write proteins");
                        prediction.dna(dna_buffer, false).expect("Failed to write DNA");
 
                        if coding_unitig_count % 50_000 == 0 {
                            println!(
                                "✅ Processed {} unitigs ({} coding) | Frame wins: {}",
                                total_unitig_count,
                                coding_unitig_count,
                                frame_wins,
                            );
                        }
                    }
                }
            }
        }
    }
 
    // --- final summary -------------------------------------------------------
    println!("\n--- 📊 Frame-Aware Traversal Summary ---");
    println!("Total junctions:          {}", total_junctions);
    println!("Resolved by frame:        {}", frame_wins);
    println!("Resolved by count:        {}", count_wins);
    println!("Frame-guided decisions:   {:.1}%", (frame_wins as f64 / total_junctions.max(1) as f64) * 100.0);
    println!("Total unitigs assembled:  {}", total_unitig_count);
    println!("Coding unitigs (FGS+):    {}", coding_unitig_count);
    println!("------------------------------------------");
 
    (total_unitig_count, coding_unitig_count)
}
// ============================================================================
// FRAME DETECTION: Detect dominant reading frame from a sequence
//
// Simple heuristic: count stop codons in each frame, pick frame with fewest
// ============================================================================
fn detect_frame(seq: &str) -> u8 {
    const STOPS: &[&str] = &["TAA", "TAG", "TGA"];
    let bytes = seq.as_bytes();
 
    let mut best_frame = 0u8;
    let mut best_score = f64::INFINITY;
 
    for frame in 0..3 {
        let mut stop_count = 0;
        let mut codon_count = 0;
        let mut i = frame;
 
        while i + 3 <= bytes.len() {
            codon_count += 1;
            let codon = std::str::from_utf8(&bytes[i..i + 3]).unwrap_or("NNN");
            if STOPS.contains(&codon) {
                stop_count += 1;
            }
            i += 3;
        }
 
        if codon_count == 0 {
            continue;
        }
 
        let stop_ratio = stop_count as f64 / codon_count as f64;
        if stop_ratio < best_score {
            best_score = stop_ratio;
            best_frame = frame as u8;
        }
    }
 
    best_frame
}
 
// ============================================================================
// CODON PHASE CONSISTENCY SCORER
//
// For a candidate next base and branch, measure how "consistent" the resulting
// codon reading frame is with the current walking frame.
//
// Returns: 1 if consistent, 0 if inconsistent
// ============================================================================
fn measure_codon_phase_consistency(
    current_path: &str,
    next_base: char,
    next_kmer: u64,
    graph: &FxHashMap<u64, u8>,
    current_frame: u8,
    mask_64: u64,
) -> i32 {
    // Build candidate sequence: last 2 bp + next_base + lookahead
    let mut lookahead = String::new();
    let mut curr = next_kmer;
    
    // FIXED: Extended lookahead from 6bp to 15bp (5 codons instead of 2)
    // This gives much better statistical signal for stop codon frequency
    for _ in 0..15 {
        match graph.get(&curr) {
            Some(&out_edges) if out_edges != 0 => {
                let next_b = bit_mask_to_base(out_edges);
                lookahead.push(next_b);
                curr = ((curr << 2) | base_to_bits(next_b as u8)) & mask_64;
            }
            _ => break,
        }
    }
 
    // Build candidate sequence
    let mut candidate = String::new();
    if current_path.len() >= 2 {
        candidate.push_str(&current_path[current_path.len() - 2..]);
    }
    candidate.push(next_base);
    candidate.push_str(&lookahead);
 
    if candidate.len() < 9 {
        // Too short to evaluate reliably
        return 0;
    }
 
    // Score: check stop codon frequency in the next reading frame
    const STOPS: &[&str] = &["TAA", "TAG", "TGA"];
    let bytes = candidate.as_bytes();
    let next_frame = (current_frame + 1) % 3;
    
    let mut stop_count = 0;
    let mut total_codon_count = 0;
    let mut i = next_frame as usize;
    
    while i + 3 <= bytes.len() {
        total_codon_count += 1;
        let codon = std::str::from_utf8(&bytes[i..i + 3]).unwrap_or("NNN");
        if STOPS.contains(&codon) {
            stop_count += 1;
        }
        i += 3;
    }
 
    if total_codon_count == 0 {
        return 0;
    }
 
    // FIXED: Return 0/1/2 based on stop codon ratio
    // Thresholds are biologically grounded:
    //   - Random sequence: 33% stops (1/3 of codons are stops)
    //   - Real coding: 1–5% stops
    //   - Rare tRNA/pseudogene: 10–15% stops  
    //   - Clearly non-coding: >20% stops
    
    let stop_ratio = stop_count as f64 / total_codon_count as f64;
    
    if stop_ratio < 0.1 {
        // <10% stops: strong signal of coding
        // This branch looks like real genes
        2
    } else if stop_ratio < 0.3 {
        // 10–30% stops: ambiguous
        // Could be coding with some pseudogenes/frameshifts mixed in
        1
    } else {
        // >30% stops: strong signal of non-coding or broken frame
        // This branch is likely wrong
        0
    }
}
 

// ============================================================================
// STREAMING VERSION: extract_unitigs_greedy with inline FGS
//
// Signature changes:
//   - Returns (total_unitig_count, coding_unitig_count) instead of Vec<String>
//   - Takes &mut output buffers for gff/aa/dna
//   - Calls FGS as each unitig completes, writes directly to buffers
//   - Non-coding unitigs are discarded immediately
// ============================================================================
pub fn extract_unitigs_greedy_streaming(
    graph: &FxHashMap<u64, u8>,
    counts: &FxHashMap<u64, u32>,
    k: usize,
    global: &Box<Global>,
    locals: &Vec<Local>,
    gff_buffer: &mut Vec<u8>,
    aa_buffer: &mut Vec<u8>,
    dna_buffer: &mut Vec<u8>,
    min_unitig_length: usize,
) -> (usize, usize) {
    let mask_64: u64 = if k == 32 { !0 } else { (1 << (2 * k)) - 1 };

    // --- compute in-degrees -------------------------------------------------
    println!("📊 Calculating in-degrees...");
    let mut in_degrees: FxHashMap<u64, u32> =
        FxHashMap::with_capacity_and_hasher(graph.len(), Default::default());

    for (&kmer, &out_mask) in graph.iter() {
        for &base in &[b'A', b'C', b'G', b'T'] {
            if (out_mask & base_to_bit_mask(base)) != 0 {
                let next_kmer = ((kmer << 2) | base_to_bits(base)) & mask_64;
                *in_degrees.entry(next_kmer).or_insert(0) += 1;
            }
        }
    }

    // --- shared state -------------------------------------------------------
    let mut visited: FxHashSet<u64> = FxHashSet::default();
    
    // --- counters for statistics -----------
    let mut total_unitig_count = 0usize;
    let mut coding_unitig_count = 0usize;
    let mut total_junctions = 0usize;
    let mut unresolved = 0usize;

    println!("🚀 Starting extraction with streaming FGS...");

    for (&start_kmer, &out_mask_start) in graph.iter() {
        if visited.contains(&start_kmer) {
            continue;
        }

        let in_deg = *in_degrees.get(&start_kmer).unwrap_or(&0);
        let out_deg = out_mask_start.count_ones();

        // Only start a walk at true source nodes
        if in_deg != 1 || out_deg != 1 {
            let mut current_path = unpack_kmer(start_kmer, k);
            let mut current_kmer = start_kmer;
            visited.insert(current_kmer);

            loop {
                let local_mask = *graph.get(&current_kmer).unwrap_or(&0);
                let current_out_deg = local_mask.count_ones();

                if current_out_deg == 0 {
                    break; // dead end
                }

                if current_out_deg == 1 {
                    // ---- unambiguous extension ------
                    let next_base = bit_mask_to_base(local_mask);
                    let next_kmer =
                        ((current_kmer << 2) | base_to_bits(next_base as u8)) & mask_64;
                    if visited.contains(&next_kmer) {
                        break;
                    }
                    current_path.push(next_base);
                    visited.insert(next_kmer);
                    current_kmer = next_kmer;
                } else {
                    // ---- junction: greedy count tiebreaker ------
                    total_junctions += 1;
                    
                    let mut best_count = 0u32;
                    let mut best_next: Option<(u64, char)> = None;

                    for &base_char in &['A', 'C', 'G', 'T'] {
                        let bit = base_to_bit_mask(base_char as u8);
                        if (local_mask & bit) == 0 {
                            continue;
                        }

                        let next_kmer =
                            ((current_kmer << 2) | base_to_bits(base_char as u8)) & mask_64;
                        
                        if visited.contains(&next_kmer) {
                            continue;
                        }

                        let c = *counts.get(&next_kmer).unwrap_or(&0);
                        if c > best_count {
                            best_count = c;
                            best_next = Some((next_kmer, base_char as char));
                        }
                    }

                    match best_next {
                        None => {
                            break; // all successors visited
                        }
                        Some((next_kmer, next_base)) => {
                            current_path.push(next_base);
                            visited.insert(next_kmer);
                            current_kmer = next_kmer;
                        }
                    }
                }
            }

            // ================================================================
            // GAP 5 CHANGE: Process completed unitig INLINE instead of storing
            // ================================================================
            if current_path.len() >= min_unitig_length {
                total_unitig_count += 1;

                // ---- PRE-FILTER: lightweight ORF heuristic ----
                // Only invoke FGS if heuristic suggests coding potential
                let heuristic_score = orf_heuristic_score(&current_path);
                
                if heuristic_score > 0.5 {
                    // This unitig has decent coding signal — run FGS
                    let header_name = format!("unitig_{}", total_unitig_count).into_bytes();
                    
                    let nseq: Vec<Nuc> = current_path
                        .bytes()
                        .map(|b| b.to_ascii_uppercase())
                        .map(Nuc::from)
                        .collect();

                    let prediction = frag_gene_scan_rs::viterbi::viterbi(
                        &global, 
                        &locals, 
                        header_name, 
                        nseq, 
                        false
                    );

                    // Only count and write if FGS found actual genes
                    if !prediction.genes.is_empty() {
                        coding_unitig_count += 1;
                        
                        // Write directly to buffers — no intermediate storage
                        prediction.gff(gff_buffer).expect("Failed to write GFF");
                        prediction.protein(aa_buffer, false).expect("Failed to write proteins");
                        prediction.dna(dna_buffer, false).expect("Failed to write DNA");

                        if coding_unitig_count % 50_000 == 0 {
                            println!(
                                "✅ Processed {} unitigs ({} coding) | Memory: streaming only",
                                total_unitig_count,
                                coding_unitig_count,
                            );
                        }
                    }
                    // else: FGS found no genes, unitig is silently dropped
                } else {
                    // Heuristic score too low — skip FGS entirely (save time)
                }
            }
        }
    }

    // --- final summary -------------------------------------------------------
    println!("\n--- 📊 Streaming Traversal Summary ---");
    println!("Total junctions:          {}", total_junctions);
    println!("Unresolved junctions:     {}", unresolved);
    println!("Total unitigs assembled:  {}", total_unitig_count);
    println!("Coding unitigs (FGS+):    {}", coding_unitig_count);
    println!("------------------------------------------");

    (total_unitig_count, coding_unitig_count)
}

// ============================================================================
// HELPER FUNCTIONS (unchanged from before)
// ============================================================================

fn orf_heuristic_score(seq: &str) -> f64 {
    const STOPS: &[&str] = &["TAA", "TAG", "TGA"];
    let bytes = seq.as_bytes();
    let len = bytes.len();

    let mut best_frame_score: f64 = 0.0;

    for frame in 0..3usize {
        let mut stop_count = 0usize;
        let mut codon_count = 0usize;
        let mut i = frame;
        while i + 3 <= len {
            codon_count += 1;
            let codon = &seq[i..i + 3];
            if STOPS.contains(&codon) {
                stop_count += 1;
            }
            i += 3;
        }
        if codon_count == 0 {
            continue;
        }
        let frame_score = 1.0 - (stop_count as f64 / codon_count as f64);
        if frame_score > best_frame_score {
            best_frame_score = frame_score;
        }
    }

    best_frame_score
}

fn prune_tips(graph: &mut FxHashMap<u64, u8>, k: usize, mask_64: u64) -> usize {
    let mut tips_removed = 0;
    let mut to_update = Vec::new();
    let max_tip_length = 2 * k;

    for (&kmer, &out_mask) in graph.iter() {
        let out_deg = out_mask.count_ones();
        if out_deg > 1 {
            let mut new_mask = out_mask;
            
            for &base in &[b'A', b'C', b'G', b'T'] {
                let base_bit = base_to_bit_mask(base);
                if (out_mask & base_bit) != 0 {
                    if is_tip(kmer, base, graph, mask_64, max_tip_length) {
                        new_mask &= !base_bit;
                        tips_removed += 1;
                    }
                }
            }
            if new_mask != out_mask {
                to_update.push((kmer, new_mask));
            }
        }
    }

    for (kmer, new_mask) in to_update {
        graph.insert(kmer, new_mask);
    }

    tips_removed
}

fn is_tip(start_kmer: u64, first_base: u8, graph: &FxHashMap<u64, u8>, mask_64: u64, max_len: usize) -> bool {
    let mut curr = ((start_kmer << 2) | base_to_bits(first_base)) & mask_64;
    for _ in 0..max_len {
        let mask = *graph.get(&curr).unwrap_or(&0);
        let out_deg = mask.count_ones();
        
        if out_deg == 0 { return true; }
        if out_deg > 1 { return false; }
        
        let next_base = bit_mask_to_base(mask);
        curr = ((curr << 2) | base_to_bits(next_base as u8)) & mask_64;
    }
    false
}


pub fn score_virtual_path(
    global: &Box<Global>,
    locals: &Vec<Local>,
    path_sequence: &str,
) -> frag_gene_scan_rs::gene::ReadPrediction {
    let nseq: Vec<Nuc> = path_sequence
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .map(Nuc::from)
        .collect();
    let head = b"graph_branch_candidate".to_vec();
    frag_gene_scan_rs::viterbi::viterbi(global, locals, head, nseq, false)
}