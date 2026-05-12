use needletail::{parse_fastx_file, Sequence};
//use rayon::prelude::*;
use dashmap::DashMap;
use std::sync::Arc;
use wyhash::wyhash;
use rustc_hash::FxHashMap;

use crate::utils::process_batch;

pub fn base_to_bits(base: u8) -> u64 {
    match base.to_ascii_uppercase() {
        b'C' => 1, b'G' => 2, b'T' => 3, _ => 0, // A=0, N/others=A
    }
}

pub fn base_to_bit_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1, b'C' => 2, b'G' => 4, b'T' => 8, _ => 0,
    }
}

// Inverse helper for the traversal
pub fn bit_mask_to_base(mask: u8) -> char {
    if mask & 1 != 0 { 'A' }
    else if mask & 2 != 0 { 'C' }
    else if mask & 4 != 0 { 'G' }
    else { 'T' }
}

pub fn unpack_kmer(mut packed: u64, k: usize) -> String {
    let mut chars = vec![0u8; k];
    for i in (0..k).rev() {
        let val = packed & 0b11;
        chars[i] = match val {
            0 => b'A', 1 => b'C', 2 => b'G', 3 => b'T', _ => unreachable!(),
        };
        packed >>= 2;
    }
    String::from_utf8(chars).unwrap()
}

pub fn count_kmers(
    input: &str,
    k: usize,
) -> (Arc<DashMap<u64, u32>>, Arc<DashMap<u64, Vec<u8>>>) {

    let counter = Arc::new(DashMap::new());
    let dna_map = Arc::new(DashMap::new());

    let mut reader = parse_fastx_file(input).expect("Invalid file");

    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid record");
        let seq = seqrec.normalize(false);

        for kmer in seq.windows(k) {
            let hash = wyhash(kmer, 0);

            *counter.entry(hash).or_insert(0) += 1;

            dna_map.entry(hash).or_insert_with(|| kmer.to_vec());
        }
    }

    println!("Num kmers before filter: {}", counter.len());
    counter.retain(|_, count| *count >= 1);
    println!("Num kmers after filter: {}", counter.len());

    (counter, dna_map)
}

pub fn count_minimizers(
        input: &str,
        k: usize,
        w: usize,
    ) -> (Arc<DashMap<u64, u32>>, Arc<DashMap<u64, Vec<u8>>>) {
    
        let counter = Arc::new(DashMap::new());
        let dna_map = Arc::new(DashMap::new());
    
        let mut reader = parse_fastx_file(input).expect("Invalid file");
        let mut batch = Vec::with_capacity(1000);
    
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid record");
            batch.push(seqrec.normalize(false).into_owned());
    
            if batch.len() >= 1000 {
                process_batch(&batch, &counter, &dna_map, k, w);
                batch.clear();
            }
        }
        process_batch(&batch, &counter, &dna_map, k, w);
    
        println!("Num minimizers before filter: {}", counter.len());
        counter.retain(|_, count| *count >= 2);
        println!("Num minimizers after filter: {}", counter.len());
        (counter, dna_map)
}
    

pub fn generate_minimizers(seq: &[u8], k: usize, w: usize) -> Vec<u64> {

    let hashes: Vec<u64> = seq.windows(k)
    .map(|kmer| {
        let canon = canonical_kmer(kmer);
        wyhash(&canon, 0)
    })
    .collect();
        //let hashes: Vec<u64> = seq.windows(k).map(|kmer| wyhash(kmer, 0)).collect();
    
        let mut minimizers: Vec<u64> = hashes
            .windows(w)
            .map(|window| *window.iter().min().unwrap())
            .collect();
    
        minimizers.dedup();
        minimizers
    }

pub fn canonical_kmer(kmer: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(kmer);

    if kmer <= rc.as_slice() {
        kmer.to_vec()
    } else {
        rc
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}


// try full k-mers instead of minimizers to see if that's the problem
// fn generate_kmer_edges(seq: &[u8], k: usize) -> Vec<(u64, u64)> {
//     let kmers: Vec<u64> = seq
//         .windows(k)
//         .map(|kmer| wyhash(kmer, 0))
//         .collect();

//     kmers
//         .windows(2)
//         .map(|pair| (pair[0], pair[1]))
//         .collect()
// }


// Added `unitig_id` so the header can trace back to the exact node in your graph
pub fn get_orfs_from_unitig(unitig_id: usize, dna: &str, min_aa_len: usize) -> Vec<(String, String)> {
    let mut orfs = Vec::new();
    let rev_dna = reverse_complement_str(dna);
    let sequences = [dna, &rev_dna];

    // Enumerate to track strand orientation (0 = Forward, 1 = Reverse)
    for (strand_idx, seq) in sequences.iter().enumerate() {
        let bytes = seq.as_bytes();
        let seq_len = bytes.len();
        
        // Define strand symbol for the header
        let strand_symbol = if strand_idx == 0 { '+' } else { '-' };

        // 3 reading frames per strand
        for frame_offset in 0..3 {
            let mut current_aa_seq = String::new();
            
            // Translate the frame
            for i in (frame_offset..seq_len).step_by(3) {
                if i + 3 <= seq_len {
                    let codon = &bytes[i..i + 3];
                    current_aa_seq.push(translate_codon(codon));
                }
            }

            // Extract valid ORFs
            let mut current_orf = String::new();
            let mut in_orf = false;
            let mut orf_counter = 0; // Track multiple ORFs in the same frame

            for ch in current_aa_seq.chars() {
                if ch == 'M' && !in_orf {
                    in_orf = true; // Found a start codon
                }
                
                if in_orf {
                    current_orf.push(ch);
                    if ch == '*' {
                        // Reached a stop codon -> Complete ORF
                        if current_orf.len() >= min_aa_len {
                            let header = format!(
                                ">unitig_{}_strand_{}_frame_{}_orf_{}", 
                                unitig_id, strand_symbol, frame_offset, orf_counter
                            );
                            orfs.push((header, current_orf.clone()));
                            orf_counter += 1;
                        }
                        current_orf.clear();
                        in_orf = false;
                    }
                }
            }

            // Capture "Open-Ended" ORFs (hit the branch before a stop codon)
            if in_orf && current_orf.len() >= min_aa_len {
                // Notice the "_partial" tag! This is the signal for your branch-resolution logic.
                let header = format!(
                    ">unitig_{}_strand_{}_frame_{}_orf_{}_partial", 
                    unitig_id, strand_symbol, frame_offset, orf_counter
                );
                orfs.push((header, current_orf));
            }
        }
    }
    orfs
}


fn translate_codon(codon: &[u8]) -> char {
    match codon {
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TGT" | b"TGC" => 'C',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TTT" | b"TTC" => 'F',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        b"CAT" | b"CAC" => 'H',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"AAA" | b"AAG" => 'K',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATG" => 'M',
        b"AAT" | b"AAC" => 'N',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"CAA" | b"CAG" => 'Q',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TGG" => 'W',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*', // Stop codons
        _ => 'X', // Unknown/Ambiguous
    }
}

fn reverse_complement_str(dna: &str) -> String {
    dna.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'N',
        })
        .collect()
}


use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Clone)]
pub struct BeamNode {
    pub kmer: u64,
    pub hmm_state: HMMState,
    pub score: f64,
    pub path: Vec<u8>,
    pub length: usize,
}

#[derive(Clone)]
pub struct ScoredNode(pub BeamNode);

impl Eq for ScoredNode {}

impl PartialEq for ScoredNode {
    fn eq(&self, other: &Self) -> bool {
        self.0.score == other.0.score
    }
}

impl Ord for ScoredNode {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.score.partial_cmp(&other.0.score).unwrap()
    }
}

impl PartialOrd for ScoredNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub struct BeamSearch {
    pub beam_width: usize,
    pub max_depth: usize,
}

// #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
// pub enum CodingState {
//     NonCoding,
//     Coding { frame: u8 }, // 0,1,2
// }

// #[derive(Clone, Debug)]
// pub struct HMMState {
//     pub coding: CodingState,
//     pub log_prob: f64,  // cumulative log-probability
// }

#[derive(Clone, Copy, Debug)]
pub enum CodingState {
    NonCoding,
    Coding { frame: u8 }, // 0,1,2
}

#[derive(Clone, Debug)]
pub struct HMMState {
    pub coding: CodingState,
    pub log_prob: f64,
    pub last_two: [u8; 2], // rolling buffer
}

// #[derive(Clone)]
// pub struct BeamNode {
//     pub kmer: u64,
//     pub hmm_state: HMMState,
//     pub score: f64,              // total score (HMM + other terms)
//     pub path: Vec<u8>,           // sequence (store as bytes for speed)
//     pub length: usize,

//     // Optional optimizations:
//     pub last_base: u8,           // avoid recomputing
// }

pub fn hmm_transition(
    prev: &HMMState,
    base: u8,
) -> HMMState {

    let mut next = prev.clone();

    // update rolling buffer
    let last_two = [prev.last_two[1], base];
    next.last_two = last_two;

    match prev.coding {

        // ---------------- NON-CODING ----------------
        CodingState::NonCoding => {
            if is_start_codon(prev.last_two, base) {
                // Enter coding frame
                next.coding = CodingState::Coding { frame: 0 };
                next.log_prob += 2.0; // reward start codon
            } else {
                next.log_prob += -0.1; // background penalty
            }
        }

        // ---------------- CODING ----------------
        CodingState::Coding { frame } => {
            let new_frame = (frame + 1) % 3;

            // Only evaluate stop codons when codon completes
            if frame == 2 && is_stop_codon(prev.last_two, base) {
                next.coding = CodingState::NonCoding;
                next.log_prob += 1.0; // reward proper stop
            } else {
                next.coding = CodingState::Coding { frame: new_frame };

                // Reward coding continuity
                next.log_prob += 0.5;

                // Penalize premature stops (out-of-frame noise)
                if is_stop_codon(prev.last_two, base) {
                    next.log_prob -= 2.0;
                }
            }
        }
    }

    next
}
pub fn expand_node(
    node: &BeamNode,
    graph: &FxHashMap<u64, u8>,
    mask: u64,
) -> Vec<BeamNode> {

    let mut children = Vec::new();

    if let Some(&out_mask) = graph.get(&node.kmer) {

        for &base in &[b'A', b'C', b'G', b'T'] {
            let bit = base_to_bit_mask(base);

            if (out_mask & bit) == 0 {
                continue;
            }

            let next_kmer = ((node.kmer << 2) | base_to_bits(base)) & mask;

            let new_hmm = hmm_transition(&node.hmm_state, base);

            let mut new_path = node.path.clone();
            new_path.push(base);

           //let new_score = node.score + new_hmm.log_prob;
            let delta = new_hmm.log_prob - node.hmm_state.log_prob;
            let new_score = node.score + delta;

            children.push(BeamNode {
                kmer: next_kmer,
                hmm_state: new_hmm,
                score: new_score,
                path: new_path,
                length: node.length + 1,
            });
        }
    }

    children
}

pub fn beam_search_from_kmer(
    start_kmer: u64,
    graph: &FxHashMap<u64, u8>,
    beam_width: usize,
    max_depth: usize,
    k: usize,
    mask: u64,
) -> Vec<BeamNode> {

    let initial_path = unpack_kmer(start_kmer, k).into_bytes();

    let initial_state = HMMState {
        coding: CodingState::NonCoding,
        log_prob: 0.0,
        last_two: [
            initial_path[initial_path.len() - 2],
            initial_path[initial_path.len() - 1],
        ],
    };

    let mut beam = vec![BeamNode {
        kmer: start_kmer,
        hmm_state: initial_state,
        score: 0.0,
        path: initial_path,
        length: k,
    }];

    // let mut beam: Vec<BeamNode> = vec![BeamNode {
    //     kmer: start_kmer,
    //     hmm_state: HMMState {
    //         coding: CodingState::NonCoding,
    //         log_prob: 0.0,
    //     },
    //     score: 0.0,
    //     path: unpack_kmer(start_kmer, k).into_bytes(),
    //     length: k,
    //     last_base: b'N',
    // }];

    for _step in 0..max_depth {
        let mut candidates: Vec<BeamNode> = Vec::new();

        for node in &beam {
            let children = expand_node(node, graph, mask);
            candidates.extend(children);
        }

        // prune to top-K
        candidates.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
        beam = candidates.into_iter().take(beam_width).collect();

        if beam.is_empty() {
            break;
        }
    }

    beam
}


#[inline]
fn is_start_codon(last_two: [u8; 2], next: u8) -> bool {
    matches!(
        (last_two[0], last_two[1], next),
        (b'A', b'T', b'G') |
        (b'G', b'T', b'G') |
        (b'T', b'T', b'G')
    )
}

#[inline]
fn is_stop_codon(last_two: [u8; 2], next: u8) -> bool {
    matches!(
        (last_two[0], last_two[1], next),
        (b'T', b'A', b'A') |
        (b'T', b'A', b'G') |
        (b'T', b'G', b'A')
    )
}