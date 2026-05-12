//! De Bruijn graph construction

use rustc_hash::FxHashMap;
use seq_io::fastq::{Reader, Record};
use crate::sequence::{base_to_bits, base_to_bit_mask};
use std::path::Path;

pub type KmerCount = FxHashMap<u64, u32>;
pub type KmerGraph = FxHashMap<u64, u8>;

/// De Bruijn graph structure
pub struct Graph {
    pub kmers: KmerGraph,
    pub counts: KmerCount,
    pub k: usize,
    pub mask: u64,
}

impl Graph {
    /// Create a new empty graph
    pub fn new(k: usize) -> Self {
        let mask = if k == 32 { !0 } else { (1u64 << (2 * k)) - 1 };
        Graph {
            kmers: FxHashMap::default(),
            counts: FxHashMap::default(),
            k,
            mask,
        }
    }

    /// Count k-mers from FASTQ file (Pass 1)
    pub fn count_kmers<P: AsRef<Path>>(&mut self, path: P, min_count: u32) -> std::io::Result<usize> {
        log::info!("Pass 1: Counting k-mers with k={}", self.k);

        let mut reader = Reader::from_path(path)?;
        let mut read_count = 0usize;

        while let Some(result) = reader.next() {
            let record = result.expect("Error reading record");
            let seq = record.seq();
            read_count += 1;

            if seq.len() < self.k {
                continue;
            }

            let mut current_packed: u64 = 0;
            for i in 0..self.k {
                current_packed = (current_packed << 2) | base_to_bits(seq[i]);
            }
            *self.counts.entry(current_packed).or_insert(0) += 1;

            for i in self.k..seq.len() {
                current_packed = ((current_packed << 2) | base_to_bits(seq[i])) & self.mask;
                *self.counts.entry(current_packed).or_insert(0) += 1;
            }
        }

        log::info!("Read {} sequences", read_count);
        log::info!("Total k-mers before filtering: {}", self.counts.len());

        self.counts.retain(|_, &mut count| count >= min_count);
        log::info!("K-mers after filtering (min_count={}): {}", min_count, self.counts.len());

        Ok(read_count)
    }

    /// Build de Bruijn graph from FASTQ file (Pass 2)
    pub fn build_graph<P: AsRef<Path>>(&mut self, path: P, min_count: u32) -> std::io::Result<()> {
        log::info!("Pass 2: Building de Bruijn graph");

        let mut reader = Reader::from_path(path)?;
        while let Some(result) = reader.next() {
            let record = result.expect("Error reading record");
            let seq = record.seq();

            if seq.len() < self.k + 1 {
                continue;
            }

            let mut current_packed: u64 = 0;
            for i in 0..self.k {
                current_packed = (current_packed << 2) | base_to_bits(seq[i]);
            }

            for i in self.k..seq.len() {
                let next_base = seq[i];
                let next_packed = ((current_packed << 2) | base_to_bits(next_base)) & self.mask;

                if *self.counts.get(&current_packed).unwrap_or(&0) >= min_count
                    && *self.counts.get(&next_packed).unwrap_or(&0) >= min_count
                {
                    let entry = self.kmers.entry(current_packed).or_insert(0);
                    *entry |= base_to_bit_mask(next_base);
                }
                current_packed = next_packed;
            }
        }

        log::info!("Graph size: {} nodes", self.kmers.len());
        Ok(())
    }

    /// Prune tips from the graph (short dead-end branches)
    pub fn prune_tips(&mut self) -> usize {
        log::info!("Pruning tips from graph...");

        let max_tip_length = 2 * self.k;
        let mut tips_removed = 0;
        let mut to_update = Vec::new();

        for (&kmer, &out_mask) in self.kmers.iter() {
            if out_mask.count_ones() > 1 {
                let mut new_mask = out_mask;

                for &base in &[b'A', b'C', b'G', b'T'] {
                    let base_bit = base_to_bit_mask(base);
                    if (out_mask & base_bit) != 0 && self.is_tip(kmer, base, max_tip_length) {
                        new_mask &= !base_bit;
                        tips_removed += 1;
                    }
                }

                if new_mask != out_mask {
                    to_update.push((kmer, new_mask));
                }
            }
        }

        for (kmer, new_mask) in to_update {
            self.kmers.insert(kmer, new_mask);
        }

        log::info!("Removed {} tips", tips_removed);
        tips_removed
    }

    /// Check if a branch from a k-mer is a tip (dead-end)
    fn is_tip(&self, start_kmer: u64, first_base: u8, max_len: usize) -> bool {
        use crate::sequence::{base_to_bits, bit_mask_to_base};

        let mut curr = ((start_kmer << 2) | base_to_bits(first_base)) & self.mask;
        for _ in 0..max_len {
            let mask = *self.kmers.get(&curr).unwrap_or(&0);
            let out_deg = mask.count_ones();

            if out_deg == 0 {
                return true;
            }
            if out_deg > 1 {
                return false;
            }

            let next_base = bit_mask_to_base(mask);
            curr = ((curr << 2) | base_to_bits(next_base as u8)) & self.mask;
        }
        false
    }

    /// Calculate in-degrees for all k-mers
    pub fn calculate_in_degrees(&self) -> FxHashMap<u64, u32> {
        use crate::sequence::{base_to_bit_mask, base_to_bits};

        log::debug!("Calculating in-degrees...");
        let mut in_degrees: FxHashMap<u64, u32> =
            FxHashMap::with_capacity_and_hasher(self.kmers.len(), Default::default());

        for (&kmer, &out_mask) in self.kmers.iter() {
            for &base in &[b'A', b'C', b'G', b'T'] {
                if (out_mask & base_to_bit_mask(base)) != 0 {
                    let next_kmer = ((kmer << 2) | base_to_bits(base)) & self.mask;
                    *in_degrees.entry(next_kmer).or_insert(0) += 1;
                }
            }
        }

        in_degrees
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_creation() {
        let graph = Graph::new(31);
        assert_eq!(graph.k, 31);
        assert_eq!(graph.kmers.len(), 0);
        assert_eq!(graph.counts.len(), 0);
    }
}