use std::fs::File;
use std::io::{Write, Result as IoResult};
use std::path::Path;


pub struct PipelineOutput {
    pub gff_buffer: Vec<u8>,
    pub aa_buffer: Vec<u8>,
    pub dna_buffer: Vec<u8>,
}

impl PipelineOutput {
    pub fn new() -> Self {
        PipelineOutput {
            gff_buffer: Vec::with_capacity(1024 * 1024), // 1MB initial
            aa_buffer: Vec::with_capacity(1024 * 1024),
            dna_buffer: Vec::with_capacity(1024 * 1024),
        }
    }

    /// Write final predictions to files
    pub fn write_predictions<P: AsRef<Path>>(
        &self,
        output_dir: P,
    ) -> IoResult<()> {
        let dir = output_dir.as_ref();
        if !dir.exists() {
            std::fs::create_dir_all(dir)?;
        }

        let mut gff_file = File::create(dir.join("predictions.gff"))?;
        gff_file.write_all(&self.gff_buffer)?;

        let mut aa_file = File::create(dir.join("proteins.faa"))?;
        aa_file.write_all(&self.aa_buffer)?;

        let mut dna_file = File::create(dir.join("genes.fna"))?;
        dna_file.write_all(&self.dna_buffer)?;

        log::info!("Predictions written to {}", dir.display());
        Ok(())
    }

    /// Write assembly results to files
    pub fn write_assembly<P: AsRef<Path>>(
        &self,
        output_dir: P,
    ) -> IoResult<()> {
        let dir = output_dir.as_ref();
        if !dir.exists() {
            std::fs::create_dir_all(dir)?;
        }

        let mut gff_file = File::create(dir.join("assembly_predictions.gff"))?;
        gff_file.write_all(&self.gff_buffer)?;

        let mut aa_file = File::create(dir.join("assembly_proteins.faa"))?;
        aa_file.write_all(&self.aa_buffer)?;

        let mut dna_file = File::create(dir.join("assembly_dna.fna"))?;
        dna_file.write_all(&self.dna_buffer)?;

        log::info!("Assembly results written to {}", dir.display());
        Ok(())
    }

    /// Write rescue results to files
    pub fn write_rescue<P: AsRef<Path>>(
        &self,
        output_dir: P,
    ) -> IoResult<()> {
        let dir = output_dir.as_ref();
        if !dir.exists() {
            std::fs::create_dir_all(dir)?;
        }

        let mut gff_file = File::create(dir.join("rescue_predictions.gff"))?;
        gff_file.write_all(&self.gff_buffer)?;

        let mut aa_file = File::create(dir.join("rescue_proteins.faa"))?;
        aa_file.write_all(&self.aa_buffer)?;

        let mut dna_file = File::create(dir.join("rescue_dna.fna"))?;
        dna_file.write_all(&self.dna_buffer)?;

        log::info!("Rescue results written to {}", dir.display());
        Ok(())
    }

    /// Check if any results were produced
    pub fn is_empty(&self) -> bool {
        self.gff_buffer.is_empty() && self.aa_buffer.is_empty() && self.dna_buffer.is_empty()
    }

    /// Get total buffer size in bytes
    pub fn total_size(&self) -> usize {
        self.gff_buffer.len() + self.aa_buffer.len() + self.dna_buffer.len()
    }
}

impl Default for PipelineOutput {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_output_creation() {
        let output = PipelineOutput::new();
        assert!(output.is_empty());
        assert_eq!(output.total_size(), 0);
    }
}