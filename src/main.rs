//! FRAME: Frame Resolved Assembly for Metagenomics

use clap::Parser;
use frame::{
    config::PipelineConfig, graph::Graph, io::PipelineOutput, prediction, traversal,
    utils::{self, ReadSummary},
    reader::open_sequence_file,
};
use log::info;
use seq_io::fastq::{Reader, Record};
use std::path::PathBuf;
use std::time::Instant;
use frag_gene_scan_rs::hmm::{Global, Local};


use rayon::prelude::*;


#[derive(Default)]
struct ThreadResult {
    gff: Vec<u8>,
    aa: Vec<u8>,
    dna: Vec<u8>,
    header_map: Vec<u8>,
    total_reads: usize,
    assembled: usize,
    rescued: usize,
}

#[derive(Default)]
struct UnitigResult {
    gff: Vec<u8>,
    aa: Vec<u8>,
    dna: Vec<u8>,
    coding_count: usize,
}


#[derive(Parser, Debug)]
#[command(
    name = "FRAME",
    about = "Frame Resolved Assembly for Metagenomics",
    long_about = "A de novo metagenomic assembler with frame-aware unitig extraction and gene prediction"
)]
struct Args {
    /// Input FASTQ/FASTA file
    // #[arg(value_name = "fasta/fastq file", help = "Path to input sequencing reads")]
    // input: PathBuf,
     /// Input FASTQ/FASTA file (R1 for paired-end)
     #[arg(value_name = "fasta/fastq file", help = "Path to input sequencing reads (can be gzipped)")]
     input: PathBuf,

     // Second input file (R2 for paired-end reads)
    #[arg(
        long,
        value_name = "fasta/fastq file",
        help = "Path to second input file for paired-end reads (can be gzipped)"
    )]
    input2: Option<PathBuf>,
    /// Output directory
    #[arg(
        short,
        long,
        value_name = "dir",
        default_value = "./frame_output",
        help = "Output directory for results"
    )]
    output: PathBuf,

    /// K-mer size
    #[arg(short, long, default_value = "31", help = "K-mer size (15-63)")]
    kmer: usize,

    /// Minimum k-mer count
    #[arg(
        short,
        long,
        default_value = "2",
        help = "Minimum k-mer count threshold"
    )]
    min_count: u32,

    /// Minimum unitig length
    #[arg(
        short,
        long,
        default_value = "100",
        help = "Minimum unitig length in bp"
    )]
    min_length: usize,

    /// HMM model directory
    #[arg(
        long,
        value_name = "dir",
        default_value = "./lib/FragGeneScanRs/train",
        help = "Path to HMM training directory"
    )]
    hmm_dir: PathBuf,

    /// HMM model name
    #[arg(
        long,
        value_name = "name",
        default_value = "illumina_5",
        help = "HMM model name"
    )]
    hmm_model: PathBuf,


    #[arg(short = 't', 
    long, 
    default_value = "16", 
    help = "Threads (default=16)")]
    threads: usize,
    }

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();

    info!("═══════════════════════════════════════════════════════");
    info!("FRAME: Frame Resolved Assembly for Metagenomics");
    info!("═══════════════════════════════════════════════════════\n");
  
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    let start_time = Instant::now();

    let config = PipelineConfig {
        k: args.kmer,
        min_count: args.min_count,
        min_unitig_length: args.min_length,
        hmm_train_dir: args.hmm_dir,
        hmm_model: args.hmm_model,
        input_file: args.input.clone(),
        input_file_2: args.input2.clone(),
        output_dir: args.output.clone(),
        rescue_batch_size: 50_000,
    };

    config.validate()?;
    info!("Configuration validated");
    info!("  k-mer size: {}", config.k);
    info!("  min count: {}", config.min_count);
    info!("  min unitig length: {}", config.min_unitig_length);
    info!("  Input file: {}", config.input_file.display());
    if let Some(ref f2) = config.input_file_2 {
        info!("  Paired-end mode: {}", f2.display());
    }
    info!("  Output directory: {}\n", config.output_dir.display());

    // Load HMM model
    let (global, locals) = prediction::load_hmm_model(
        config.hmm_train_dir.clone(),
        config.hmm_model.clone(),
    )?;

    // Phase 1: Count k-mers
    let phase1_start = Instant::now();
    let mut graph = Graph::new(config.k);
    graph.count_kmers(&config.input_file, config.input_file_2.as_ref(), config.min_count)?;
    //graph.count_kmers(&config.input_file, config.min_count)?;
    let phase1_time = phase1_start.elapsed();
    info!(
        "Phase 1 completed in {:.2}s: {} k-mers after filtering",
        phase1_time.as_secs_f64(),
        graph.counts.len()
    );

    // Phase 2: Build de Bruijn graph
    let phase2_start = Instant::now();
    //graph.build_graph(&config.input_file, config.min_count)?;
    graph.build_graph(&config.input_file, config.input_file_2.as_ref(), config.min_count)?;
    let phase2_time = phase2_start.elapsed();
    info!(
        "Phase 2 completed in {:.2}s: {} nodes in graph",
        phase2_time.as_secs_f64(),
        graph.kmers.len()
    );

    // Phase 3: Tip pruning and bubble removal
    let phase3_start = Instant::now();
    let tips_removed = graph.prune_tips();
    let phase3_time = phase3_start.elapsed();
    info!(
        "Phase 3 completed in {:.2}s: {} tips removed",
        phase3_time.as_secs_f64(),
        tips_removed
    );

    // Phase 3b: Bubble removal
    let phase3b_start = Instant::now();
    let bubbles_removed = graph.remove_bubbles();
    let phase3b_time = phase3b_start.elapsed();
    info!(
        "Bubble removal completed in {:.2}s: {} bubbles removed",
        phase3b_time.as_secs_f64(),
        bubbles_removed
    );

    // // Phase 4: Unitig extraction
    // let phase4_start = Instant::now();
    // let (unitigs, traversal_stats) =
    //     traversal::extract_unitigs_frame_aware(&graph.kmers, &graph.counts, config.k, config.min_unitig_length);

    // let mut output = PipelineOutput::new();
    // let mut coding_unitigs = 0;

    // for (idx, unitig) in unitigs.iter().enumerate() {
    //     let unitig_id = idx + 1;
    //     let heuristic_score = traversal::orf_heuristic_score(unitig);

    //     if heuristic_score > 0.5 {
    //         if let Ok(is_coding) = prediction::predict_and_write_unitig(
    //             unitig_id,
    //             unitig,
    //             &global,
    //             &locals,
    //             &mut output.gff_buffer,
    //             &mut output.aa_buffer,
    //             &mut output.dna_buffer,
    //         ) {
    //             if is_coding {
    //                 coding_unitigs += 1;
    //             }
    //         }
    //     }

    //     if coding_unitigs % 100_000 == 0 && coding_unitigs > 0 {
    //         info!("Processed {} coding unitigs...", coding_unitigs);
    //     }
    // }

    // let phase4_time = phase4_start.elapsed();
    // info!(
    //     "Phase 4 completed in {:.2}s: {} coding unitigs found",
    //     phase4_time.as_secs_f64(),
    //     coding_unitigs
    // );

    let phase4_start = Instant::now();

    let (unitigs, traversal_stats) =
        traversal::extract_unitigs_frame_aware(
            &graph.kmers,
            &graph.counts,
            config.k,
            config.min_unitig_length,
        );

    let results = unitigs
    .par_iter()
    .enumerate()
    .fold(
        || UnitigResult::default(),
        |mut local, (idx, unitig)| {
            let unitig_id = idx + 1;

            let heuristic_score =
                traversal::orf_heuristic_score(unitig);

            if heuristic_score <= 0.5 {
                return local;
            }

            if let Ok(is_coding) =
                prediction::predict_and_write_unitig(
                    unitig_id,
                    unitig,
                    &global,
                    &locals,
                    &mut local.gff,
                    &mut local.aa,
                    &mut local.dna,
                )
            {
                if is_coding {
                    local.coding_count += 1;
                }
            }

            local
        },
    )
    .reduce(
        || UnitigResult::default(),
        |mut a, b| {
            a.gff.extend(b.gff);
            a.aa.extend(b.aa);
            a.dna.extend(b.dna);
            a.coding_count += b.coding_count;
            a
        },
    );

    // Merge final outputs
    let mut output = PipelineOutput::new();

    output.gff_buffer = results.gff;
    output.aa_buffer = results.aa;
    output.dna_buffer = results.dna;

    // Write output
    if !output.is_empty() {
        output.write_assembly(&config.output_dir)?;
    }

    output.gff_buffer.clear();
    output.aa_buffer.clear();
    output.dna_buffer.clear();


    let coding_unitigs = results.coding_count;

    let phase4_time = phase4_start.elapsed();

    info!(
        "Phase 4 completed in {:.2}s: {} coding unitigs found",
        phase4_time.as_secs_f64(),
        coding_unitigs
    );

    // Phase 5: Read rescue
    let mut rescue_summary = ReadSummary {
        total_reads: 0,
        reads_assembled: 0,
        reads_rescued: 0,
    };

    
    let phase5_start = Instant::now();
    info!("Phase 5: Rescuing unmapped reads...");

    // let mut reader = Reader::from_path(&config.input_file)?;
    // let all_seqs: Vec<_> = reader
    //     .into_records()  // Collect all
    //     .collect::<Result<Vec<_>, _>>()?;

    // all_seqs
    //     .par_iter()  // Process in parallel
    //     .for_each(|record| {
    //         // Same prediction code
    //         predict_genes(record.seq(), ...);
    //     });

    // let mut reader = Reader::from_path(&config.input_file)?;
    // let mut unitig_counter = coding_unitigs + 1;

    // while let Some(result) = reader.next() {
    //     let record = result.expect("Error reading record");
    //     let seq = record.seq();
    //     rescue_summary.total_reads += 1;

    //     let is_assembled = utils::is_read_assembled(
    //         seq,
    //         config.k,
    //         config.get_mask(),
    //         &graph.counts,
    //         config.min_count,
    //     );

    //     if is_assembled {
    //         rescue_summary.reads_assembled += 1;
    //     } else {
    //         if let Ok(is_coding) = prediction::predict_and_write_read(
    //             unitig_counter,
    //             seq,
    //             &global,
    //             &locals,
    //             &mut output.gff_buffer,
    //             &mut output.aa_buffer,
    //             &mut output.dna_buffer,
    //         ) {
    //             if is_coding {
    //                 rescue_summary.reads_rescued += 1;
    //                 unitig_counter += 1;
    //             }
    //         }
    //     }
    // }


    //let mut reader = Reader::from_path(&config.input_file)?;
    let reader_io = open_sequence_file(&config.input_file)?;
    let reader = Reader::new(reader_io);
    let mut all_records: Vec<_> = reader
        .into_records()
        .collect::<Result<Vec<_>, _>>()?;

    // Load second file if provided (paired-end)
    if let Some(ref input2) = config.input_file_2 {
        let reader_io = open_sequence_file(input2)?;
        let reader = Reader::new(reader_io);
        let records: Vec<_> = reader
            .into_records()
            .collect::<Result<Vec<_>, _>>()?;
        all_records.extend(records);
    }

    
    // let all_records: Vec<_> = reader
    //     .into_records()
    //     .collect::<Result<Vec<_>, _>>()?;

    let results = all_records
        .par_iter()
        .enumerate()
        .fold(
        || ThreadResult::default(),
        |mut local, (i, record)| {
            let seq = record.seq();

            local.total_reads += 1;

            let is_assembled = utils::is_read_assembled(
                seq,
                config.k,
                config.get_mask(),
                &graph.counts,
                config.min_count,
            );

            if is_assembled {
                local.assembled += 1;
                return local;
            }

            let read_id = coding_unitigs + i;

            //let _ = prediction::predict_and_write_read(
            let _ = prediction::predict_and_write_read_with_header(
                read_id,
                &record.head,
                seq,
                &global,
                &locals,
                &mut local.gff,
                &mut local.aa,
                &mut local.dna,
                &mut local.header_map,
            )
            .map(|is_coding| {
                if is_coding {
                    local.rescued += 1;
                }
            });

            local
        },
    )
    .reduce(
        || ThreadResult::default(),
        |mut a, b| {
            a.gff.extend(b.gff);
            a.aa.extend(b.aa);
            a.dna.extend(b.dna);
            a.header_map.extend(b.header_map);

            a.total_reads += b.total_reads;
            a.assembled += b.assembled;
            a.rescued += b.rescued;

            a
        },
    );

    // Final merge
    output.gff_buffer = results.gff;
    output.aa_buffer = results.aa;
    output.dna_buffer = results.dna;
    output.header_map_buffer = results.header_map;

    // write rescue reads
    if !output.is_empty() {
        output.write_rescue(&config.output_dir)?;
    }



    rescue_summary.total_reads = results.total_reads;
    rescue_summary.reads_assembled = results.assembled;
    rescue_summary.reads_rescued = results.rescued;

    let phase5_time = phase5_start.elapsed();
    info!(
        "Phase 5 completed in {:.2}s: {} reads rescued",
        phase5_time.as_secs_f64(),
        rescue_summary.reads_rescued
    );

    // Write output
    if !output.is_empty() {
        output.write_predictions(&config.output_dir)?;
    }


    let elapsed = start_time.elapsed();

    // Print summary
    //utils::print_statistics(&rescue_summary, traversal_stats.total_unitigs, coding_unitigs);

    info!("═══════════════════════════════════════════════════════");
    info!("Total execution time: {:.2}s", elapsed.as_secs_f64());
    info!("Output size: {}", utils::format_file_size(output.total_size() as u64));
    info!("═══════════════════════════════════════════════════════\n");

    Ok(())
}