use std::path::PathBuf;

use clap::Parser;

#[derive(Debug, Parser)]
pub struct Args {
    #[arg(long, help = "Path to FASTQ file (reads).")]
    pub fastq: PathBuf,

    #[arg(long, help = "Path to FASTA file (to find coverage for).")]
    pub fasta: PathBuf,

    #[arg(short, long, help = "Kmer size", default_value_t = 15)]
    pub kmer_size: usize,

    #[arg(
        short,
        long,
        help = "Min kmer coverage to count as valid.",
        default_value_t = 1
    )]
    pub min_kmer_coverage: usize,

    #[arg(short, long, help = "Output file")]
    pub outfile: PathBuf,
}
