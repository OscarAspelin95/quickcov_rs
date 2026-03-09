use crate::errors::AppError;
use bio_utils_rs::io::needletail_reader;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::slice::ParallelSliceMut;
use simd_minimizers::canonical_minimizers;
use std::path::PathBuf;
use std::time::Duration;

pub fn build_minimizer_index(
    fasta: PathBuf,
    kmer_size: usize,
) -> Result<(Vec<(u64, usize, u32)>, Vec<usize>), AppError> {
    let mut fasta_reader = needletail_reader(Some(fasta))?;

    //
    let mut minimizer_index: Vec<(u64, usize, u32)> = vec![];
    let mut contig_index = 0;
    let mut contig_lengths: Vec<usize> = vec![];

    // Extract minimizers.
    let spinner = ProgressBar::new_spinner();
    spinner.enable_steady_tick(Duration::from_millis(200));
    spinner.set_style(ProgressStyle::with_template(
        "Loading fasta and building index {spinner:.blue} [{elapsed_precise}]",
    )?);

    while let Some(record) = fasta_reader.next() {
        let record = match record {
            Ok(record) => record,
            Err(_) => continue,
        };

        contig_index += 1;
        contig_lengths.push(record.num_bases());

        // extract minimizers.
        let mut _simd_positions = vec![];

        // NOTE - using w=1 means extract ALL canonical kmers, which is what we want.
        for (pos, minimizer) in canonical_minimizers(kmer_size, 1)
            .run(record.seq().iter().as_slice(), &mut _simd_positions)
            .pos_and_values_u64()
        {
            minimizer_index.push((minimizer, contig_index, pos));
        }
    }
    spinner.finish_and_clear();

    // Sort and reformat.
    let spinner = ProgressBar::new_spinner();
    spinner.enable_steady_tick(Duration::from_millis(200));
    spinner.set_style(ProgressStyle::with_template(
        "Sorting and re-formatting index {spinner:.blue} [{elapsed_precise}]",
    )?);

    // -- sort
    minimizer_index.par_sort_by_key(|a| a[0]);

    let total = minimizer_index.len();
    let mut kmer_hashes: Vec<u64> = Vec::new();
    let mut offsets: Vec<usize> = Vec::new();
    let mut entries: Vec<u64> = Vec::with_capacity(total);

    let mut prev_hash: u64 = 0;

    for (i, (minimizer, contig_id, contig_pos)) in minimizer_index.into_iter().enumerate() {
        if i == 0 || minimizer != prev_hash {
            kmer_hashes.push(minimizer);
            offsets.push(i);
            prev_hash = minimizer;
        }
        entries.push((contig_id as u64) << 32 | contig_pos as u64);
    }
    offsets.push(entries.len());

    spinner.finish_and_clear();
    Ok((kmer_hashes, offsets, entries, contig_lengths))
}
