use crate::errors::AppError;
use bio_utils_rs::io::bio_fasta_reader;
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use packed_seq::{PackedSeqVec, SeqVec};
use rayon::prelude::*;
use rustc_hash::FxBuildHasher;
use simd_minimizers::canonical_minimizers;
use std::path::PathBuf;
use std::time::Duration;

/// There are different ways to build a minimizer index.
/// This implementation uses a DashMap<u64, Vec<(usize, usize)>> to store
/// information as DashMap<kmer_hash, Vec<(contig_id, contig_location)>>.
pub fn build_minimizer_index(
    fasta: PathBuf,
    kmer_size: usize,
) -> Result<
    (
        DashMap<u64, Vec<(usize, usize)>, FxBuildHasher>,
        DashMap<usize, (String, Vec<usize>), FxBuildHasher>,
    ),
    AppError,
> {
    let fasta_reader = bio_fasta_reader(Some(fasta))?;

    let minimizer_index: DashMap<u64, Vec<(usize, usize)>, FxBuildHasher> =
        DashMap::with_hasher(FxBuildHasher);

    let contig_index: DashMap<usize, (String, Vec<usize>), FxBuildHasher> =
        DashMap::with_hasher(FxBuildHasher);

    // Extract minimizers.
    let spinner = ProgressBar::new_spinner();
    spinner.enable_steady_tick(Duration::from_millis(200));
    spinner.set_style(ProgressStyle::with_template(
        "Loading fasta and building index {spinner:.blue} [{elapsed_precise}]",
    )?);

    fasta_reader
        .records()
        .enumerate()
        .par_bridge()
        .for_each(|(i, record)| {
            let record = match record {
                Ok(record) => record,
                Err(_) => return,
            };

            // If we use window_size = 1, we include all canonical kmers. For a contig of length l,
            // we can generate a maximum of l-k+1 kmers.
            let pos_vec = vec![0; record.seq().len() - kmer_size + 1];
            contig_index.insert(i, (record.id().to_string(), pos_vec));

            //
            let mut _simd_positions = vec![];

            let seq = PackedSeqVec::from_ascii(record.seq());

            // NOTE - using w=1 means extract ALL canonical kmers, which is what we want.
            for (pos, minimizer) in canonical_minimizers(kmer_size, 1)
                .run(seq.as_slice(), &mut _simd_positions)
                .pos_and_values_u64()
            {
                minimizer_index
                    .entry(minimizer)
                    .and_modify(|contig_vec| contig_vec.push((i, pos as usize)))
                    .or_insert(vec![(i, pos as usize)]);
            }
        });

    spinner.finish_and_clear();

    Ok((minimizer_index, contig_index))
}
