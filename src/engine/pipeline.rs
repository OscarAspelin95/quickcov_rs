use crate::{args::Args, errors::AppError};
use bio_utils_rs::io::bio_fastq_reader;
use bio_utils_rs::io::needletail_reader;
use indicatif::{ProgressBar, ProgressStyle};
use packed_seq::{PackedSeqVec, SeqVec};
use simd_minimizers::canonical_minimizers;
use std::time::Duration;

/// We'll use the minimap2 approach for storing kmers.
/// This method avoids a hashmap all together and instead uses vectors.
/// 1. We'll iterate over the contigs and extract tuples of type (kmer_hash, contig_id, contig_loc).
/// e.g., (234329867, 1, 34) means kmer hash 234329867 exists in contig 1 at location 34.
///
/// Doing this for all contigs and kmers results in a vec like
/// [(kmer_hash, contig_id, contig_location), ...,].
///
/// We'll sort this by kmer_hash first (expensive, but we can do it with rayon par_sort)
/// Then, we'll initialize three arrays/vecs:
/// 	[h1, h2, ..., hn] // sorted, deduplicated kmer hashes
/// 	[0, 5, ..., num_unique_(contig_id, location)] // offsets
/// 	[(contig_id, contig_location), ...] // actual contig_ids, locations
///
/// the way this works is that the offsets array maps directly to the kmer_hash array index and
/// offsets can be used to extract the relevant (contig_id, contig_location) tuples.
///
/// Assume we have kmer_hashes 20, 50 and 99, where:
/// - hash 20 is found in (1, 10), (2, 3)
/// - hash 50 is found in (1, 2), (5, 50), (3, 20)
/// - hash 99 is found in (1, 25).
///
/// the hash vec will be [20, 50, 99].
/// the offset vec will be [0, 2, 5, 6].
/// the entry vec will be [(1, 10), (2, 3), (1, 2), (5, 50), (3, 20), (1, 25)]
///
/// The offset vec tells us that
/// 	kmer hash 20 starts at index 0 in the entry vec,
/// 	kmer hash 50 starts at index 2 in the entry vec,
/// 	kmer hash 99 starts at index 5 in the entry vec.
/// the last value in the offset vec in len(entry_vec).
///
/// Pretend have a read that generates a kmer hash of 20. We search (binary search) in our hash vec and
/// find its index to be 0. We check the offset vec for offset[i] and offset[i+1] = offset[0] & offset[1] = 0, 2.
/// We extract the relevant contig information as entry_vec[offset[i]..offset[i+1]] = entry_vec[0..2] = [(1, 10), (2, 3)].
/// We now know that contigs 1 and 2 at location 10 and 3, respectively, contain this kmer hash.
///
pub fn run(args: Args) -> Result<(), AppError> {
    let fq_reader = bio_fastq_reader(Some(args.fastq))?;
    let mut fs_reader = needletail_reader(Some(args.fasta))?;

    // -- minimizer index
    let spinner = ProgressBar::new_spinner();
    spinner.enable_steady_tick(Duration::from_millis(200));
    spinner.set_style(ProgressStyle::with_template(
        "Loading reads and building index {spinner:.blue} [{elapsed_precise}]",
    )?);

    let mut minimizer_index: Vec<(u64, usize, u32)> = vec![];

    let mut contig_index = 0;

    // We might need to use a HashMap or tuple if we want both contig name and contig length.
    let mut contig_lengths: Vec<usize> = vec![];

    while let Some(record) = fs_reader.next() {
        let record = match record {
            Ok(record) => record,
            Err(_) => continue,
        };

        contig_index += 1;
        contig_lengths.push(record.num_bases());

        // extract minimizers.
        let mut _simd_positions = vec![];

        // NOTE - using w=1 means extract ALL canonical kmers, which is what we want.
        for (pos, minimizer) in canonical_minimizers(args.kmer_size, 1)
            .run(record.seq().iter().as_slice(), &mut _simd_positions)
            .pos_and_values_u64()
        {
            minimizer_index.push((minimizer, contig_index, pos));
        }
    }
    spinner.finish_and_clear();

    // println!("{}", map.len());

    // // classify each contig.
    // //
    // // let mut writer = BufWriter::new(File::create(&args.outfile)?);

    // let spinner = ProgressBar::new_spinner();
    // spinner.enable_steady_tick(Duration::from_millis(200));
    // spinner.set_style(ProgressStyle::with_template(
    //     "Classifying query sequences {spinner:.blue} [{elapsed_precise}]",
    // )?);
    // // classify(
    // //     &reverse_index,
    // //     &valid_records,
    // //     query_reader,
    // //     &mut writer,
    // //     &*sketcher,
    // //     args.num_hits,
    // //     args.min_score,
    // // )?;

    // spinner.finish();

    Ok(())
}
