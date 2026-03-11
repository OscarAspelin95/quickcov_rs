use crate::engine::index::build_minimizer_index;
use crate::{args::Args, errors::AppError};
use bio_utils_rs::io::{bio_fastq_reader, get_bufwriter};
use packed_seq::{PackedSeqVec, SeqVec};
use rayon::prelude::*;
use simd_minimizers::canonical_minimizers;

pub fn run(args: Args) -> Result<(), AppError> {
    let fq_reader = bio_fastq_reader(Some(args.fastq))?;

    let (minimizer_index, contig_index) = build_minimizer_index(args.fasta, args.kmer_size)?;

    fq_reader.records().par_bridge().for_each(|record| {
        let record = match record {
            Ok(record) => record,
            Err(_) => return,
        };

        let capacity = record.seq().len() * 2 / (1 + 1);
        let packed_seq = PackedSeqVec::from_ascii(record.seq());
        let mut minimizer_positions = Vec::with_capacity(capacity);

        // extract minimizers
        for minimizer in canonical_minimizers(args.kmer_size, 1)
            .run(packed_seq.as_slice(), &mut minimizer_positions)
            .values_u64()
        {
            if let Some(contig_info) = minimizer_index.get(&minimizer) {
                contig_info.iter().for_each(|(id, pos)| {
                    contig_index
                        .entry(*id)
                        .and_modify(|positions| positions.1[*pos] += 1);
                });
            }
        }
    });

    //
    let mut writer = get_bufwriter(Some(args.outfile))?;

    contig_index.iter().for_each(|entry| {
        let (contig_name, contig_positions) = entry.value();

        let mut total_coverage: usize = 0;
        let mut num_valid_positions: usize = 0;

        for coverage in contig_positions {
            total_coverage += coverage;

            if *coverage >= args.min_kmer_coverage {
                num_valid_positions += 1
            }
        }

        let mean_coverage: f64 = total_coverage as f64 / contig_positions.len() as f64;
        let frac_coverage: f64 = num_valid_positions as f64 / contig_positions.len() as f64;

        writer
            .write_all(format!("{contig_name}\t{mean_coverage}\t{frac_coverage}\n").as_bytes())
            .expect("");
    });

    Ok(())
}
