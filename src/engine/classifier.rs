use crate::errors::AppError;
use bio::io::fasta::{Reader, Record};
use bio_utils_rs::simd_sketch::Sketcher;
use dashmap::DashMap;
use fixedbitset::FixedBitSet;
use rayon::prelude::*;
use rustc_hash::FxBuildHasher;
use std::collections::HashSet;
use std::io::Write;
use std::{
    fs::File,
    io::{BufReader, BufWriter},
};

struct QueryResult {
    query_id: String,
    hits: Vec<HitResult>,
}

struct HitResult {
    db_id: String,
    shared_hashes: usize,
    score: f64,
}

pub fn classify(
    reverse_index: &DashMap<u64, FixedBitSet, FxBuildHasher>,
    valid_records: &[Record],
    query_reader: Reader<BufReader<File>>,
    writer: &mut BufWriter<File>,
    sketcher: &dyn Sketcher,
    num_hits: usize,
    min_score: f64,
) -> Result<(), AppError> {
    // Pre-collect query records so we can use par_iter (more efficient than par_bridge)
    let query_records: Vec<Record> = query_reader.records().filter_map(|r| r.ok()).collect();

    // Process queries in parallel, collect results (no mutex needed)
    let results: Vec<QueryResult> = query_records
        .par_iter()
        .filter_map(|r| {
            let query_hashes: HashSet<u64> = sketcher.sketch(r.seq());
            let num_query_hashes = query_hashes.len();
            if num_query_hashes == 0 {
                return None;
            }

            let mut hits: Vec<usize> = vec![0; valid_records.len()];

            for hash in &query_hashes {
                // Use .get() for a read lock instead of .entry().and_modify() which takes a write lock
                if let Some(bitset) = reverse_index.get(hash) {
                    for idx in bitset.ones() {
                        hits[idx] += 1;
                    }
                }
            }

            // Extract top N hits, filtering zeros and scores below min allowed score.
            let mut scored: Vec<(usize, usize)> = hits
                .into_iter()
                .enumerate()
                .filter(|(_, count)| {
                    *count > 0 && *count as f64 / num_query_hashes as f64 >= min_score
                })
                .map(|(idx, count)| (count, idx))
                .collect();

            if scored.is_empty() {
                return None;
            }

            let n = num_hits.min(scored.len());
            scored.select_nth_unstable_by(n - 1, |a, b| b.0.cmp(&a.0));
            scored.truncate(n);
            scored.sort_unstable_by(|(count_a, _), (count_b, _)| count_b.cmp(count_a));

            let hit_results: Vec<HitResult> = scored
                .into_iter()
                .map(|(count, idx)| HitResult {
                    db_id: valid_records[idx].id().to_string(),
                    shared_hashes: count,
                    score: count as f64 / num_query_hashes as f64,
                })
                .collect();

            Some(QueryResult {
                query_id: r.id().to_string(),
                hits: hit_results,
            })
        })
        .collect();

    // Write TSV header
    writeln!(writer, "query_id\tsubject_id\tshared_hashes\tscore")?;

    // Write all results sequentially to avoid mutex contention
    for result in &results {
        for hit in &result.hits {
            writeln!(
                writer,
                "{}\t{}\t{}\t{:.6}",
                result.query_id, hit.db_id, hit.shared_hashes, hit.score
            )?;
        }
    }
    writer.flush()?;

    Ok(())
}
