//! A variant annotater.

#![warn(missing_docs)]

/* std use */

/* crate use */

use anyhow::Context as _;
use clap::Parser as _;

/* project use */
use variant_myth::annotations_db;
use variant_myth::cli;
use variant_myth::error;
use variant_myth::sequences_db;
use variant_myth::variant;

fn main() -> error::Result<()> {
    // parse cli
    let params = cli::Command::parse();

    // Setup logger
    stderrlog::new()
        .module(module_path!())
        .quiet(params.quiet())
        .verbosity(params.verbosity())
        .timestamp(params.timestamp())
        .init()
        .context("stderrlog already create a logger")?;

    log::info!("Start read annotations");
    let annotations = annotations_db::AnnotationsDataBase::from_reader(params.annotations()?)?;
    log::info!("End read annotations");

    log::info!("Start read genome reference");
    let sequences = sequences_db::SequencesDataBase::from_reader(params.reference()?)?;
    log::info!("End read genome reference");

    log::info!("Start read vcf");
    let vcf_reader = variant::VcfReader::from_reader(params.variant()?);
    let mut annotations_counter = ahash::AHashMap::new();
    for result in vcf_reader {
        let variant = result?;
        let (seqname, interval) = variant.get_interval();
        for annot in annotations.get_annotation(seqname, interval.clone()) {
            let (s, i) = annot.get_interval();
            std::hint::black_box(sequences.get_interval(s, &i));
        }
        annotations_counter
            .entry(annotations.get_annotation(seqname, interval).len() as u64)
            .and_modify(|c| *c += 1)
            .or_insert(1u64);
    }
    log::info!("End read vcf");

    let mut res: Vec<(u64, u64)> = annotations_counter.iter().map(|x| (*x.0, *x.1)).collect();
    res.sort();
    let total: u64 = res.iter().map(|x| x.1).sum();
    for (key, value) in res {
        println!(
            "{},{},{:.4}",
            key,
            value,
            value as f64 / total as f64 * 100.0
        );
    }

    Ok(())
}
