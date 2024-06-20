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

    let (annotations, sequences) = get_database(&params)?;

    log::info!("Start annotate variant");
    let vcf_reader = variant::VcfReader::from_reader(params.variant()?);
    variant_myth::vcf2json(&annotations, &sequences, vcf_reader, params.output()?)?;
    log::info!("End annotate variant");

    Ok(())
}

#[cfg(not(feature = "parallel"))]
#[inline(always)]
fn get_database(
    params: &cli::Command,
) -> error::Result<(
    annotations_db::AnnotationsDataBase,
    sequences_db::SequencesDataBase,
)> {
    log::info!("Start read annotations");
    let annotations = annotations_db::AnnotationsDataBase::from_reader(params.annotations()?)?;
    log::info!("End read annotations");

    log::info!("Start read genome reference");
    let sequences = sequences_db::SequencesDataBase::from_reader(params.reference()?)?;
    log::info!("End read genome reference");

    Ok((annotations, sequences))
}

#[cfg(feature = "parallel")]
#[inline(always)]
fn get_database(
    params: &cli::Command,
) -> error::Result<(
    annotations_db::AnnotationsDataBase,
    sequences_db::SequencesDataBase,
)> {
    let annot_reader = params.annotations()?;
    let annot_thread = std::thread::spawn(|| {
        log::info!("Start read annotations");
        let annotations = annotations_db::AnnotationsDataBase::from_reader(annot_reader)?;
        log::info!("End read annotations");

        Ok::<annotations_db::AnnotationsDataBase, anyhow::Error>(annotations)
    });

    let seq_reader = params.reference()?;
    let seq_thread = std::thread::spawn(|| {
        log::info!("Start read genome reference");
        let sequences = sequences_db::SequencesDataBase::from_reader(seq_reader)?;
        log::info!("End read genome reference");

        Ok::<sequences_db::SequencesDataBase, anyhow::Error>(sequences)
    });

    Ok((annot_thread.join().unwrap()?, seq_thread.join().unwrap()?))
}
