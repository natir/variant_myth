//! A variant annotater.

#![warn(missing_docs)]

/* std use */

/* crate use */

use anyhow::Context as _;
use clap::Parser as _;

use variant_myth::annotate;
/* project use */
use variant_myth::annotations_db;
use variant_myth::cli;
use variant_myth::error;
use variant_myth::output;
use variant_myth::sequences_db;
use variant_myth::translate;
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

    #[cfg(feature = "parallel")]
    rayon::ThreadPoolBuilder::new()
        .num_threads(params.threads())
        .build_global()?;

    let (annotations, sequences, translate) = get_database(&params)?;

    log::info!("Start annotate variant");
    let vcf_reader = variant::VcfReader::from_reader(params.variant()?);

    let output_type = crate::output::OutputFileType::Parquet;

    let block_size = 1 << 13;
    match output_type {
        output::OutputFileType::Parquet => {
            let writer = output::parquet::ParquetWriter::new(params.output()?, block_size)?;
            annotate(
                &annotations,
                &sequences,
                &translate,
                vcf_reader,
                params.no_annotation(),
                block_size,
                writer,
            )?;
        }
        output::OutputFileType::JSON => {
            let writer = output::json::JsonWriter::new(params.output()?)?;
            annotate(
                &annotations,
                &sequences,
                &translate,
                vcf_reader,
                params.no_annotation(),
                block_size,
                writer,
            )?;
        }
    }
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
    translate::Translate,
)> {
    log::info!("Start read genome reference");
    let sequences = sequences_db::SequencesDataBase::from_reader(params.reference()?)?;
    log::info!("End read genome reference");

    log::info!("Start read annotations");
    let annotations = annotations_db::AnnotationsDataBase::from_reader(
        params.annotations()?,
        params.updown_distance(),
    )?;
    log::info!("End read annotations");

    log::info!("Start read translation table");
    let translate = if let Some(reader) = params.translate()? {
        translate::Translate::from_reader(reader)?
    } else {
        translate::Translate::default()
    };
    log::info!("End read translation table");

    Ok((annotations, sequences, translate))
}

#[cfg(feature = "parallel")]
#[inline(always)]
fn get_database(
    params: &cli::Command,
) -> error::Result<(
    annotations_db::AnnotationsDataBase,
    sequences_db::SequencesDataBase,
    translate::Translate,
)> {
    let seq_reader = params.reference()?;
    let seq_thread = std::thread::spawn(|| {
        log::info!("Start read genome reference");
        let sequences = sequences_db::SequencesDataBase::from_reader(seq_reader)?;
        log::info!("End read genome reference");

        Ok::<sequences_db::SequencesDataBase, anyhow::Error>(sequences)
    });

    let annot_reader = params.annotations()?;
    let updown_distance = params.updown_distance();
    let annot_thread = std::thread::spawn(move || {
        log::info!("Start read annotations");
        let annotations =
            annotations_db::AnnotationsDataBase::from_reader(annot_reader, updown_distance)?;

        log::info!("End read annotations");

        Ok::<annotations_db::AnnotationsDataBase, anyhow::Error>(annotations)
    });

    let translate_reader = params.translate()?;
    let translate_thread = std::thread::spawn(|| {
        log::info!("Start read translation table");
        let translate: Result<translate::Translate, anyhow::Error> =
            if let Some(reader) = translate_reader {
                translate::Translate::from_reader(reader)
            } else {
                Ok(translate::Translate::default())
            };
        log::info!("End read translation table");

        translate
    });

    Ok((
        annot_thread.join().unwrap()?,
        seq_thread.join().unwrap()?,
        translate_thread.join().unwrap()?,
    ))
}
