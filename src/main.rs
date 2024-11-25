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

    match params.subcommands {
        cli::SubCommands::Var2Gene(ref subparams) => main_var2gene(&params, subparams),
        cli::SubCommands::Var2Full(ref subparams) => main_var2full(&params, subparams),
    }
}

fn main_var2gene(params: &cli::Command, subparams: &cli::Var2Gene) -> error::Result<()> {
    log::info!("Start read annotations");
    let annotations = get_annotations(subparams.annotations()?, params.updown_distance())?;
    log::info!("End read annotations");

    log::info!("Start annotate variant");
    let vcf_reader = variant::VcfReader::from_reader(subparams.variant()?);
    variant_myth::vcf2gene(&annotations, vcf_reader, params.output()?)?;
    log::info!("End annotate variant");

    Ok(())
}

fn main_var2full(params: &cli::Command, subparams: &cli::Var2Full) -> error::Result<()> {
    let (annotations, sequences, translate) = get_database(params, subparams)?;

    log::info!("Start annotate variant");
    let vcf_reader = variant::VcfReader::from_reader(subparams.variant()?);
    variant_myth::vcf2json(
        &annotations,
        &sequences,
        &translate,
        vcf_reader,
        params.output()?,
    )?;
    log::info!("End annotate variant");

    Ok(())
}

fn get_annotations(
    reader: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>,
    updown_distance: u64,
) -> error::Result<annotations_db::AnnotationsDataBase> {
    let annotations = annotations_db::AnnotationsDataBase::from_reader(reader, updown_distance)?;

    Ok::<annotations_db::AnnotationsDataBase, anyhow::Error>(annotations)
}

#[cfg(not(feature = "parallel"))]
#[inline(always)]
fn get_database(
    params: &cli::Command,
    subparams: &cli::Var2Full,
) -> error::Result<(
    annotations_db::AnnotationsDataBase,
    sequences_db::SequencesDataBase,
    translate::Translate,
)> {
    log::info!("Start read genome reference");
    let sequences = sequences_db::SequencesDataBase::from_reader(subparams.reference()?)?;
    log::info!("End read genome reference");

    log::info!("Start read annotations");
    let annotations = get_annotations(subparams.annotations()?, params.updown_distance())?;
    log::info!("End read annotations");

    log::info!("Start read translation table");
    let translate = if let Some(reader) = subparams.translate()? {
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
    subparams: &cli::Var2Full,
) -> error::Result<(
    annotations_db::AnnotationsDataBase,
    sequences_db::SequencesDataBase,
    translate::Translate,
)> {
    let seq_reader = subparams.reference()?;
    let seq_thread = std::thread::spawn(|| {
        log::info!("Start read genome reference");
        let sequences = sequences_db::SequencesDataBase::from_reader(seq_reader)?;
        log::info!("End read genome reference");

        Ok::<sequences_db::SequencesDataBase, anyhow::Error>(sequences)
    });

    let annot_reader = subparams.annotations()?;
    let updown_distance = params.updown_distance();
    let annot_thread = std::thread::spawn(move || {
        log::info!("Start read annotations");
        let annotations = get_annotations(annot_reader, updown_distance)?;
        log::info!("End read annotations");

        Ok::<annotations_db::AnnotationsDataBase, anyhow::Error>(annotations)
    });

    let translate_reader = subparams.translate()?;
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
