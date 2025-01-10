//! Command Line Interface declaration of project variant_myth

/* std use */

/* crate use */
use std::io::Read;

use enumflags2::BitFlag as _;

/* project use */
use crate::error;
use crate::variant2myth;

fn get_reader(
    path: &std::path::PathBuf,
) -> error::Result<Box<dyn std::io::Read + std::marker::Send>> {
    let file = std::fs::File::open(path)?;
    let boxed = Box::new(file);
    let (reader, _compression) = niffler::send::get_reader(boxed)?;

    Ok(reader)
}

/// A variant annotater.
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(
    name = "variant_myth",
    bin_name = "variant_myth",
    version = "0.1",
    author = "Pierre Marijon <pierre@marijon.fr>"
)]
#[command(propagate_version = true)]
pub struct Command {
    // Specific option
    /// Variant path
    #[clap(short = 'i', long = "input")]
    variant_path: std::path::PathBuf,

    /// Reference genome path
    #[clap(short = 'r', long = "reference")]
    reference_path: std::path::PathBuf,

    /// Annotation path
    #[clap(short = 'a', long = "annotations", required = true)]
    annotations_path: Vec<std::path::PathBuf>,

    /// Translate table path, if not set use human
    #[clap(short = 't', long = "translate")]
    translate_path: Option<std::path::PathBuf>,

    /// Output path
    #[clap(short = 'o', long = "output")]
    output_path: std::path::PathBuf,

    /// [Up|Down]stream transcript distance, default: 5,000
    #[clap(short = 'd', long = "updown-distance")]
    updown_distance: Option<u64>,

    /// Select which type of annotation you want run
    #[clap(short = 'c', long = "annotators-choices")]
    annotators_choices: Vec<variant2myth::AnnotatorsChoicesRaw>,

    // Generic option
    #[cfg(feature = "parallel")]
    /// Number of theard use 0 use all avaible core, default value 0
    #[clap(long = "threads")]
    threads: Option<usize>,

    /// Silence all output
    #[clap(short = 'q', long = "quiet")]
    quiet: bool,

    /// Verbose mode (-v, -vv, -vvv, etc)
    #[clap(short = 'v', long = "verbosity", action = clap::ArgAction::Count)]
    verbosity: u8,

    /// Timestamp (sec, ms, ns, none)
    #[clap(short = 'T', long = "timestamp")]
    ts: Option<stderrlog::Timestamp>,
}

impl Command {
    /// Get variant reader
    pub fn variant(
        &self,
    ) -> error::Result<std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>> {
        get_reader(&self.variant_path).map(std::io::BufReader::new)
    }

    /// Get reference reader
    pub fn reference(
        &self,
    ) -> error::Result<std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>> {
        get_reader(&self.reference_path).map(std::io::BufReader::new)
    }

    /// Get annotations reader
    pub fn annotations(
        &self,
    ) -> error::Result<std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>> {
        let mut handle: Box<dyn std::io::Read + std::marker::Send> =
            Box::new(std::io::Cursor::new(vec![]));

        for path in &self.annotations_path {
            handle = Box::new(handle.chain(get_reader(path)?));
        }

        Ok(std::io::BufReader::new(handle))
    }

    /// Get translate reader
    pub fn translate(
        &self,
    ) -> error::Result<Option<std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>>> {
        if let Some(path) = &self.translate_path {
            get_reader(path).map(std::io::BufReader::new).map(Some)
        } else {
            Ok(None)
        }
    }

    /// Get output writer
    pub fn output(
        &self,
    ) -> error::Result<std::io::BufWriter<Box<dyn crate::WriteSeek + std::marker::Send>>> {
        let file = std::fs::File::create(&self.output_path)?;
        let boxed = Box::new(file);

        Ok(std::io::BufWriter::new(boxed))
    }

    /// Get [Up|Down]stream transcript distance
    pub fn updown_distance(&self) -> u64 {
        self.updown_distance.unwrap_or(5000)
    }

    /// Get no annotation
    pub fn annotators_choices(&self) -> variant2myth::AnnotatorsChoices {
        self.annotators_choices
            .iter()
            .fold(variant2myth::AnnotatorsChoicesRaw::empty(), |acc, x| {
                acc | *x
            })
    }

    /// Get number of thread
    #[cfg(feature = "parallel")]
    pub fn threads(&self) -> usize {
        self.threads.unwrap_or(0)
    }

    /// Get verbosity level
    pub fn verbosity(&self) -> usize {
        self.verbosity as usize
    }

    /// Get quiet
    pub fn quiet(&self) -> bool {
        self.quiet
    }

    /// Get timestamp granularity
    pub fn timestamp(&self) -> stderrlog::Timestamp {
        self.ts.unwrap_or(stderrlog::Timestamp::Off)
    }
}
