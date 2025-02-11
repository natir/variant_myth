//! Command Line Interface declaration of project variant_myth

/* std use */

/* crate use */
use std::io::Read;

use enumflags2::BitFlag as _;

/* project use */
use crate::error;
use crate::output;
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
    /// Variants path
    #[clap(short = 'i', long = "input", required = true)]
    variant_paths: Vec<std::path::PathBuf>,

    /// Reference genome path
    #[clap(short = 'r', long = "reference")]
    reference_path: std::path::PathBuf,

    /// Annotation path
    #[clap(short = 'a', long = "annotations", required = true)]
    annotations_path: Vec<std::path::PathBuf>,

    /// Translate table path, if not set use human
    #[clap(short = 't', long = "translate")]
    translate_path: Option<std::path::PathBuf>,

    /// [Up|Down]stream transcript distance, default: 5,000
    #[clap(short = 'd', long = "updown-distance")]
    updown_distance: Option<u64>,

    /// Select which type of annotation you want run
    #[clap(short = 'c', long = "annotators-choices")]
    annotators_choices: Vec<variant2myth::AnnotatorsChoicesRaw>,

    /// Output subcommand
    #[clap(subcommand)]
    pub output: OutputSubCommand,

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
    ) -> error::Result<Vec<std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>>> {
        self.variant_paths
            .iter()
            .map(|p| get_reader(p).map(std::io::BufReader::new))
            .collect()
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

/// Subcommand to control how output are write
#[derive(clap::Subcommand, std::fmt::Debug)]
pub enum OutputSubCommand {
    /// Output are write in parquet format
    #[cfg(feature = "parquet")]
    Parquet(Parquet),
    /// Output are write in json format
    #[cfg(feature = "json")]
    Json(Json),
}

impl OutputSubCommand {
    /// Create myth writer
    pub fn writers(&self) -> error::Result<Vec<Box<dyn output::MythWriter + std::marker::Send>>> {
        match self {
            #[cfg(feature = "parquet")]
            OutputSubCommand::Parquet(obj) => obj.writers(),
            #[cfg(feature = "json")]
            OutputSubCommand::Json(obj) => obj.writers(),
        }
    }
}

/// Output are write in parquet format
#[derive(clap::Args, std::fmt::Debug)]
#[cfg(feature = "parquet")]
pub struct Parquet {
    /// Output path
    #[clap(short = 'p', long = "path", required = true)]
    paths: Vec<std::path::PathBuf>,

    /// Size of parquet block, default value 65536
    #[clap(short = 'b', long = "block-size")]
    block_size: Option<usize>,
}

#[cfg(feature = "parquet")]
impl Parquet {
    /// Create myth writer
    pub fn writers(&self) -> error::Result<Vec<Box<dyn output::MythWriter + std::marker::Send>>> {
        let mut result = Vec::new();

        for p in &self.paths {
            result.push(Box::new(output::ParquetWriter::new(
                std::fs::File::create(p).map(std::io::BufWriter::new)?,
                self.block_size(),
            )?)
                as Box<dyn output::MythWriter + std::marker::Send>);
        }

        Ok(result)
    }

    /// Get block_size
    pub fn block_size(&self) -> usize {
        self.block_size.unwrap_or(1 << 16)
    }
}

/// Output are write in json format
#[derive(clap::Args, std::fmt::Debug)]
#[cfg(feature = "json")]
pub struct Json {
    /// Output path
    #[clap(short = 'p', long = "path", required = true)]
    paths: Vec<std::path::PathBuf>,

    /// Json format, default json
    #[clap(short = 'f', long = "format")]
    json_format: Option<output::JsonFormat>,
}

#[cfg(feature = "json")]
impl Json {
    /// Create myth writer
    pub fn writers(&self) -> error::Result<Vec<Box<dyn output::MythWriter + std::marker::Send>>> {
        let mut result = Vec::new();

        for p in &self.paths {
            result.push(Box::new(output::JsonWriter::new(
                std::fs::File::create(p).map(std::io::BufWriter::new)?,
                self.format(),
            )?)
                as Box<dyn output::MythWriter + std::marker::Send>);
        }

        Ok(result)
    }

    /// Get format value
    pub fn format(&self) -> output::JsonFormat {
        self.json_format.unwrap_or_default()
    }
}
