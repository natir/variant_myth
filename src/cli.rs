//! Command Line Interface declaration of project variant_myth

/* std use */

/* crate use */

use std::io::Read;

/* project use */
use crate::error;

/// A variant annotater.
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(
    name = "variant_myth",
    version = "0.1",
    author = "Pierre Marijon <pierre@marijon.fr>"
)]
pub struct Command {
    // Specific option
    /// Variant path
    #[clap(short = 'i', long = "input")]
    variant_path: std::path::PathBuf,

    /// Reference genome path
    #[clap(short = 'r', long = "reference")]
    reference_path: std::path::PathBuf,

    /// Annotation path
    #[clap(short = 'a', long = "annotations")]
    annotations_path: Vec<std::path::PathBuf>,

    /// Output path
    #[clap(short = 'o', long = "output")]
    output_path: std::path::PathBuf,

    // Generic option
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
    fn get_reader(path: &std::path::PathBuf) -> error::Result<Box<dyn std::io::Read>> {
        let file = std::fs::File::open(path)?;
        let boxed = Box::new(file);
        let (reader, _compression) = niffler::send::get_reader(boxed)?;

        Ok(reader)
    }

    /// Get variant reader
    pub fn variant(&self) -> error::Result<std::io::BufReader<Box<dyn std::io::Read>>> {
        Command::get_reader(&self.variant_path).map(std::io::BufReader::new)
    }

    /// Get reference reader
    pub fn reference(&self) -> error::Result<std::io::BufReader<Box<dyn std::io::Read>>> {
        Command::get_reader(&self.variant_path).map(std::io::BufReader::new)
    }

    /// Get annotations reader
    pub fn annotations(&self) -> error::Result<std::io::BufReader<Box<dyn std::io::Read>>> {
        let mut handle: Box<dyn std::io::Read> = Box::new(std::io::Cursor::new(vec![]));

        for path in &self.annotations_path {
            handle = Box::new(handle.chain(Command::get_reader(path)?));
        }

        Ok(std::io::BufReader::new(handle))
    }

    /// Get output writer
    pub fn output(&self) -> error::Result<std::io::BufWriter<Box<dyn std::io::Write>>> {
        let file = std::fs::File::create(&self.output_path)?;
        let boxed = Box::new(file);

        Ok(std::io::BufWriter::new(boxed))
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
