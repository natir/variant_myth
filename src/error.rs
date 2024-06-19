//! Error struct of project variant_myth

/* crate use */
use anyhow;
use thiserror;

/// Enum to manage error
#[derive(std::fmt::Debug, thiserror::Error)]
pub enum Error {
    /// std Input Output error
    #[error(transparent)]
    StdIO(#[from] std::io::Error),

    /// Error in logging system configuration
    #[error(transparent)]
    Log(#[from] log::SetLoggerError),

    /// Error in int parsing
    #[error(transparent)]
    IntParsing(#[from] std::num::ParseIntError),

    /// Error in float parsing
    #[error(transparent)]
    FloatParsing(#[from] std::num::ParseFloatError),

    /// Not a valid strand
    #[error("Gff contain an invalid strand")]
    GffBadStrand,

    /// Not a valid Frame
    #[error("Gff contain an invalid frame")]
    GffBadFrame,

    /// Bad vcf record
    #[error("Bad vcf record")]
    VcfBadRecord,
}

/// Alias of result
pub type Result<T> = anyhow::Result<T>;
