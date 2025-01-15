//! Manage output format

/* std use */

/* crate use */

/* module declaration */
mod json;
mod parquet;

/* project use */
use crate::error;
use crate::myth;

/* reexport */
pub use json::JsonWriter;
pub use parquet::ParquetWriter;

/// Common metadata to all output
pub fn get_metadata() -> Vec<(&'static str, &'static str)> {
    vec![
        ("chr", "chromosome name same ase original vcf"),
        ("pos", "position of variant"),
        ("ref", "reference sequence"),
        ("alt", "alternative sequence"),
        ("source", "source of variant in gff3 file"),
        (
            "feature",
            "type of feature affected by variant gene/transcript",
        ),
        (
            "name",
            "name of feature, same value of Name gff3 attributes",
        ),
        ("id", "id of feature, same value of Id gff3 attributes"),
        ("effect", "List of sequence ontology terms"),
        (
            "impact",
            "0: UNKOWN, 1:LOW, 2:MODIFIER, 3: MODERATE, 4:HIGH",
        ),
    ]
}

/// A trait to implement your own writers
/// You only need to implement [`add_myth`], [`batch_full`], [`write_batch`], and [`close`]
/// [`add_myth`]: ./fn.add_myth.html
/// [`batch_full`]: ./fn.batch_full.html
/// [`write_batch`]: ./fn.write_batch.html
/// [`close`]: ./fn.close.html
pub trait MythWriter {
    /// This method is called for each variant for which a Myth object was found.
    /// Do not implement this method!
    fn write_myth(&mut self, myth: myth::Myth) -> error::Result<()> {
        self.add_myth(myth)?;
        if self.batch_full() {
            self.write_batch()?;
        }

        Ok(())
    }

    /// This method is called when we are done with all the input variants.
    /// Do not implement this method!
    fn close(&mut self) -> error::Result<()> {
        self.write_batch()?;
        self.finalize()
    }

    /// Implement this method to add a Myth object to your batch
    fn add_myth(&mut self, myth: myth::Myth) -> error::Result<()>;

    /// Return true once your batch is full. This will trigger [`write_batch`] at the next call to [`write_myth`].
    /// [`write_batch`]: ./fn.write_batch.html
    /// [`write_myth`]: ./fn.write_myth.html
    fn batch_full(&self) -> bool;

    /// Specialized method to write a whole batch
    fn write_batch(&mut self) -> error::Result<()>;

    /// Specialized method to safely close your writer
    fn finalize(&mut self) -> error::Result<()>;
}
