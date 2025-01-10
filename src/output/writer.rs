use crate::myth::Myth;

use crate::error::Result;

/// A trait to implement your own writers
/// You only need to implement [`add_myth`], [`batch_full`], [`write_batch`], and [`close`]
/// [`add_myth`]: ./fn.add_myth.html
/// [`batch_full`]: ./fn.batch_full.html
/// [`write_batch`]: ./fn.write_batch.html
/// [`close`]: ./fn.close.html
pub trait MythWriter {
    /// This method is called for each variant for which a Myth object was found.
    /// Do not implement this method!
    fn write_myth(&mut self, myth: Myth) -> Result<()> {
        self.add_myth(myth)?;
        if self.batch_full() {
            self.write_batch()?;
        }

        Ok(())
    }

    /// This method is called when we are done with all the input variants.
    /// Do not implement this method!
    fn close(&mut self) -> Result<()> {
        self.write_batch()?;
        self.finalize()
    }

    /// Implement this method to add a Myth object to your batch
    fn add_myth(&mut self, myth: Myth) -> Result<()>;

    /// Return true once your batch is full. This will trigger [`write_batch`] at the next call to [`write_myth`].
    /// [`write_batch`]: ./fn.write_batch.html
    /// [`write_myth`]: ./fn.write_myth.html
    fn batch_full(&self) -> bool;

    /// Specialized method to write a whole batch
    fn write_batch(&mut self) -> Result<()>;

    /// Specialized method to safely close your writer
    fn finalize(&mut self) -> Result<()>;
}
