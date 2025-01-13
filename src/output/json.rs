//! The JSON writer module. Provides implementation for writing [`Myth`] objects to JSON.

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::myth;
use crate::output;

/// Struct to write Myth in json format
pub struct JsonWriter<W> {
    _serializer: serde_json::Serializer<W>,
}

impl<W: std::io::Write> JsonWriter<W> {
    /// Create a new JsonWriter
    pub fn new(writer: W) -> error::Result<Self> {
        Ok(Self {
            _serializer: serde_json::Serializer::new(writer),
        })
    }
}

impl<W: std::io::Write> output::MythWriter for JsonWriter<W> {
    fn add_myth(&mut self, _myth: myth::Myth) -> error::Result<()> {
        Ok(())
    }
    fn batch_full(&self) -> bool {
        true
    }
    fn finalize(&mut self) -> error::Result<()> {
        Ok(())
    }
    fn write_batch(&mut self) -> error::Result<()> {
        Ok(())
    }
}
