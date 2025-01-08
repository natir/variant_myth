//! The JSON writer module. Provides implementation for writing [`Myth`] objects to JSON.
use crate::error::Result;
use crate::output::writer::MythWriter;
use serde_json::Serializer;

/// This is the JSON Writer struct
pub struct JsonWriter<W> {
    serializer: Serializer<W>,
}

impl<W: std::io::Write> JsonWriter<W> {
    pub fn new(writer: W) -> Result<Self> {
        Ok(Self {
            serializer: serde_json::Serializer::new(writer),
        })
    }
}

impl<W: std::io::Write> MythWriter for JsonWriter<W> {
    fn add_myth(&mut self, myth: crate::myth::Myth) -> Result<()> {
        Ok(())
    }
    fn batch_full(&self) -> bool {
        true
    }
    fn close(&mut self) -> Result<()> {
        Ok(())
    }
    fn write_batch(&mut self) -> Result<()> {
        Ok(())
    }
}
