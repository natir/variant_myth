use crate::error::Result;
use crate::output::writer::MythWriter;
use serde_json::Serializer;

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
    fn write_myth(&mut self, myth: crate::myth::Myth) -> crate::error::Result<()> {
        Ok(())
    }
    fn end_batch(&mut self) -> crate::error::Result<()> {
        Ok(())
    }
    fn close(self) -> crate::error::Result<()> {
        Ok(())
    }
}
