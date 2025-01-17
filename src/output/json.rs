//! The JSON writer module. Provides implementation for writing [`Myth`] objects to JSON.

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::myth;
use crate::output;

/// Struct to write Myth in json format
pub struct JsonWriter<W> {
    output_stream: W,
    json_format: JsonFormat,
}

/// Choose output JSON format
#[derive(Debug, Clone, Default, PartialEq, clap::ValueEnum)]
pub enum JsonFormat {
    /// Standard JSON, comma-delimited, with opening and closing brackets
    #[default]
    #[clap(name = "json")]
    Json,

    /// Newline-delimited JSON
    #[clap(name = "nd-json")]
    NdJson,
}

fn get_metadata_pretty() -> Vec<String> {
    crate::output::get_metadata()
        .iter()
        .filter_map(|(k, v)| serde_json::to_string_pretty(&(k, v)).ok())
        .collect()
}

fn get_metadata() -> Vec<String> {
    crate::output::get_metadata()
        .iter()
        .filter_map(|(k, v)| serde_json::to_string(&(k, v)).ok())
        .collect()
}

impl<W: std::io::Write> JsonWriter<W> {
    /// Create a new JsonWriter
    pub fn new(mut output_stream: W, json_format: JsonFormat) -> error::Result<Self> {
        if json_format == JsonFormat::Json {
            output_stream.write(b"{\n")?;
            let metadata_repr = serde_json::to_string_pretty(&get_metadata_pretty())?;
            output_stream.write(metadata_repr.as_bytes())?;
        } else {
            let metadata_repr = serde_json::to_string(&get_metadata())?;
            output_stream.write(metadata_repr.as_bytes())?;
        }
        Ok(Self {
            output_stream,
            json_format,
        })
    }
}

/// TODO
/// Work in progress
impl<W: std::io::Write> output::MythWriter for JsonWriter<W> {
    fn add_myth(&mut self, myth: myth::Myth) -> error::Result<()> {
        match self.json_format {
            JsonFormat::Json => {
                let variant_repr = serde_json::to_string_pretty(&myth.variant)?;
                let one_line_per_record = true;
                if one_line_per_record {
                    self.output_stream.write(b"")?;
                }
            }
            JsonFormat::NdJson => {}
        }
        Ok(())
    }
    fn batch_full(&self) -> bool {
        false
    }
    fn finalize(&mut self) -> error::Result<()> {
        if self.json_format == JsonFormat::Json {
            self.output_stream.write(b"\n}")?;
        }
        self.output_stream.flush()?;
        Ok(())
    }
    fn write_batch(&mut self) -> error::Result<()> {
        Ok(())
    }
}
