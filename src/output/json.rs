//! The JSON writer module. Provides implementation for writing [`Myth`] objects to JSON.

/* std use */

/* crate use */

use std::collections::BTreeMap;

use serde_json::json;

/* project use */
use crate::error;
use crate::myth;
use crate::output;

/// Struct to write Myth in json format
pub struct JsonWriter<W> {
    output_stream: W,
    json_format: JsonFormat,
    write_state: WriteState,
}

enum WriteState {
    WroteMetadata,
    WroteRecord,
}

/// Choose output JSON format
#[derive(Debug, Clone, Default, PartialEq, clap::ValueEnum)]
pub enum JsonFormat {
    /// Standard JSON, comma-delimited, with opening and closing brackets
    #[default]
    Json,

    /// Newline-delimited JSON
    NdJson,
}

fn get_metadata() -> serde_json::Value {
    let mut map = BTreeMap::new();
    for (k, v) in crate::output::get_metadata() {
        map.entry(k).or_insert(v);
    }
    json!(map)
}

impl<W: std::io::Write> JsonWriter<W> {
    /// Create a new JsonWriter
    pub fn new(mut output_stream: W, json_format: JsonFormat) -> error::Result<Self> {
        let metadata_repr = json!(get_metadata());
        if json_format == JsonFormat::Json {
            output_stream.write_fmt(format_args!("{{\n\"metadata\":\n{:#}", metadata_repr))?;
        } else {
            output_stream.write_fmt(format_args!("{{\"metadata\":{}}}", metadata_repr))?;
        }
        if json_format == JsonFormat::Json {
            output_stream.write_all(b",\n\"variants\": [\n")?;
        }
        if json_format == JsonFormat::NdJson {
            output_stream.write_all(b"\n")?;
        }
        Ok(Self {
            output_stream,
            json_format,
            write_state: WriteState::WroteMetadata,
        })
    }
}

impl<W: std::io::Write> output::MythWriter for JsonWriter<W> {
    fn add_myth(&mut self, myth: myth::Myth) -> error::Result<()> {
        let myth_repr = json!(
            {
                "variant": &myth.variant,
                "myth": &myth.annotations
            }
        );
        let separator = match self.write_state {
            WriteState::WroteMetadata => "".to_string(),
            WriteState::WroteRecord if self.json_format == JsonFormat::Json => ",\n".to_string(),
            WriteState::WroteRecord if self.json_format == JsonFormat::NdJson => "\n".to_string(),
            _ => unreachable!(),
        };
        match self.json_format {
            JsonFormat::Json => {
                self.output_stream
                    .write_fmt(format_args!("{}{:#}", separator, &myth_repr))?;
            }
            JsonFormat::NdJson => {
                self.output_stream
                    .write_fmt(format_args!("{}{}", separator, &myth_repr))?;
            }
        }
        self.write_state = match self.write_state {
            WriteState::WroteMetadata => WriteState::WroteRecord,
            WriteState::WroteRecord => WriteState::WroteRecord,
        };
        Ok(())
    }
    fn batch_full(&self) -> bool {
        false
    }
    fn finalize(&mut self) -> error::Result<()> {
        if self.json_format == JsonFormat::Json {
            self.output_stream.write_all(b"]\n}")?;
        } else {
            // Do nothing
        }
        self.output_stream.flush()?;
        Ok(())
    }

    fn write_batch(&mut self) -> error::Result<()> {
        log::debug!("Calling write_batch on JsonWriter makes no sense!");
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::myth::Myth;

    fn get_two_myths() -> (Myth, Myth) {
        use crate::effect;
        use crate::myth::{AnnotationMyth, Myth};
        use crate::variant;

        let mut annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .feature(b"gene".to_vec())
            .name(b"gene1".to_vec())
            .id(b"1111".to_vec());

        annotation.add_effect(effect::Effect::GeneVariant);
        annotation.add_effect(effect::Effect::ExonRegion);

        let mut myth = Myth::from_variant(variant::Variant {
            seqname: b"93".to_vec(),
            position: 2036067340,
            ref_seq: b"T".to_vec(),
            alt_seq: b".".to_vec(),
        });
        myth.add_annotation(annotation.build().unwrap());

        let mut annotation2 = AnnotationMyth::builder()
            .source(b"test2".to_vec())
            .feature(b"gene2".to_vec())
            .name(b"gene51".to_vec())
            .id(b"7777".to_vec());

        annotation2.add_effect(effect::Effect::DisruptiveInframeDeletion);
        annotation2.add_effect(effect::Effect::ExonRegion);

        let mut myth2 = Myth::from_variant(variant::Variant {
            seqname: b"21".to_vec(),
            position: 1970,
            ref_seq: b"C".to_vec(),
            alt_seq: b"T".to_vec(),
        });

        myth2.add_annotation(annotation2.build().unwrap());
        (myth, myth2)
    }

    #[test]
    fn test_write_myth_json() {
        use super::*;

        use crate::output::MythWriter;

        let (myth, myth2) = get_two_myths();

        let output_stream: Vec<u8> = Vec::new();
        let mut annot_writer = JsonWriter::new(output_stream, JsonFormat::Json).unwrap();

        annot_writer.write_myth(myth).unwrap();
        annot_writer.write_myth(myth2).unwrap();
        annot_writer.finalize().unwrap();

        eprintln!(
            "{}",
            std::str::from_utf8(&annot_writer.output_stream).unwrap()
        );

        assert_eq!(
            std::str::from_utf8(&annot_writer.output_stream).unwrap(),
            r#"{
"metadata":
{
  "alt": "alternative sequence",
  "chr": "chromosome name same ase original vcf",
  "effect": "List of sequence ontology terms",
  "feature": "type of feature affected by variant gene/transcript",
  "id": "id of feature, same value of Id gff3 attributes",
  "impact": "0: UNKOWN, 1:LOW, 2:MODIFIER, 3: MODERATE, 4:HIGH",
  "name": "name of feature, same value of Name gff3 attributes",
  "pos": "position of variant",
  "ref": "reference sequence",
  "source": "source of variant in gff3 file"
},
"variants": [
{
  "variant": {
    "seqname": "93",
    "position": 2036067340,
    "ref_seq": "T",
    "alt_seq": "."
  },
  "myth": [
    {
      "source": "test",
      "feature": "gene",
      "id": "1111",
      "name": "gene1",
      "effects": [
        "GeneVariant",
        "ExonRegion"
      ],
      "impact": "Modifier"
    }
  ]
},
{
  "variant": {
    "seqname": "21",
    "position": 1970,
    "ref_seq": "C",
    "alt_seq": "T"
  },
  "myth": [
    {
      "source": "test2",
      "feature": "gene2",
      "id": "7777",
      "name": "gene51",
      "effects": [
        "DisruptiveInframeDeletion",
        "ExonRegion"
      ],
      "impact": "Moderate"
    }
  ]
}]
}"#
        );
    }

    #[test]
    fn test_write_myth_nd_json() {
        use crate::output::{JsonFormat, JsonWriter, MythWriter};
        let output_stream: Vec<u8> = Vec::new();
        let mut annot_writer = JsonWriter::new(output_stream, JsonFormat::NdJson).unwrap();

        let (myth, myth2) = get_two_myths();

        annot_writer.write_myth(myth).unwrap();
        annot_writer.write_myth(myth2).unwrap();
        annot_writer.finalize().unwrap();

        eprintln!(
            "{}",
            std::str::from_utf8(&annot_writer.output_stream).unwrap()
        );
        assert_eq!(
            std::str::from_utf8(&annot_writer.output_stream).unwrap(),
            r#"{"metadata":{"alt":"alternative sequence","chr":"chromosome name same ase original vcf","effect":"List of sequence ontology terms","feature":"type of feature affected by variant gene/transcript","id":"id of feature, same value of Id gff3 attributes","impact":"0: UNKOWN, 1:LOW, 2:MODIFIER, 3: MODERATE, 4:HIGH","name":"name of feature, same value of Name gff3 attributes","pos":"position of variant","ref":"reference sequence","source":"source of variant in gff3 file"}}
{"variant":{"seqname":"93","position":2036067340,"ref_seq":"T","alt_seq":"."},"myth":[{"source":"test","feature":"gene","id":"1111","name":"gene1","effects":["GeneVariant","ExonRegion"],"impact":"Modifier"}]}
{"variant":{"seqname":"21","position":1970,"ref_seq":"C","alt_seq":"T"},"myth":[{"source":"test2","feature":"gene2","id":"7777","name":"gene51","effects":["DisruptiveInframeDeletion","ExonRegion"],"impact":"Moderate"}]}"#
        );
    }
}
