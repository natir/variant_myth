//! Module to manage output

/* std use */

/* crate use */

/* project use */

use crate::error::Result;
use crate::myth::Myth;
use std::io::Write;
use std::marker::Send;

use crate::output::writer::MythWriter;

/// Get schema of parquet output
pub fn schema() -> arrow::datatypes::Schema {
    let mut fields = vec![
        arrow::datatypes::Field::new("chr", arrow::datatypes::DataType::Utf8, false),
        arrow::datatypes::Field::new("pos", arrow::datatypes::DataType::UInt64, false),
        arrow::datatypes::Field::new("ref", arrow::datatypes::DataType::Utf8, false),
        arrow::datatypes::Field::new("alt", arrow::datatypes::DataType::Utf8, false),
    ];

    fields.extend(vec![
        arrow::datatypes::Field::new("source", arrow::datatypes::DataType::Utf8, true),
        arrow::datatypes::Field::new("transcript_id", arrow::datatypes::DataType::Utf8, true),
        arrow::datatypes::Field::new("gene_name", arrow::datatypes::DataType::Utf8, true),
        arrow::datatypes::Field::new("effects", arrow::datatypes::DataType::Utf8, true),
        arrow::datatypes::Field::new("impact", arrow::datatypes::DataType::UInt8, true),
    ]);

    arrow::datatypes::Schema::new(fields)
}

pub struct ParquetWriter<W: Write + std::marker::Send + std::io::Seek + 'static> {
    writer: parquet::arrow::arrow_writer::ArrowWriter<W>,
    schema: std::sync::Arc<arrow::datatypes::Schema>,
    chrs: Vec<String>,
    poss: Vec<u64>,
    refs: Vec<String>,
    alts: Vec<String>,
    source: Vec<String>,
    transcript_id: Vec<String>,
    gene_name: Vec<String>,
    effects: Vec<String>,
    impact: Vec<u8>,
}

impl<W: Write + Send + std::io::Seek + 'static> ParquetWriter<W> {
    pub fn new(writer: W, block_size: usize) -> Result<Self> {
        let schema = std::sync::Arc::new(schema());
        let writer = parquet::arrow::arrow_writer::ArrowWriter::try_new(
            writer,
            schema.clone(),
            Default::default(),
        )?;
        Ok(ParquetWriter {
            writer,
            schema,
            chrs: Vec::with_capacity(block_size),
            poss: Vec::with_capacity(block_size),
            refs: Vec::with_capacity(block_size),
            alts: Vec::with_capacity(block_size),
            source: Vec::with_capacity(block_size),
            transcript_id: Vec::with_capacity(block_size),
            gene_name: Vec::with_capacity(block_size),
            effects: Vec::with_capacity(block_size),
            impact: Vec::with_capacity(block_size),
        })
    }
}

impl<W: Write + Send + std::io::Seek + 'static> MythWriter for ParquetWriter<W> {
    fn write_myth(&mut self, myth: Myth) -> Result<()> {
        for annotation in myth.annotations {
            self.chrs
                .push(unsafe { String::from_utf8_unchecked(myth.variant.seqname.clone()) });
            self.poss.push(myth.variant.position);
            self.refs
                .push(unsafe { String::from_utf8_unchecked(myth.variant.ref_seq.clone()) });
            self.alts
                .push(unsafe { String::from_utf8_unchecked(myth.variant.alt_seq.clone()) });
            self.source
                .push(unsafe { String::from_utf8_unchecked(annotation.source) });
            self.transcript_id
                .push(unsafe { String::from_utf8_unchecked(annotation.transcript_id) });
            self.gene_name
                .push(unsafe { String::from_utf8_unchecked(annotation.gene_name) });
            self.effects.push(
                annotation
                    .effects
                    .iter()
                    .map(|e| unsafe { String::from_utf8_unchecked(e.clone().into()) })
                    .collect::<Vec<String>>()
                    .join(";"),
            );
            self.impact.push(annotation.impact as u8);
        }
        Ok(())
    }
    fn end_batch(&mut self) -> Result<()> {
        let batch = arrow::record_batch::RecordBatch::try_new(
            self.schema.clone(),
            vec![
                std::sync::Arc::new(arrow::array::StringArray::from(std::mem::take(
                    &mut self.chrs,
                ))),
                std::sync::Arc::new(arrow::array::UInt64Array::from(std::mem::take(
                    &mut self.poss,
                ))),
                std::sync::Arc::new(arrow::array::StringArray::from(std::mem::take(
                    &mut self.refs,
                ))),
                std::sync::Arc::new(arrow::array::StringArray::from(std::mem::take(
                    &mut self.alts,
                ))),
                std::sync::Arc::new(arrow::array::StringArray::from(std::mem::take(
                    &mut self.source,
                ))),
                std::sync::Arc::new(arrow::array::StringArray::from(std::mem::take(
                    &mut self.transcript_id,
                ))),
                std::sync::Arc::new(arrow::array::StringArray::from(std::mem::take(
                    &mut self.gene_name,
                ))),
                std::sync::Arc::new(arrow::array::StringArray::from(std::mem::take(
                    &mut self.effects,
                ))),
                std::sync::Arc::new(arrow::array::UInt8Array::from(std::mem::take(
                    &mut self.impact,
                ))),
            ],
        )?;

        self.writer.write(&batch)?;
        Ok(())
    }
    fn close(self) -> Result<()> {
        self.writer.close()?;
        Ok(())
    }
}
