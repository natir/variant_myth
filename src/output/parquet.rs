//! Module to manage output

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::myth;
use crate::output;

fn get_metadata() -> Vec<parquet::file::metadata::KeyValue> {
    crate::output::get_metadata()
        .iter()
        .map(|(k, v)| {
            parquet::file::metadata::KeyValue::new(String::from(*k), Some(String::from(*v)))
        })
        .collect()
}

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
        arrow::datatypes::Field::new("feature", arrow::datatypes::DataType::Utf8, true),
        arrow::datatypes::Field::new("name", arrow::datatypes::DataType::Utf8, true),
        arrow::datatypes::Field::new("id", arrow::datatypes::DataType::Utf8, true),
        arrow::datatypes::Field::new("effects", arrow::datatypes::DataType::Utf8, true),
        arrow::datatypes::Field::new("impact", arrow::datatypes::DataType::UInt8, true),
    ]);

    arrow::datatypes::Schema::new(fields)
}

/// Struct to write Myth in parquet format
pub struct ParquetWriter<W: std::io::Write + std::marker::Send + std::io::Seek + 'static> {
    writer: parquet::arrow::arrow_writer::ArrowWriter<W>,
    schema: std::sync::Arc<arrow::datatypes::Schema>,
    chrs: Vec<String>,
    poss: Vec<u64>,
    refs: Vec<String>,
    alts: Vec<String>,
    source: Vec<String>,
    feature: Vec<String>,
    name: Vec<String>,
    id: Vec<String>,
    effects: Vec<String>,
    impact: Vec<u8>,
    block_size: usize,
}

impl<W: std::io::Write + Send + std::io::Seek + 'static> ParquetWriter<W> {
    /// Create a new ParquetWriter
    pub fn new(writer: W, block_size: usize) -> error::Result<Self> {
        let schema = std::sync::Arc::new(schema());

        let columns_metadata = get_metadata();

        let writer = parquet::arrow::arrow_writer::ArrowWriter::try_new(
            writer,
            schema.clone(),
            Some(
                parquet::file::properties::WriterProperties::builder()
                    .set_key_value_metadata(Some(columns_metadata))
                    .build(),
            ),
        )?;
        Ok(ParquetWriter {
            writer,
            schema,
            chrs: Vec::with_capacity(block_size),
            poss: Vec::with_capacity(block_size),
            refs: Vec::with_capacity(block_size),
            alts: Vec::with_capacity(block_size),
            source: Vec::with_capacity(block_size),
            feature: Vec::with_capacity(block_size),
            name: Vec::with_capacity(block_size),
            id: Vec::with_capacity(block_size),
            effects: Vec::with_capacity(block_size),
            impact: Vec::with_capacity(block_size),
            block_size,
        })
    }
}

impl<W: std::io::Write + std::marker::Send + std::io::Seek + 'static> output::MythWriter
    for ParquetWriter<W>
{
    fn add_myth(&mut self, myth: myth::Myth) -> error::Result<()> {
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
            self.feature
                .push(unsafe { String::from_utf8_unchecked(annotation.feature) });
            self.name
                .push(unsafe { String::from_utf8_unchecked(annotation.name) });
            self.id
                .push(unsafe { String::from_utf8_unchecked(annotation.id) });
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

    fn batch_full(&self) -> bool {
        self.chrs.len() > self.block_size
    }
    fn finalize(&mut self) -> error::Result<()> {
        self.writer.finish()?;
        Ok(())
    }
    fn write_batch(&mut self) -> error::Result<()> {
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
                    &mut self.feature,
                ))),
                std::sync::Arc::new(arrow::array::StringArray::from(std::mem::take(
                    &mut self.name,
                ))),
                std::sync::Arc::new(arrow::array::StringArray::from(std::mem::take(
                    &mut self.id,
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
}
