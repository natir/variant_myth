//! Module to manage output

/* std use */

/* crate use */

/* project use */

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
