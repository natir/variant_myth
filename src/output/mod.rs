pub mod json;
pub mod parquet;

pub mod writer;

pub enum OutputFileType {
    JSON,
    Parquet,
}

pub fn get_metadata() -> Vec<(&'static str, &'static str)> {
    vec![
        (
            "impact",
            "0: UNKOWN, 1:LOW, 2:MODIFIER, 3: MODERATE, 4:HIGH",
        ),
        ("effect", "List of sequence ontology terms"),
        ("gene_name", "HUGO symbol of affected gene"),
    ]
}
