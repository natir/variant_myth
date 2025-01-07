pub mod json;
pub mod parquet;

pub mod writer;

pub enum OutputFileType {
    JSON,
    Parquet,
}
