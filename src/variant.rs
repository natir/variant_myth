//! Variant database

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Store Variant content
pub struct Variant {
    seqname: Vec<u8>,
    position: u64,
    ref_seq: Vec<u8>,
    alt_seq: Vec<u8>,
}

impl Variant {
    /// Build a Variant from csv::ByteRecord
    pub fn from_byte_record(record: csv::ByteRecord) -> error::Result<Self> {
        Ok(Self {
            seqname: record.get(0).ok_or(error::Error::VcfBadRecord)?.to_vec(),
            position: unsafe {
                String::from_utf8_unchecked(
                    record.get(1).ok_or(error::Error::VcfBadRecord)?.to_vec(),
                )
                .parse::<u64>()?
            },
            ref_seq: record.get(3).ok_or(error::Error::VcfBadRecord)?.to_vec(),
            alt_seq: record.get(4).ok_or(error::Error::VcfBadRecord)?.to_vec(),
        })
    }

    /// Create interval associate with variant
    pub fn get_interval(&self) -> (&[u8], core::ops::Range<u64>) {
        (
            &self.seqname,
            self.position..self.position + self.ref_seq.len() as u64,
        )
    }
}

/// Struct to generate Variant iterator from vcf
pub struct VcfReader<R>
where
    R: std::io::BufRead,
{
    inner: csv::Reader<R>,
}

impl<R> VcfReader<R>
where
    R: std::io::BufRead,
{
    /// Create a VcfReader from a std::io::BufRead
    pub fn from_reader(inner: R) -> Self {
        Self {
            inner: csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'#'))
                .from_reader(inner),
        }
    }
}

impl<R> std::iter::Iterator for VcfReader<R>
where
    R: std::io::BufRead,
{
    type Item = error::Result<Variant>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = csv::ByteRecord::new();
        match self.inner.read_byte_record(&mut record) {
            Ok(true) => Some(Variant::from_byte_record(record)),
            Ok(false) => None,
            Err(e) => Some(Err(e.into())),
        }
    }
}
