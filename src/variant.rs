//! Variant database

/* std use */

/* crate use */

/* project use */
use crate::error;

#[derive(Clone, serde::Serialize, PartialEq)]
/// Store Variant content
pub struct Variant {
    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Sequence name associate with variant
    pub seqname: Vec<u8>,

    /// Position of the variant
    pub position: u64,

    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Reference sequence associate with variant
    pub ref_seq: Vec<u8>,

    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Alternative sequence associate with variant
    pub alt_seq: Vec<u8>,
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
    pub fn get_interval(&self) -> core::ops::Range<u64> {
        self.position..self.position + self.ref_seq.len() as u64
    }
}

impl std::fmt::Debug for Variant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "Variant {{ seqname: b\"{}\".to_vec(), position: {}, ref_seq: b\"{}\".to_vec(), , alt_seq: b\"{}\".to_vec }}" ,
            String::from_utf8(self.seqname.to_vec()).unwrap(),
            self.position,
            String::from_utf8(self.ref_seq.to_vec()).unwrap(),
            String::from_utf8(self.alt_seq.to_vec()).unwrap()
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

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */
    use biotest::Format as _;

    /* project use */
    use super::*;

    #[test]
    fn variant() -> error::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Vcf::default();

        let mut input = vec![];
        generator.records(&mut input, &mut rng, 5)?;

        let reader = VcfReader::from_reader(std::io::Cursor::new(input));

        let records = reader
            .into_iter()
            .collect::<error::Result<Vec<Variant>>>()?;

        assert_eq!(
            records,
            vec![
                Variant {
                    seqname: b"YAR028W".to_vec(),
                    position: 509242864,
                    ref_seq: b"A".to_vec(),
                    alt_seq: b".".to_vec(),
                },
                Variant {
                    seqname: b"93".to_vec(),
                    position: 2036067340,
                    ref_seq: b"T".to_vec(),
                    alt_seq: b".".to_vec()
                },
                Variant {
                    seqname: b"X".to_vec(),
                    position: 2138516245,
                    ref_seq: b"A".to_vec(),
                    alt_seq: b".".to_vec()
                },
                Variant {
                    seqname: b"NC_000015.10".to_vec(),
                    position: 1204106469,
                    ref_seq: b"c".to_vec(),
                    alt_seq: b".".to_vec()
                },
                Variant {
                    seqname: b"NC_016845.1".to_vec(),
                    position: 1745241132,
                    ref_seq: b"c".to_vec(),
                    alt_seq: b".".to_vec()
                }
            ]
        );

        assert_eq!(records[1].get_interval(), (2036067340u64..2036067341u64));

        assert_eq!(format!("{:?}", records[2]), "Variant { seqname: b\"X\".to_vec(), position: 2138516245, ref_seq: b\"A\".to_vec(), , alt_seq: b\".\".to_vec }".to_string());
        Ok(())
    }
}
