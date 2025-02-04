//! Variant database

/* std use */

/* crate use */
use bstr::ByteSlice as _;

/* project use */
use crate::error;

#[derive(Clone, PartialEq, Debug)]
/// Store type of variant
pub enum Type {
    /// Variant are small (alt_seq store a sequence)
    Small,
    /// Variant are an insertion (alt_seq contains <INS>)
    Ins(u64),
    /// Variant are a deletion (alt_seq contains <DEL>)
    Del(u64),
    /// Variant are a duplication (alt_seq contains <DUP>)
    Dup(u64),
    /// Variant are an inversion (alt_seq contains <INV>)
    Inv(u64),
    /// Variant are a copy number variation (alt_seq contains <CNV>)
    Cnv(u64),
}

impl Type {
    /// Create a Type from alt_seq
    pub fn from_alt(alt: &[u8], optional_infos: Option<&[u8]>) -> error::Result<Self> {
        let mut opt_svlen = None;
        if let Some(infos) = optional_infos {
            for info in infos.split_str(";") {
                if let [b'S', b'V', b'L', b'E', b'N', b'=', value @ ..] = info {
                    opt_svlen =
                        unsafe { Some(String::from_utf8_unchecked(value.to_vec()).parse::<u64>()?) }
                }
            }
        }

        match alt {
            b"<INS>" => {
                if let Some(svlen) = opt_svlen {
                    Ok(Type::Ins(svlen))
                } else {
                    Err(error::Error::VcfStructVariantNoSvLen.into())
                }
            }
            b"<DEL>" => {
                if let Some(svlen) = opt_svlen {
                    Ok(Type::Del(svlen))
                } else {
                    Err(error::Error::VcfStructVariantNoSvLen.into())
                }
            }
            b"<DUP>" => {
                if let Some(svlen) = opt_svlen {
                    Ok(Type::Dup(svlen))
                } else {
                    Err(error::Error::VcfStructVariantNoSvLen.into())
                }
            }
            b"<INV>" => {
                if let Some(svlen) = opt_svlen {
                    Ok(Type::Inv(svlen))
                } else {
                    Err(error::Error::VcfStructVariantNoSvLen.into())
                }
            }
            b"<CNV>" => {
                if let Some(svlen) = opt_svlen {
                    Ok(Type::Cnv(svlen))
                } else {
                    Err(error::Error::VcfStructVariantNoSvLen.into())
                }
            }
            _ => Ok(Type::Small),
        }
    }
}

#[derive(Clone, PartialEq)]
#[cfg_attr(feature = "out_json", derive(serde::Serialize, serde::Deserialize))]
/// Store Variant content
pub struct Variant {
    /// Sequence name associate with variant
    #[cfg_attr(feature = "out_json", serde(serialize_with = "crate::serialize_bstr"))]
    pub seqname: Vec<u8>,

    /// Position of the variant
    pub position: u64,

    /// Reference sequence associate with variant
    #[cfg_attr(feature = "out_json", serde(serialize_with = "crate::serialize_bstr"))]
    pub ref_seq: Vec<u8>,

    /// Alternative sequence associate with variant
    #[cfg_attr(feature = "out_json", serde(serialize_with = "crate::serialize_bstr"))]
    pub alt_seq: Vec<u8>,

    /// Store the type of variant
    #[cfg_attr(feature = "out_json", serde(skip_serializing))]
    pub variant_type: Type,
}

impl Variant {
    /// Build a Variant from csv::ByteRecord
    pub fn from_byte_record(record: csv::ByteRecord) -> error::Result<Self> {
        let seqname = record.get(0).ok_or(error::Error::VcfBadRecord)?.to_vec();
        let position = unsafe {
            String::from_utf8_unchecked(record.get(1).ok_or(error::Error::VcfBadRecord)?.to_vec())
                .parse::<u64>()?
        } - 1;
        let ref_seq = record.get(3).ok_or(error::Error::VcfBadRecord)?.to_vec();
        let alt_seq = record.get(4).ok_or(error::Error::VcfBadRecord)?.to_vec();
        let info = record.get(7);

        let variant_type = Type::from_alt(&alt_seq, info)?;

        Ok(Self {
            seqname,
            position,
            ref_seq,
            alt_seq,
            variant_type,
        })
    }

    /// Create interval associate with variant
    pub fn get_interval(&self) -> core::ops::Range<u64> {
        match self.variant_type {
            Type::Small => self.position..self.position + self.ref_seq.len() as u64,
            Type::Ins(_size) => self.position..self.position + 1,
            Type::Del(size) => self.position..self.position + size,
            Type::Dup(size) => self.position..self.position + size * 2,
            Type::Inv(size) => self.position..self.position + size,
            Type::Cnv(size) => self.position..self.position + size,
        }
    }

    /// Variant can be annotate by variant_myth
    ///
    /// ref_seq and alt_seq must contains at least one A,C,T,G,a,c,t or g
    pub fn valid(&self) -> bool {
        let mut nuc_iter = if self.structural() {
            Box::new(self.ref_seq.iter()) as Box<dyn std::iter::Iterator<Item = &u8>>
        } else {
            Box::new(self.ref_seq.iter().chain(self.alt_seq.iter()))
                as Box<dyn std::iter::Iterator<Item = &u8>>
        };

        (!self.ref_seq.is_empty())
            && (!self.alt_seq.is_empty())
            && nuc_iter.all(|c| {
                c == &b'A'
                    || c == &b'C'
                    || c == &b'T'
                    || c == &b'G'
                    || c == &b'a'
                    || c == &b'c'
                    || c == &b't'
                    || c == &b'g'
            })
    }

    /// Variant is structural variant
    pub fn structural(&self) -> bool {
        self.variant_type != Type::Small
    }

    #[cfg(test)]
    /// Generate a fake variant with a seqname, position, ref_seq, alt_seqdb
    pub fn test_variant(
        seqname: &[u8],
        position: u64,
        ref_seq: &[u8],
        alt_seq: &[u8],
        opt_info: Option<&[u8]>,
    ) -> error::Result<Variant> {
        Ok(Self {
            seqname: seqname.to_vec(),
            position,
            ref_seq: ref_seq.to_vec(),
            alt_seq: alt_seq.to_vec(),
            variant_type: Type::from_alt(alt_seq, opt_info)?,
        })
    }
}

impl std::fmt::Debug for Variant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "Variant {{ seqname: b\"{}\".to_vec(), position: {}, ref_seq: b\"{}\".to_vec(), alt_seq: b\"{}\".to_vec }}" ,
            String::from_utf8(self.seqname.to_vec()).unwrap(),
            self.position,
            String::from_utf8(self.ref_seq.to_vec()).unwrap(),
            String::from_utf8(self.alt_seq.to_vec()).unwrap()
        )
    }
}

impl std::fmt::Display for Variant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}-{}-{}-{}",
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
                    position: 509242863,
                    ref_seq: b"A".to_vec(),
                    alt_seq: b".".to_vec(),
                    variant_type: Type::Small,
                },
                Variant {
                    seqname: b"93".to_vec(),
                    position: 2036067339,
                    ref_seq: b"T".to_vec(),
                    alt_seq: b".".to_vec(),
                    variant_type: Type::Small,
                },
                Variant {
                    seqname: b"X".to_vec(),
                    position: 2138516244,
                    ref_seq: b"A".to_vec(),
                    alt_seq: b".".to_vec(),
                    variant_type: Type::Small,
                },
                Variant {
                    seqname: b"NC_000015.10".to_vec(),
                    position: 1204106468,
                    ref_seq: b"c".to_vec(),
                    alt_seq: b".".to_vec(),
                    variant_type: Type::Small,
                },
                Variant {
                    seqname: b"NC_016845.1".to_vec(),
                    position: 1745241131,
                    ref_seq: b"c".to_vec(),
                    alt_seq: b".".to_vec(),
                    variant_type: Type::Small,
                }
            ]
        );

        assert_eq!(records[1].get_interval(), (2036067339u64..2036067340u64));

        assert_eq!(format!("{:?}", records[2]), "Variant { seqname: b\"X\".to_vec(), position: 2138516244, ref_seq: b\"A\".to_vec(), alt_seq: b\".\".to_vec }".to_string());

        assert_eq!(format!("{}", records[2]), "X-2138516244-A-.".to_string());

        Ok(())
    }

    #[test]
    fn test_variant() -> error::Result<()> {
        let variant = Variant::test_variant(b"chr1", 62103, b"ACT", b"A", None)?; // 0-based

        assert_eq!(variant.seqname, b"chr1");
        assert_eq!(variant.position, 62103); // 0-based
        assert_eq!(variant.ref_seq, b"ACT");
        assert_eq!(variant.alt_seq, b"A");
        assert_eq!(variant.variant_type, Type::Small);

        Ok(())
    }

    #[test]
    fn validate_variant() -> error::Result<()> {
        let variant = Variant::test_variant(b"chr1", 62103, b"ACT", b"A", None)?;
        assert!(variant.valid());
        assert!(!variant.structural());

        let variant = Variant::test_variant(b"chr1", 62103, b"act", b"a", None)?;
        assert!(variant.valid());
        assert!(!variant.structural());

        // empty alt
        let variant = Variant::test_variant(b"chr1", 62103, b"a", b"", None)?;
        assert!(!variant.valid());
        assert!(!variant.structural());

        // empty ref
        let variant = Variant::test_variant(b"chr1", 62103, b"", b"a", None)?;
        assert!(!variant.valid());
        assert!(!variant.structural());

        // not valid ref
        let variant = Variant::test_variant(b"chr1", 62103, b"AAA,A", b"a", None)?;
        assert!(!variant.valid());
        assert!(!variant.structural());

        // not valid alt
        let variant = Variant::test_variant(b"chr1", 62103, b"A", b"a,A", None)?;
        assert!(!variant.valid());
        assert!(!variant.structural());

        Ok(())
    }

    #[test]
    fn struct_variant() -> error::Result<()> {
        let variant = Variant::from_byte_record(csv::ByteRecord::from(vec![
            "chr1", "62104", ".", "ACT", "<INS>", ".", ".", "SVLEN=40",
        ]))?;
        assert_eq!(variant.seqname, b"chr1");
        assert_eq!(variant.position, 62103); // 0-based
        assert_eq!(variant.ref_seq, b"ACT");
        assert_eq!(variant.alt_seq, b"<INS>");
        assert_eq!(variant.variant_type, Type::Ins(40));
        assert_eq!(variant.get_interval(), 62103..62104);
        assert!(variant.valid());
        assert!(variant.structural());

        let variant = Variant::test_variant(b"chr1", 62103, b"ACT", b"<INS>", Some(b"SVLEN=40"))?;
        assert_eq!(variant.seqname, b"chr1");
        assert_eq!(variant.position, 62103); // 0-based
        assert_eq!(variant.ref_seq, b"ACT");
        assert_eq!(variant.alt_seq, b"<INS>");
        assert_eq!(variant.variant_type, Type::Ins(40));
        assert_eq!(variant.get_interval(), 62103..62104);
        assert!(variant.valid());
        assert!(variant.structural());

        let variant = Variant::test_variant(b"chr1", 62103, b"ACT", b"<DEL>", Some(b"SVLEN=40"))?;
        assert_eq!(variant.seqname, b"chr1");
        assert_eq!(variant.position, 62103); // 0-based
        assert_eq!(variant.ref_seq, b"ACT");
        assert_eq!(variant.alt_seq, b"<DEL>");
        assert_eq!(variant.variant_type, Type::Del(40));
        assert_eq!(variant.get_interval(), 62103..62143);
        assert!(variant.valid());
        assert!(variant.structural());

        let variant = Variant::test_variant(b"chr1", 62103, b"ACT", b"<DUP>", Some(b"SVLEN=40"))?;
        assert_eq!(variant.seqname, b"chr1");
        assert_eq!(variant.position, 62103); // 0-based
        assert_eq!(variant.ref_seq, b"ACT");
        assert_eq!(variant.alt_seq, b"<DUP>");
        assert_eq!(variant.variant_type, Type::Dup(40));
        assert_eq!(variant.get_interval(), 62103..62183);
        assert!(variant.valid());
        assert!(variant.structural());

        let variant = Variant::test_variant(b"chr1", 62103, b"ACT", b"<INV>", Some(b"SVLEN=40"))?;
        assert_eq!(variant.seqname, b"chr1");
        assert_eq!(variant.position, 62103); // 0-based
        assert_eq!(variant.ref_seq, b"ACT");
        assert_eq!(variant.alt_seq, b"<INV>");
        assert_eq!(variant.variant_type, Type::Inv(40));
        assert_eq!(variant.get_interval(), 62103..62143);
        assert!(variant.valid());
        assert!(variant.structural());

        let variant = Variant::test_variant(b"chr1", 62103, b"ACT", b"<CNV>", Some(b"SVLEN=40"))?;
        assert_eq!(variant.seqname, b"chr1");
        assert_eq!(variant.position, 62103); // 0-based
        assert_eq!(variant.ref_seq, b"ACT");
        assert_eq!(variant.alt_seq, b"<CNV>");
        assert_eq!(variant.variant_type, Type::Cnv(40));
        assert_eq!(variant.get_interval(), 62103..62143);
        assert!(variant.valid());
        assert!(variant.structural());

        let miss_svlen = Variant::from_byte_record(csv::ByteRecord::from(vec![
            "chr1", "62104", ".", "ACT", "<INS>", ".", ".", "SVLE=40",
        ]));
        assert!(miss_svlen.is_err());

        let miss_svlen = Variant::from_byte_record(csv::ByteRecord::from(vec![
            "chr1", "62104", ".", "ACT", "<DEL>", ".", ".", "SVLE=40",
        ]));
        assert!(miss_svlen.is_err());

        let miss_svlen = Variant::from_byte_record(csv::ByteRecord::from(vec![
            "chr1", "62104", ".", "ACT", "<DUP>", ".", ".", "SVLE=40",
        ]));
        assert!(miss_svlen.is_err());

        let miss_svlen = Variant::from_byte_record(csv::ByteRecord::from(vec![
            "chr1", "62104", ".", "ACT", "<INV>", ".", ".", "SVLE=40",
        ]));
        assert!(miss_svlen.is_err());

        let miss_svlen = Variant::from_byte_record(csv::ByteRecord::from(vec![
            "chr1", "62104", ".", "ACT", "<CNV>", ".", ".", "SVLE=40",
        ]));
        assert!(miss_svlen.is_err());

        let bad_alt = Variant::from_byte_record(csv::ByteRecord::from(vec![
            "chr1",
            "62104",
            ".",
            "ACT",
            "<ANYTHING>",
            ".",
            ".",
            "SVLEN=40",
        ]))?;
        assert!(!bad_alt.valid());

        Ok(())
    }
}
