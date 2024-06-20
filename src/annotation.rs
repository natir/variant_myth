//! Struct to store annotation

/* std use */

/* crate use */
use bstr::ByteSlice as _;

/* project use */
use crate::error;

/// Define annotation strand
#[derive(Debug, Clone, Copy)]
pub enum Strand {
    /// Annotation is in forward strand
    Forward,
    /// Annotation is in reverse strand
    Reverse,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

/// Define annotation frame
#[derive(Debug, Clone, Copy)]
pub enum Frame {
    /// Annotation frame is Unknow
    Unknow,
    /// Annotation frame is 0
    Zero,
    /// Annotation frame is 1
    One,
    /// Annotation frame is 2
    Two,
}

impl std::fmt::Display for Frame {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        match self {
            Frame::Unknow => write!(f, "."),
            Frame::Zero => write!(f, "0"),
            Frame::One => write!(f, "1"),
            Frame::Two => write!(f, "2"),
        }
    }
}

/// Store attribute of gff record
#[derive(Debug, Clone, std::default::Default)]
pub struct Attribute {
    gene_id: Vec<u8>,
    transcript_id: Vec<u8>,
    gene_name: Option<Vec<u8>>,
    exon_number: Option<u64>,
    exon_id: Option<Vec<u8>>,
}

impl Attribute {
    /// Create an attribute from u8 slice
    pub fn from_u8_slice(slice: &[u8]) -> error::Result<Self> {
        let mut obj = Attribute::default();

        for attribute in slice.split_str(";") {
            match attribute {
                [b'g', b'e', b'n', b'e', b'_', b'i', b'd', b'=', value @ ..] => {
                    obj.gene_id = value.to_vec()
                }
                [b't', b'r', b'a', b'n', b's', b'c', b'r', b'i', b'p', b't', b'_', b'i', b'd', b'=', value @ ..] => {
                    obj.transcript_id = value.to_vec()
                }
                [b'g', b'e', b'n', b'e', b'_', b'n', b'a', b'm', b'e', b'=', value @ ..] => {
                    obj.gene_name = Some(value.to_vec())
                }
                [b'e', b'x', b'o', b'n', b'_', b'n', b'u', b'm', b'b', b'e', b'r', b'=', value @ ..] => {
                    obj.exon_number =
                        Some(unsafe { String::from_utf8_unchecked(value.to_vec()).parse::<u64>()? })
                }
                [b'e', b'x', b'o', b'n', b'_', b'i', b'd', b'=', value @ ..] => {
                    obj.exon_id = Some(value.to_vec())
                }
                _ => panic!(
                    "Not support attribute type {}",
                    String::from_utf8(attribute.to_vec()).unwrap()
                ),
            }
        }

        Ok(obj)
    }

    /// Get transcript_id
    pub fn get_transcript_id(&self) -> &[u8] {
        &self.transcript_id
    }

    /// Get exon_number
    pub fn get_exon_number(&self) -> u64 {
        self.exon_number.unwrap_or(0)
    }
}

impl std::fmt::Display for Attribute {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "gene_id={};transcript_id={}",
            String::from_utf8(self.gene_id.to_vec()).unwrap(),
            String::from_utf8(self.transcript_id.to_vec()).unwrap(),
        )?;
        if let Some(gn) = &self.gene_name {
            write!(f, ";gene_name={}", String::from_utf8(gn.to_vec()).unwrap())?;
        }
        if let Some(en) = &self.exon_number {
            write!(f, ";exon_number={}", en)?;
        }
        if let Some(ei) = &self.exon_id {
            write!(f, ";exon_id={}", String::from_utf8(ei.to_vec()).unwrap())?;
        }

        Ok(())
    }
}

#[derive(Debug, Clone)]
/// Store information of a Gff3 field
pub struct Annotation {
    seqname: Vec<u8>,
    source: Vec<u8>,
    feature: Vec<u8>,
    start: u64,
    end: u64,
    score: f64,
    strand: Strand,
    frame: Frame,
    attribute: Attribute,
}

impl Annotation {
    /// Build a new annotations from a csv::ByteRecord
    pub fn from_byte_record(record: &csv::ByteRecord) -> error::Result<Self> {
        Ok(Self {
            seqname: record.get(0).unwrap().to_vec(),
            source: record.get(1).unwrap().to_vec(),
            feature: record.get(2).unwrap().to_vec(),
            start: String::from_utf8(record.get(3).unwrap().to_vec())
                .unwrap()
                .parse::<u64>()
                .unwrap(),
            end: String::from_utf8(record.get(4).unwrap().to_vec())
                .unwrap()
                .parse::<u64>()
                .unwrap(),
            score: match record.get(5).unwrap() {
                b"." => f64::NAN,
                other => String::from_utf8(other.to_vec())
                    .unwrap()
                    .parse::<f64>()
                    .unwrap(),
            },
            strand: match record.get(6).unwrap() {
                b"+" => Strand::Forward,
                b"-" => Strand::Reverse,
                _ => return Err(error::Error::GffBadStrand.into()),
            },
            frame: match record.get(7).unwrap() {
                b"." => Frame::Unknow,
                b"0" => Frame::Zero,
                b"1" => Frame::One,
                b"2" => Frame::Two,
                _ => return Err(error::Error::GffBadFrame.into()),
            },
            attribute: Attribute::from_u8_slice(record.get(8).unwrap())?,
        })
    }

    /// Create interval associate with annotation
    pub fn get_interval(&self) -> core::ops::Range<u64> {
        self.start..self.end
    }

    /// Get seqname
    pub fn get_seqname(&self) -> &[u8] {
        &self.seqname
    }

    /// Get seqname
    pub fn get_source(&self) -> &[u8] {
        &self.source
    }

    /// Get strand
    pub fn get_strand(&self) -> &Strand {
        &self.strand
    }

    /// Get attribute annotation
    pub fn get_attribute(&self) -> &Attribute {
        &self.attribute
    }

    /// Get transcript id
    pub fn get_transcript_id(&self) -> &[u8] {
        self.attribute.get_transcript_id()
    }

    /// Get feature
    pub fn get_feature(&self) -> &[u8] {
        &self.feature
    }
}

impl std::fmt::Display for Annotation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{} {} {} {} {} {} {} {} {}",
            String::from_utf8(self.seqname.clone()).unwrap(),
            String::from_utf8(self.source.clone()).unwrap(),
            String::from_utf8(self.feature.clone()).unwrap(),
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            self.attribute,
        )
    }
}
