//! A module to store and prepare test data

/* std use */

use core::f64;

/* crate use */
use bstr::ByteSlice;

/* project use */
use crate::{annotation, variant};

/// GFF file
pub const GFF: &[u8] = std::include_bytes!("test_data/annotations.gff3");

/// Sequence file
pub const SEQUENCE: &[u8] = std::include_bytes!("test_data/references.fasta");

/// Variant file
pub const VARIANT: &[u8] = std::include_bytes!("test_data/variants.vcf");

/// GFF file split by line
pub static GFF_BY_LINE: std::sync::LazyLock<Vec<Vec<u8>>> =
    std::sync::LazyLock::new(|| GFF.split_str("\n").map(|line| line.to_vec()).collect());
/// GFF csv record
pub static GFF_CSV_RECORD: std::sync::LazyLock<Vec<csv::ByteRecord>> =
    std::sync::LazyLock::new(|| {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(GFF);

        reader.byte_records().map(|r| r.unwrap()).collect()
    });
/// GFF annotation
pub static GFF_ANNOTATION: std::sync::LazyLock<Vec<annotation::Annotation>> =
    std::sync::LazyLock::new(|| {
        vec![
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"gene".to_vec(),
                start: 51,
                stop: 30235,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Forward,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"ENSG00000286586.2".to_vec(),
                    name: b"".to_vec(),
                    parent: b"".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"transcript".to_vec(),
                start: 51,
                stop: 30235,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Forward,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"ENST00000797271.1".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENSG00000286586.2".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"five_prime_UTR".to_vec(),
                start: 51,
                stop: 61,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Forward,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"UTR5:ENST00000797271.1".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENST00000797271.1".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"exon".to_vec(),
                start: 61,
                stop: 261,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Forward,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"exon:ENST00000797271.1:1".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENST00000797271.1".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"exon".to_vec(),
                start: 13202,
                stop: 13359,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Forward,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"exon:ENST00000797271.1:2".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENST00000797271.1".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"exon".to_vec(),
                start: 25008,
                stop: 25211,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Forward,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"exon:ENST00000797271.1:3".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENST00000797271.1".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"exon".to_vec(),
                start: 26180,
                stop: 26358,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Forward,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"exon:ENST00000797271.1:4".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENST00000797271.1".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"exon".to_vec(),
                start: 29612,
                stop: 30235,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Forward,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"exon:ENST00000797271.1:5".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENST00000797271.1".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"gene".to_vec(),
                start: 121694345,
                stop: 121695599,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Reverse,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"ENSG00000309035.1".to_vec(),
                    name: b"".to_vec(),
                    parent: b"".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"transcript".to_vec(),
                start: 121694345,
                stop: 121695599,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Reverse,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"ENST00000837983.1".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENSG00000309035.1".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"exon".to_vec(),
                start: 121695271,
                stop: 121695599,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Reverse,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"exon:ENST00000837983.1:1".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENST00000837983.1".to_vec(),
                },
            },
            annotation::Annotation {
                seqname: b"chrA".to_vec(),
                source: b"HAVANA".to_vec(),
                feature: b"exon".to_vec(),
                start: 121694345,
                stop: 121694556,
                score: annotation::Score(f64::INFINITY),
                strand: annotation::Strand::Reverse,
                frame: annotation::Frame::Unknow,
                attribute: annotation::Attribute {
                    id: b"exon:ENST00000837983.1:2".to_vec(),
                    name: b"".to_vec(),
                    parent: b"ENST00000837983.1".to_vec(),
                },
            },
        ]
    });

/// Variant record
pub static VARIANT_RECORD: std::sync::LazyLock<Vec<variant::Variant>> =
    std::sync::LazyLock::new(|| {
        vec![
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 5772,
                ref_seq: b"A".to_vec(),
                alt_seq: b"T".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 18217,
                ref_seq: b"C".to_vec(),
                alt_seq: b"CT".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 12883,
                ref_seq: b"A".to_vec(),
                alt_seq: b"G".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 9035,
                ref_seq: b"C".to_vec(),
                alt_seq: b"G".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 6736,
                ref_seq: b"A".to_vec(),
                alt_seq: b"G".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 121694676,
                ref_seq: b"A".to_vec(),
                alt_seq: b"C".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 121694675,
                ref_seq: b"C".to_vec(),
                alt_seq: b"G".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 121694690,
                ref_seq: b"T".to_vec(),
                alt_seq: b"TCTC".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 121694674,
                ref_seq: b"C".to_vec(),
                alt_seq: b"T".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 121695123,
                ref_seq: b"G".to_vec(),
                alt_seq: b"A".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 28567,
                ref_seq: b"T".to_vec(),
                alt_seq: b"g".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 9685,
                ref_seq: b"A".to_vec(),
                alt_seq: b"c".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 15351,
                ref_seq: b"G".to_vec(),
                alt_seq: b"Gtgcg".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 19167,
                ref_seq: b"T".to_vec(),
                alt_seq: b"Tgg".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 14832,
                ref_seq: b"ACTT".to_vec(),
                alt_seq: b"A".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 29505,
                ref_seq: b"AA".to_vec(),
                alt_seq: b"A".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 16413,
                ref_seq: b"TC".to_vec(),
                alt_seq: b"T".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 22846,
                ref_seq: b"GCTCT".to_vec(),
                alt_seq: b"G".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 13205,
                ref_seq: b"T".to_vec(),
                alt_seq: b"Tggat".to_vec(),
                variant_type: variant::Type::Small,
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 6650,
                ref_seq: b"A".to_vec(),
                alt_seq: b"<CNV>".to_vec(),
                variant_type: variant::Type::Cnv(589),
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 11777,
                ref_seq: b"C".to_vec(),
                alt_seq: b"<INS>".to_vec(),
                variant_type: variant::Type::Ins(613),
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 28302,
                ref_seq: b"C".to_vec(),
                alt_seq: b"<DEL>".to_vec(),
                variant_type: variant::Type::Del(974),
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 3334,
                ref_seq: b"A".to_vec(),
                alt_seq: b"<DUP>".to_vec(),
                variant_type: variant::Type::Dup(753),
            },
            variant::Variant {
                seqname: b"chrA".to_vec(),
                position: 3279,
                ref_seq: b"A".to_vec(),
                alt_seq: b"<CNV>".to_vec(),
                variant_type: variant::Type::Ins(907),
            },
        ]
    });
