//! A variant annotator.

#![warn(missing_docs)]

/* std use */

use std::io::BufWriter;

/* crate use */
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */

/* mod declaration */
pub mod annotation;
pub mod annotations_db;
pub mod cli;
pub mod effect;
pub mod error;
pub mod myth;
pub mod output;
pub mod sequences_db;
pub mod translate;
pub mod variant;
pub mod variant2myth;

use crate::output::writer::MythWriter;

/// Just a trait to combine Write and Seek trait
pub trait WriteSeek: std::io::Write + std::io::Seek {}
impl<T> WriteSeek for T where T: std::io::Write + std::io::Seek {}

#[cfg(not(feature = "parallel"))]
/// For each variants found matching annotations
pub fn vcf2json<R>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    mut vcf_reader: variant::VcfReader<R>,
    no_annotation: bool,
    block_size: usize,
    output: BufWriter<Box<dyn WriteSeek + Send>>,
    output_type: output::OutputFileType,
) -> error::Result<()>
where
    R: std::io::BufRead,
{
    let variant2myth =
        variant2myth::Variant2Myth::new(annotations, translate, sequences, no_annotation);

    let mut writer: Box<dyn MythWriter> = match output_type {
        output::OutputFileType::Parquet => {
            Box::new(output::parquet::ParquetWriter::new(output, block_size)?)
        }
        output::OutputFileType::JSON => Box::new(output::json::JsonWriter::new(output)?),
    };

    let mut end = true;
    while end {
        for _ in 0..block_size {
            if let Some(result) = vcf_reader.next() {
                let variant = result?;

                let myth = variant2myth.myth(variant);

                writer.write_myth(myth)?;
            } else {
                end = false;
                break;
            }
        }
    }
    writer.close()?;

    Ok(())
}

#[cfg(feature = "parallel")]
pub fn vcf2json<R>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    vcf_reader: variant::VcfReader<R>,
    no_annotation: bool,
    block_size: usize,
    output: BufWriter<Box<dyn WriteSeek + Send>>,
    output_type: output::OutputFileType,
) -> error::Result<()>
where
    R: std::io::BufRead + std::marker::Send,
{
    let (tx, rx) = std::sync::mpsc::channel::<myth::Myth>();

    let write_thread = std::thread::spawn(move || -> error::Result<()> {
        let mut writer: Box<dyn MythWriter> = match output_type {
            output::OutputFileType::Parquet => {
                Box::new(output::parquet::ParquetWriter::new(output, block_size)?)
            }
            output::OutputFileType::JSON => Box::new(output::json::JsonWriter::new(output)?),
        };

        'outer: loop {
            for _ in 0..block_size {
                if let Some(myth) = rx.iter().next() {
                    MythWriter::write_myth(writer.as_mut(), myth)?;
                } else {
                    break 'outer;
                }
            }
        }
        MythWriter::close(writer.as_mut())?;
        Ok(())
    });

    let variant2myth =
        variant2myth::Variant2Myth::new(annotations, translate, sequences, no_annotation);

    let results = vcf_reader
        .par_bridge()
        .filter(Result::is_ok)
        .map(error::Result::unwrap)
        .map(|variant| tx.send(variant2myth.myth(variant)))
        .filter(|r| r.is_err())
        .collect::<Vec<core::result::Result<(), std::sync::mpsc::SendError<myth::Myth>>>>();

    for result in results {
        result?
    }

    drop(tx);

    write_thread.join().unwrap()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */
    use biotest::Format as _;
    use bstr::ByteSlice as _;

    /* project use */
    use super::*;

    pub const GFF: &[u8] =
        b"{0}	ensembl_havana	gene	445	4439	.	+	.	ID=gene:ENSG00000112473;Name=SLC39A7
{0}	ensembl_havana	mRNA	826	4439	.	+	.	ID=transcript:ENST00000374675;Parent=gene:ENSG00000112473;Name=SLC39A7-201
{0}	ensembl_havana	exon	826	938	.	+	.	Parent=transcript:ENST00000374675;Name=ENSE00001746858
{0}	ensembl_havana	five_prime_UTR	826	938	.	+	.	Parent=transcript:ENST00000374675
{0}	ensembl_havana	five_prime_UTR	1242	1245	.	+	.	Parent=transcript:ENST00000374675
{0}	ensembl_havana	exon	1242	1656	.	+	.	Parent=transcript:ENST00000374675;Name=ENSE00001732423
{0}	ensembl_havana	CDS	1246	1656	.	+	0	ID=CDS:ENSP00000363807;Parent=transcript:ENST00000374675
{0}	ensembl_havana	exon	1745	1913	.	+	.	Parent=transcript:ENST00000374675;Name=ENSE00003321302
{0}	ensembl_havana	CDS	1745	1913	.	+	0	ID=CDS:ENSP00000363807;Parent=transcript:ENST00000374675
{0}	ensembl_havana	exon	2072	2125	.	+	.	Parent=transcript:ENST00000374675;Name=ENSE00003349240
{0}	ensembl_havana	CDS	2072	2125	.	+	2	ID=CDS:ENSP00000363807;Parent=transcript:ENST00000374675
{0}	ensembl_havana	exon	2263	2427	.	+	.	Parent=transcript:ENST00000374675;Name=ENSE00001700955
{0}	ensembl_havana	CDS	2263	2427	.	+	2	ID=CDS:ENSP00000363807;Parent=transcript:ENST00000374675
{0}	ensembl_havana	exon	2560	2700	.	+	.	Parent=transcript:ENST00000374675;Name=ENSE00003546616
{0}	ensembl_havana	CDS	2560	2700	.	+	2	ID=CDS:ENSP00000363807;Parent=transcript:ENST00000374675
{0}	ensembl_havana	exon	2910	3106	.	+	.	Parent=transcript:ENST00000374675;Name=ENSE00003693188
{0}	ensembl_havana	CDS	2910	3106	.	+	2	ID=CDS:ENSP00000363807;Parent=transcript:ENST00000374675
{0}	ensembl_havana	CDS	3541	3813	.	+	0	ID=CDS:ENSP00000363807;Parent=transcript:ENST00000374675
{0}	ensembl_havana	exon	3541	4439	.	+	.	Parent=transcript:ENST00000374675;Name=ENSE00003460070
{0}	ensembl_havana	three_prime_UTR	3814	4439	.	+	.	Parent=transcript:ENST00000374675";

    fn _setup() -> error::Result<(
        translate::Translate,
        sequences_db::SequencesDataBase,
        annotations_db::AnnotationsDataBase,
        Vec<u8>,
    )> {
        // generate fasta
        let mut rng = biotest::rand();
        let seq_gen = biotest::Fasta::builder().sequence_len(40_000).build()?;
        let mut fasta_reader = vec![];
        seq_gen.record(&mut fasta_reader, &mut rng)?;
        let reader: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>> =
            std::io::BufReader::new(Box::new(std::io::Cursor::new(fasta_reader.clone())));

        // produce gff
        let gff_reader: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>> =
            std::io::BufReader::new(Box::new(std::io::Cursor::new(
                GFF.replace(b"{0}", &fasta_reader[1..11]),
            )));
        let annotations_db = annotations_db::AnnotationsDataBase::from_reader(gff_reader, 100)?;

        // produce sequence db
        let seq_db = sequences_db::SequencesDataBase::from_reader(reader)?;
        // produce translate worker
        let translate = translate::Translate::default();

        Ok((
            translate,
            seq_db,
            annotations_db,
            fasta_reader[1..11].to_vec(),
        ))
    }
}
