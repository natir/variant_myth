//! A variant annotator.

#![warn(missing_docs)]

/* std use */

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
pub fn vcf2json<R, W>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    mut vcf_reader: variant::VcfReader<R>,
    no_annotation: bool,
    mut output: W,
) -> error::Result<()>
where
    R: std::io::BufRead,
    W: std::io::Write + std::io::Seek + std::marker::Send + 'static,
{
    let variant2myth =
        variant2myth::Variant2Myth::new(annotations, translate, sequences, no_annotation);

    let schema = std::sync::Arc::new(output::parquet::schema());
    let mut writer = parquet::arrow::arrow_writer::ArrowWriter::try_new(
        &mut output,
        schema.clone(),
        Default::default(),
    )?;

    let block_size = 1 << 13;

    let mut end = true;
    while end {
        let mut chrs = Vec::with_capacity(block_size);
        let mut poss = Vec::with_capacity(block_size);
        let mut refs = Vec::with_capacity(block_size);
        let mut alts = Vec::with_capacity(block_size);
        let mut source = Vec::with_capacity(block_size);
        let mut transcript_id = Vec::with_capacity(block_size);
        let mut gene_name = Vec::with_capacity(block_size);
        let mut effects = Vec::with_capacity(block_size);
        let mut impact: Vec<u8> = Vec::with_capacity(block_size);

        for _ in 0..block_size {
            if let Some(result) = vcf_reader.next() {
                let variant = result?;

                let myth = variant2myth.myth(variant);

                for annotation in myth.annotations {
                    chrs.push(unsafe { String::from_utf8_unchecked(myth.variant.seqname.clone()) });
                    poss.push(myth.variant.position);
                    refs.push(unsafe { String::from_utf8_unchecked(myth.variant.ref_seq.clone()) });
                    alts.push(unsafe { String::from_utf8_unchecked(myth.variant.alt_seq.clone()) });
                    source.push(unsafe { String::from_utf8_unchecked(annotation.source) });
                    transcript_id
                        .push(unsafe { String::from_utf8_unchecked(annotation.transcript_id) });
                    gene_name.push(unsafe { String::from_utf8_unchecked(annotation.gene_name) });
                    effects.push(
                        annotation
                            .effects
                            .iter()
                            .map(|e| unsafe { String::from_utf8_unchecked(e.clone().into()) })
                            .collect::<Vec<String>>()
                            .join(";"),
                    );
                    impact.push(annotation.impact as u8);
                }
            } else {
                end = false;
                break;
            }
        }

        let batch = arrow::record_batch::RecordBatch::try_new(
            schema.clone(),
            vec![
                std::sync::Arc::new(arrow::array::StringArray::from(chrs)),
                std::sync::Arc::new(arrow::array::UInt64Array::from(poss)),
                std::sync::Arc::new(arrow::array::StringArray::from(refs)),
                std::sync::Arc::new(arrow::array::StringArray::from(alts)),
                std::sync::Arc::new(arrow::array::StringArray::from(source)),
                std::sync::Arc::new(arrow::array::StringArray::from(transcript_id)),
                std::sync::Arc::new(arrow::array::StringArray::from(gene_name)),
                std::sync::Arc::new(arrow::array::StringArray::from(effects)),
                std::sync::Arc::new(arrow::array::UInt8Array::from(impact)),
            ],
        )?;

        writer.write(&batch)?;
    }

    writer.close()?;
    Ok(())
}

#[cfg(not(feature = "parallel"))]
/// For each variants found matching annotations
pub fn annotate<R>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    mut vcf_reader: variant::VcfReader<R>,
    no_annotation: bool,
    block_size: usize,
    writer: &mut dyn MythWriter,
) -> error::Result<()>
where
    R: std::io::BufRead,
{
    let variant2myth =
        variant2myth::Variant2Myth::new(annotations, translate, sequences, no_annotation);

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
/// For each variants found matching annotations
pub fn vcf2json<R, W>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    vcf_reader: variant::VcfReader<R>,
    no_annotation: bool,
    mut output: W,
) -> error::Result<()>
where
    R: std::io::BufRead + std::marker::Send,
    W: std::io::Write + std::io::Seek + std::marker::Send + 'static,
{
    let (tx, rx) = std::sync::mpsc::channel::<myth::Myth>();

    let write_thread = std::thread::spawn(move || {
        let schema = std::sync::Arc::new(output::schema());
        let mut writer = parquet::arrow::arrow_writer::ArrowWriter::try_new(
            &mut output,
            schema.clone(),
            Default::default(),
        )
        .unwrap();

        let block_size = 1 << 13;

        let mut end = true;
        while end {
            let mut chrs = Vec::with_capacity(block_size);
            let mut poss = Vec::with_capacity(block_size);
            let mut refs = Vec::with_capacity(block_size);
            let mut alts = Vec::with_capacity(block_size);
            let mut source = Vec::with_capacity(block_size);
            let mut transcript_id = Vec::with_capacity(block_size);
            let mut gene_name = Vec::with_capacity(block_size);
            let mut effects = Vec::with_capacity(block_size);
            let mut impact: Vec<u8> = Vec::with_capacity(block_size);

            for _ in 0..block_size {
                if let Some(myth) = rx.iter().next() {
                    for annotation in myth.annotations {
                        chrs.push(unsafe {
                            String::from_utf8_unchecked(myth.variant.seqname.clone())
                        });
                        poss.push(myth.variant.position);
                        refs.push(unsafe {
                            String::from_utf8_unchecked(myth.variant.ref_seq.clone())
                        });
                        alts.push(unsafe {
                            String::from_utf8_unchecked(myth.variant.alt_seq.clone())
                        });
                        source.push(unsafe { String::from_utf8_unchecked(annotation.source) });
                        transcript_id
                            .push(unsafe { String::from_utf8_unchecked(annotation.transcript_id) });
                        gene_name
                            .push(unsafe { String::from_utf8_unchecked(annotation.gene_name) });
                        effects.push(
                            annotation
                                .effects
                                .iter()
                                .map(|e| unsafe { String::from_utf8_unchecked(e.clone().into()) })
                                .collect::<Vec<String>>()
                                .join(";"),
                        );
                        impact.push(annotation.impact as u8);
                    }
                } else {
                    end = false;
                    break;
                }
            }

            let batch = arrow::record_batch::RecordBatch::try_new(
                schema.clone(),
                vec![
                    std::sync::Arc::new(arrow::array::StringArray::from(chrs)),
                    std::sync::Arc::new(arrow::array::UInt64Array::from(poss)),
                    std::sync::Arc::new(arrow::array::StringArray::from(refs)),
                    std::sync::Arc::new(arrow::array::StringArray::from(alts)),
                    std::sync::Arc::new(arrow::array::StringArray::from(source)),
                    std::sync::Arc::new(arrow::array::StringArray::from(transcript_id)),
                    std::sync::Arc::new(arrow::array::StringArray::from(gene_name)),
                    std::sync::Arc::new(arrow::array::StringArray::from(effects)),
                    std::sync::Arc::new(arrow::array::UInt8Array::from(impact)),
                ],
            )
            .unwrap();

            writer.write(&batch).unwrap();
        }
    });

    let variant2myth =
        variant2myth::Variant2Myth::new(annotations, translate, sequences, no_annotation);

    let results = vcf_reader
        .par_bridge()
        .filter(Result::is_ok)
        .map(error::Result::unwrap)
        .map(|variant| tx.send(variant2myth.myth(variant)))
        .collect::<Vec<core::result::Result<(), std::sync::mpsc::SendError<myth::Myth>>>>();

    for result in results {
        result?
    }

    drop(tx);

    write_thread.join().unwrap();

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
