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
pub mod effect;
pub mod error;
pub mod myth;
pub mod output;
pub mod sequences_db;
pub mod translate;
pub mod variant;
pub mod variant2myth;

#[cfg(feature = "cli")]
pub mod cli;

#[cfg(feature = "json")]
fn serialize_bstr<T, S>(v: T, serializer: S) -> Result<S::Ok, S::Error>
where
    T: AsRef<[u8]>,
    S: serde::Serializer,
{
    serializer.serialize_str(unsafe { std::str::from_utf8_unchecked(v.as_ref()) })
}

/// For each variants found matching annotations
#[cfg(not(feature = "parallel"))]
pub fn vcf2myth<R>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    vcf_reader: variant::VcfReader<R>,
    annotators_choices: variant2myth::AnnotatorsChoices,
    mut writer: Box<dyn output::MythWriter>,
) -> error::Result<()>
where
    R: std::io::BufRead,
{
    let variant2myth =
        variant2myth::Variant2Myth::new(annotations, translate, sequences, annotators_choices);

    for result in vcf_reader {
        let variant = result?;

        let myth = variant2myth.myth(variant);

        writer.write_myth(myth)?;
    }

    writer.close()?;

    Ok(())
}

/// For each variants found matching annotations
#[cfg(feature = "parallel")]
pub fn vcf2myth<R>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    vcf_reader: variant::VcfReader<R>,
    annotators_choices: variant2myth::AnnotatorsChoices,
    mut writer: Box<dyn output::MythWriter + std::marker::Send>,
) -> error::Result<()>
where
    R: std::io::BufRead + std::marker::Send,
{
    let (tx, rx) = std::sync::mpsc::channel::<myth::Myth>();

    let write_thread = std::thread::spawn(move || -> error::Result<()> {
        for myth in rx {
            writer.write_myth(myth)?;
        }
        writer.close()?;
        Ok(())
    });

    let variant2myth =
        variant2myth::Variant2Myth::new(annotations, translate, sequences, annotators_choices);

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

    write_thread.join().unwrap() // Err only if panic! if write thread panic all should panic !
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
{0}	ensembl_havana	transcript	826	4439	.	+	.	ID=transcript:ENST00000374675;Parent=gene:ENSG00000112473;Name=SLC39A7-201
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
