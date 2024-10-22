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
pub mod sequences_db;
pub mod translate;
pub mod variant;
pub mod variant2myth;

#[cfg(not(feature = "parallel"))]
/// For each variants found matching annotations
pub fn vcf2json<R, W>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    vcf_reader: variant::VcfReader<R>,
    mut output: W,
) -> error::Result<()>
where
    R: std::io::BufRead,
    W: std::io::Write,
{
    let variant2myth = variant2myth::Variant2Myth::new(annotations, translate, sequences);

    for result in vcf_reader {
        serde_json::to_writer(&mut output, &variant2myth.myth(result?))?;
    }

    Ok(())
}

#[cfg(feature = "parallel")]
/// For each variants found matching annotations
pub fn vcf2json<R, W>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    vcf_reader: variant::VcfReader<R>,
    mut output: W,
) -> error::Result<()>
where
    R: std::io::BufRead + std::marker::Send,
    W: std::io::Write + std::marker::Send + 'static,
{
    let (tx, rx) = std::sync::mpsc::channel();

    let write_thread = std::thread::spawn(move || {
        rx.iter()
            .map(|message| serde_json::to_writer(&mut output, &message))
            .collect::<Vec<serde_json::Result<()>>>()
    });

    let variant2myth = variant2myth::Variant2Myth::new(annotations, translate, sequences);

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

    let results = write_thread.join().unwrap();
    for result in results {
        result?
    }

    Ok(())
}

fn serialize_bstr<T, S>(v: T, serializer: S) -> Result<S::Ok, S::Error>
where
    T: AsRef<[u8]>,
    S: serde::Serializer,
{
    serializer.serialize_str(unsafe { std::str::from_utf8_unchecked(v.as_ref()) })
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */
    use biotest::Format as _;
    use bstr::ByteSlice as _;

    /* project use */
    use super::*;

    const _GFF: &[u8] =
        b"{0}	knownGene	transcript	11869	14409	.	+	.	gene_id=gene1;transcript_id=gene1
{0}	knownGene	exon	11869	12227	.	+	.	gene_id=gene1;transcript_id=gene1;exon_number=1;exon_id=gene1.1
{0}	knownGene	exon	12613	12721	.	+	.	gene_id=gene1;transcript_id=gene1;exon_number=2;exon_id=gene1.2
{0}	knownGene	exon	13221	14409	.	+	.	gene_id=gene1;transcript_id=gene1;exon_number=3;exon_id=gene1.3
{0}	knownGene	transcript	17369	17436	.	-	.	gene_id=gene2;transcript_id=gene2
{0}	knownGene	exon	17369	17436	.	-	.	gene_id=gene2;transcript_id=gene2;exon_number=1;exon_id=gene2.1
{0}	knownGene	transcript	29544	31107	.	+	.	gene_id=gene3;transcript_id=gene3
{0}	knownGene	exon	29554	30039	.	+	.	gene_id=gene3;transcript_id=gene3;exon_number=1;exon_id=gene3.1
{0}	knownGene	exon	30564	30667	.	+	.	gene_id=gene3;transcript_id=gene3;exon_number=2;exon_id=gene3.2
{0}	knownGene	exon	30976	31097	.	+	.	gene_id=gene3;transcript_id=gene3;exon_number=3;exon_id=gene3.3
";

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
                _GFF.replace(b"{0}", &fasta_reader[1..11]),
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
