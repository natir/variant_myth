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
pub mod error;
pub mod interval_tree;
pub mod myth;
pub mod sequences_db;
pub mod translate;
pub mod variant;

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
    for result in vcf_reader {
        serde_json::to_writer(
            &mut output,
            &variants2myth(annotations, sequences, translate, result?),
        )?;
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

    let results = vcf_reader
        .par_bridge()
        .filter(Result::is_ok)
        .map(error::Result::unwrap)
        .map(|variant| tx.send(variants2myth(annotations, sequences, translate, variant)))
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

/// Get myth about variants
pub fn variants2myth(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    variant: variant::Variant,
) -> myth::Myth {
    let (seqname, interval) = variant.get_interval();

    log::debug!("Variant: {:?}", variant);

    let mut myth = myth::Myth::from_variant(variant.clone());

    for transcript in annotations
        .get_annotation(seqname, interval)
        .iter()
        .filter(|a| a.get_feature() == b"transcript")
    {
        let mut annotation_myth = myth::AnnotationMyth::builder()
            .source(transcript.get_source().to_vec())
            .transcript_id(transcript.get_transcript_id().to_vec());

        let related_annot = annotations
            .get_annotation(transcript.get_seqname(), transcript.get_interval())
            .into_iter()
            .filter(|a| {
                a.get_attribute().get_transcript_id()
                    == transcript.get_attribute().get_transcript_id()
            })
            .collect::<Vec<&annotation::Annotation>>();

        // 5' UTR

        // Get exon sequence
        let mut pos_exons = related_annot
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .map(|a| {
                (
                    a.get_attribute().get_exon_number(),
                    (a.get_interval(), *a.get_strand()),
                )
            })
            .collect::<Vec<(u64, (interval_tree::Interval<u64>, annotation::Strand))>>();

        #[cfg(feature = "parallel")]
        pos_exons.par_sort_unstable_by_key(|a| a.0);
        #[cfg(not(feature = "parallel"))]
        pos_exons.sort_unstable_by_key(|a| a.0);

        let exons = pos_exons
            .iter()
            .map(|a| a.1.clone())
            .collect::<Vec<(interval_tree::Interval<u64>, annotation::Strand)>>();

        // Edit sequence with variant
        let sequence = sequences.get_transcript(transcript.get_seqname(), &exons);

        // Found position

        let aa = translate.translate(&sequence);

        // 3' UTR

        myth.add_annotation(annotation_myth.build().unwrap()); // Build error could never append
    }

    myth
}

fn serialize_bstr<T, S>(v: T, serializer: S) -> Result<S::Ok, S::Error>
where
    T: AsRef<[u8]>,
    S: serde::Serializer,
{
    serializer.serialize_str(unsafe { std::str::from_utf8_unchecked(v.as_ref()) })
}
