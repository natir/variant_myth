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
pub mod variant;

#[cfg(not(feature = "parallel"))]
/// For each variants found matching annotations
pub fn vcf2json<R, W>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    vcf_reader: variant::VcfReader<R>,
    mut output: W,
) -> error::Result<()>
where
    R: std::io::BufRead,
    W: std::io::Write,
{
    for result in vcf_reader {
        serde_json::to_writer(&mut output, &variants2myth(annotations, sequences, result?))?;
    }

    Ok(())
}

#[cfg(feature = "parallel")]
/// For each variants found matching annotations
pub fn vcf2json<R, W>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    mut vcf_reader: variant::VcfReader<R>,
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
        .map(|variant| tx.send(variants2myth(annotations, sequences, variant)))
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

pub fn variants2myth(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    variant: variant::Variant,
) -> myth::Myth {
    let (seqname, interval) = variant.get_interval();

    log::debug!("Variant: {:?}", variant);

    let mut myth = myth::Myth::from_variant(variant.clone());

    for annotation in annotations.get_annotation(seqname, interval) {
        if annotation.get_feature() == b"transcript" {
            let annotation_myth = myth::AnnotationMyth::builder()
                .source(annotation.get_source().to_vec())
                .transcript_id(annotation.get_transcript_id().to_vec());

            let transcript_annot_interval = annotation.get_interval();
            for subannot in
                annotations.get_annotation(transcript_annot_interval.0, transcript_annot_interval.1)
            {
                if subannot.get_attribute().get_transcript_id()
                    == annotation.get_attribute().get_transcript_id()
                {
                    log::debug!("\tannotation: {}", subannot)
                }
            }

            myth.add_annotation(annotation_myth.build());
        }
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
