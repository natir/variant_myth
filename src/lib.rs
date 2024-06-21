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

    log::debug!("Start Variant: {:?}", variant);

    let mut myth = myth::Myth::from_variant(variant.clone());

    if variant.alt_seq.contains(&b'*') {
        return myth;
    }

    let var_annot = annotations.get_annotation(seqname, interval);

    if var_annot.is_empty() {
        myth.add_annotation(
            myth::AnnotationMyth::builder()
                .source(b"variant_myth".to_vec())
                .transcript_id(vec![])
                .effects(vec![myth::Effect::Intergenic])
                .build()
                .unwrap(), // No possible error in build
        )
    }

    for up_annot in var_annot.iter().filter(|a| a.get_feature() == b"upstream") {
        myth.add_annotation(
            myth::AnnotationMyth::builder()
                .source(up_annot.get_source().to_vec())
                .transcript_id(up_annot.get_transcript_id().to_vec())
                .effects(vec![myth::Effect::Upstream])
                .build()
                .unwrap(), // No possible error in build
        )
    }

    for down_annot in var_annot
        .iter()
        .filter(|a| a.get_feature() == b"downstream")
    {
        myth.add_annotation(
            myth::AnnotationMyth::builder()
                .source(down_annot.get_source().to_vec())
                .transcript_id(down_annot.get_transcript_id().to_vec())
                .effects(vec![myth::Effect::Downstream])
                .build()
                .unwrap(), // No possible error in build
        )
    }

    for transcript in var_annot
        .iter()
        .filter(|a| a.get_feature() == b"transcript")
    {
        let mut transcript_myth = myth::AnnotationMyth::builder()
            .source(transcript.get_source().to_vec())
            .transcript_id(transcript.get_transcript_id().to_vec())
            .effects(vec![myth::Effect::Transcript]);

        let transcript_annot = annotations
            .get_annotation(transcript.get_seqname(), transcript.get_interval())
            .iter()
            .filter(|a| {
                a.get_attribute().get_transcript_id()
                    == transcript.get_attribute().get_transcript_id()
            })
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        if transcript_annot.iter().any(|a| a.get_feature() == b"CDS") {
            transcript_myth.add_effect(myth::Effect::Cds)
        }

        let utr5_annot = transcript_annot
            .iter()
            .filter(|a| a.get_feature() == b"5UTR")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        let exon_annot = transcript_annot
            .iter()
            .filter(|a| a.get_feature() == b"5UTR")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        let utr3_annot = transcript_annot
            .iter()
            .filter(|a| a.get_feature() == b"3UTR")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // 5' UTR
        if !utr5_annot.is_empty() {
            transcript_myth.add_effect(myth::Effect::Utr5Prime)
        }

        // Get exon sequence
        if !exon_annot.is_empty() {
            transcript_myth.extend_effect(
                exon_effect(translate, sequences, exon_annot.as_ref(), &variant).as_ref(),
            );
        }

        // 3' UTR
        if !utr3_annot.is_empty() {
            transcript_myth.add_effect(myth::Effect::Utr5Prime)
        }

        myth.add_annotation(transcript_myth.build().unwrap()); // Build error could never append
    }

    log::debug!("End Variant: {:?}", variant);
    myth
}

fn exon_effect(
    translate: &translate::Translate,
    sequences: &sequences_db::SequencesDataBase,
    exon_annot: &[&annotation::Annotation],
    variant: &variant::Variant,
) -> Vec<myth::Effect> {
    let mut pos_exons = exon_annot
        .iter()
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
    let original_seq = sequences.get_transcript(exon_annot[0].get_seqname(), &exons);

    // Found position
    let options = exons
        .iter()
        .position(|e| variant.position > e.0.start && variant.position < e.0.end);
    if options.is_none() {
        return vec![];
    }
    let exon_target = options.unwrap();

    let transcript_pos = exons[..exon_target]
        .iter()
        .map(|e: &(interval_tree::Interval<u64>, annotation::Strand)| e.0.end - e.0.start)
        .sum::<u64>()
        + variant.position
        - exons[exon_target].0.start;

    let after_variant = original_seq
        .len()
        .min(transcript_pos as usize + variant.ref_seq.len());

    let variant_seq = original_seq[..transcript_pos as usize]
        .iter()
        .chain(&variant.alt_seq)
        .chain(&original_seq[after_variant..])
        .cloned()
        .collect::<Vec<u8>>();

    let _aa = translate.translate(&variant_seq);

    vec![]
}

fn serialize_bstr<T, S>(v: T, serializer: S) -> Result<S::Ok, S::Error>
where
    T: AsRef<[u8]>,
    S: serde::Serializer,
{
    serializer.serialize_str(unsafe { std::str::from_utf8_unchecked(v.as_ref()) })
}
