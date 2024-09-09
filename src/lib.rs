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
pub mod exons2sequences;
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

/// Get myth about variants
pub fn variants2myth(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    variant: variant::Variant,
) -> myth::Myth {
    let seqname = variant.seqname.to_vec();
    let interval = variant.get_interval();

    log::debug!("Start Variant: {:?}", variant);

    let mut myth = myth::Myth::from_variant(variant.clone());

    if variant.alt_seq.contains(&b'*') {
        return myth;
    }

    let var_annot = annotations.get_annotation(&seqname, interval);

    if var_annot.is_empty() {
        myth.add_annotation(
            myth::AnnotationMyth::builder()
                .source(b"variant_myth".to_vec())
                .transcript_id(vec![])
                .effects(vec![effect::Effect::IntergenicRegion])
                .build()
                .unwrap(), // No possible error in build
        )
    }

    for up_annot in var_annot.iter().filter(|a| a.get_feature() == b"upstream") {
        myth.add_annotation(
            myth::AnnotationMyth::builder()
                .source(up_annot.get_source().to_vec())
                .transcript_id(up_annot.get_transcript_id().to_vec())
                .effects(vec![effect::Effect::UpstreamGeneVariant])
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
                .effects(vec![effect::Effect::DownstreamGeneVariant])
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
            .gene_name(transcript.get_attribute().get_gene_name().to_vec())
            .effects(vec![]);

        let transcript_annot = annotations
            .get_annotation(transcript.get_seqname(), transcript.get_interval())
            .iter()
            .filter(|a| {
                a.get_attribute().get_transcript_id()
                    == transcript.get_attribute().get_transcript_id()
            })
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        let utr5_annot = transcript_annot
            .iter()
            .filter(|a| a.get_feature() == b"5UTR")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        let start_position = transcript_annot
            .iter()
            .find(|&&x| x.get_feature() == b"start_codon")
            .map(|x| x.get_start());

        let exon_annot = transcript_annot
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        let stop_position = transcript_annot
            .iter()
            .find(|&&x| x.get_feature() == b"stop_codon")
            .map(|x| x.get_stop());

        let utr3_annot = transcript_annot
            .iter()
            .filter(|a| a.get_feature() == b"3UTR")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // 5' UTR
        if !utr5_annot.is_empty() {
            transcript_myth.add_effect(effect::Effect::FivePrimeUtrVariant)
        }

        // Exon effect
        if !exon_annot.is_empty() {
            transcript_myth.extend_effect(
                exon_effect(
                    translate,
                    sequences,
                    exon_annot.as_ref(),
                    &variant,
                    start_position,
                    stop_position,
                )
                .as_ref(),
            );
        }

        // 3' UTR
        if !utr3_annot.is_empty() {
            transcript_myth.add_effect(effect::Effect::ThreePrimeUtrVariant)
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
    start_position: std::option::Option<u64>,
    stop_position: std::option::Option<u64>,
) -> Vec<effect::Effect> {
    let mut effects = vec![];

    let mut pos_exons = exon_annot
        .iter()
        .map(|a| {
            (
                a.get_attribute().get_exon_number(),
                (a.get_interval(), *a.get_strand()),
            )
        })
        .collect::<Vec<(u64, (core::ops::Range<u64>, annotation::Strand))>>();

    #[cfg(feature = "parallel")]
    pos_exons.par_sort_unstable_by_key(|a| a.0);
    #[cfg(not(feature = "parallel"))]
    pos_exons.sort_unstable_by_key(|a| a.0);

    let exons = pos_exons
        .iter()
        .map(|a| a.1.clone())
        .collect::<Vec<(core::ops::Range<u64>, annotation::Strand)>>();

    // Get transcript sequence
    let original_seq = sequences.get_transcript(exon_annot[0].get_seqname(), &exons);

    // Found position of variant
    let exon_target = if let Some(e) = exons
        .iter()
        .position(|e| variant.position > e.0.start && variant.position < e.0.end)
    {
        e
    } else {
        // Variant isn't present in a an exon end of investigation
        return effects;
    };

    // Check if variant is on a splice junction
    let variant_range = variant.position
        ..(variant.position + variant.ref_seq.len().max(variant.alt_seq.len()) as u64);

    if exon_target != 0
        && (exons[exon_target].0.start..exons[exon_target].0.start + 2)
            .any(|x| variant_range.contains(&x))
    {
        effects.push(effect::Effect::SpliceAcceptorVariant)
    }
    if (exons[exon_target].0.end - 2..exons[exon_target].0.end).any(|x| variant_range.contains(&x))
    {
        effects.push(effect::Effect::SpliceDonorVariant)
    }

    // Check if variant produce frameshift
    let len_diff = variant.ref_seq.len().abs_diff(variant.alt_seq.len());
    if len_diff != 0 {
        if len_diff % 3 != 0 {
            effects.push(effect::Effect::FrameshiftVariant)
        } else if variant.ref_seq.len() < variant.alt_seq.len() {
            effects.push(effect::Effect::ConservativeInframeInsertion)
        } else {
            effects.push(effect::Effect::ConservativeInframeDeletion)
        }
    }

    // Get sequence
    let transcript_pos = exons[..exon_target]
        .iter()
        .map(|e: &(core::ops::Range<u64>, annotation::Strand)| e.0.end - e.0.start)
        .sum::<u64>()
        + variant.position
        - exons[exon_target].0.start;

    let after_variant = original_seq
        .len()
        .min(transcript_pos as usize + variant.ref_seq.len() + 1);

    let variant_seq = original_seq[..=transcript_pos as usize]
        .iter()
        .chain(&variant.alt_seq)
        .chain(&original_seq[after_variant..])
        .cloned()
        .collect::<Vec<u8>>();

    println!();
    dbg!(exon_annot);
    dbg!(variant);
    dbg!(start_position, stop_position);
    let start = if let Some(s) = start_position {
        if let Some(exon_index) = exons.iter().position(|e| s > e.0.start && s < e.0.end) {
            if exons[exon_index].1 == annotation::Strand::Reverse {
                return effects;
            }

            dbg!(exon_index);
            dbg!(&exons[exon_index]);
            if variant.position < s {
                dbg!("start upper than variant position");
                (s - exons[exon_index].0.start) as usize + len_diff
            } else {
                dbg!("start lower than variant position");
                (s - exons[exon_index].0.start) as usize
            }
        } else {
            dbg!("start not found exon_index");
            0
        }
    } else {
        dbg!("no start_position");
        0
    };

    let stop = if let Some(s) = stop_position {
        if let Some(exon_index) = exons.iter().position(|e| s > e.0.start && s < e.0.end) {
            if exons[exon_index].1 == annotation::Strand::Reverse {
                return effects;
            }

            dbg!(exon_index);
            dbg!(&exons[exon_index]);
            if variant.position < s {
                dbg!("stop lower than variant position");
                variant_seq.len() - (exons[exon_index].0.end - s) as usize + len_diff
            } else {
                dbg!("stop upper than variant position");
                variant_seq.len() - (exons[exon_index].0.end - s) as usize
            }
        } else {
            dbg!("stop not found exon_index");
            variant_seq.len()
        }
    } else {
        dbg!("no stop_position");
        variant_seq.len()
    };

    dbg!(variant_seq.len(), start, stop);
    let var_coding_seq = &variant_seq[start..stop];

    let aa = translate.translate(var_coding_seq);

    effects
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

    const GFF: &[u8] = b"{0}	knownGene	transcript	11869	14409	.	+	.	gene_id=gene1;transcript_id=gene1
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

    fn setup() -> error::Result<(
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

        return Ok((
            translate,
            seq_db,
            annotations_db,
            fasta_reader[1..11].to_vec(),
        ));
    }

    #[test]
    fn basic_variant_annotation() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 29654,
            ref_seq: b"c".to_vec(),
            alt_seq: b"a".to_vec(),
        };

        let effects = exon_effect(&translate, &seq_db, &exon_annot, &variant, None, None);

        assert_eq!(effects, vec![]);

        Ok(())
    }

    #[test]
    fn splice_acceptor_annotation() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 30565,
            ref_seq: b"c".to_vec(),
            alt_seq: b"a".to_vec(),
        };

        let effects = exon_effect(&translate, &seq_db, &exon_annot, &variant, None, None);

        assert_eq!(effects, vec![effect::Effect::SpliceAcceptorVariant,]);

        Ok(())
    }

    #[test]
    fn splice_donor_annotation() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 30665,
            ref_seq: b"c".to_vec(),
            alt_seq: b"a".to_vec(),
        };

        let effects = exon_effect(&translate, &seq_db, &exon_annot, &variant, None, None);

        assert_eq!(effects, vec![effect::Effect::SpliceDonorVariant,]);

        Ok(())
    }

    #[test]
    fn framshift_annotation() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 29654,
            ref_seq: b"c".to_vec(),
            alt_seq: b"aa".to_vec(),
        };

        let effects = exon_effect(&translate, &seq_db, &exon_annot, &variant, None, None);

        assert_eq!(effects, vec![effect::Effect::FrameshiftVariant,]);

        Ok(())
    }

    #[test]
    fn codon_insertion_annotation() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 29654,
            ref_seq: b"c".to_vec(),
            alt_seq: b"caaa".to_vec(),
        };

        let effects = exon_effect(&translate, &seq_db, &exon_annot, &variant, None, None);

        assert_eq!(effects, vec![effect::Effect::ConservativeInframeInsertion,]);

        Ok(())
    }

    #[test]
    fn codon_deletion_annotation() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 29654,
            ref_seq: b"caaa".to_vec(),
            alt_seq: b"c".to_vec(),
        };

        let effects = exon_effect(&translate, &seq_db, &exon_annot, &variant, None, None);

        assert_eq!(effects, vec![effect::Effect::ConservativeInframeDeletion,]);

        Ok(())
    }

    #[test]
    fn not_a_variant_annotation() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 30100,
            ref_seq: b"c".to_vec(),
            alt_seq: b"c".to_vec(),
        };

        let effects = exon_effect(&translate, &seq_db, &exon_annot, &variant, None, None);

        assert_eq!(effects, vec![]);

        Ok(())
    }

    #[test]
    fn before_start() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 29560,
            ref_seq: b"c".to_vec(),
            alt_seq: b"c".to_vec(),
        };

        let effects = exon_effect(
            &translate,
            &seq_db,
            &exon_annot,
            &variant,
            Some(29565),
            Some(31042),
        );

        assert_eq!(effects, vec![]);

        Ok(())
    }

    #[test]
    fn after_start() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 29570,
            ref_seq: b"c".to_vec(),
            alt_seq: b"c".to_vec(),
        };

        let effects = exon_effect(
            &translate,
            &seq_db,
            &exon_annot,
            &variant,
            Some(29565),
            Some(31042),
        );

        assert_eq!(effects, vec![]);

        Ok(())
    }

    #[test]
    fn before_stop() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 31020,
            ref_seq: b"c".to_vec(),
            alt_seq: b"c".to_vec(),
        };

        let effects = exon_effect(
            &translate,
            &seq_db,
            &exon_annot,
            &variant,
            Some(29565),
            Some(31042),
        );

        assert_eq!(effects, vec![]);

        Ok(())
    }

    #[test]
    fn after_stop() -> error::Result<()> {
        let (translate, seq_db, annotations_db, contig_name) = setup()?;

        let exon_annot = annotations_db
            .get_annotation(&contig_name, 29554..31097)
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        // target variant
        let variant = variant::Variant {
            seqname: contig_name,
            position: 31020,
            ref_seq: b"c".to_vec(),
            alt_seq: b"c".to_vec(),
        };

        let effects = exon_effect(
            &translate,
            &seq_db,
            &exon_annot,
            &variant,
            Some(29565),
            Some(31012),
        );

        assert_eq!(effects, vec![]);

        Ok(())
    }
}
