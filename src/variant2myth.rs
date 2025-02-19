//! A collection of Annotator apply for each (variant, transcript)

/* std use */

/* crate use */

/* module declaration */
mod feature_presence;
mod sequence_analysis;

/* project use */
use crate::annotations_db;
use crate::effect;
use crate::memoizor;
use crate::myth;
use crate::sequences_db;
use crate::translate;
use crate::variant;

#[enumflags2::bitflags]
#[repr(u8)]
#[derive(std::clone::Clone, std::marker::Copy, std::fmt::Debug)]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
/// Choose how variant are annotate
pub enum AnnotatorsChoicesRaw {
    /// Variant are annotate with just gene information
    Gene = 1 << 1,
    /// Variant are annotate with feature
    Feature = 1 << 2,
    /// Variant are annotate with effect
    Effect = 1 << 3,
    /// Variant are annotate with Hgvs
    Hgvs = 1 << 4,
}

/// Choose how variant are annotate
pub type AnnotatorsChoices = enumflags2::BitFlags<AnnotatorsChoicesRaw>;

trait Annotator {
    fn annotate(
        &self,
        variant: &variant::Variant,
        memoizor: &mut memoizor::Memoizor,
    ) -> Vec<effect::Effect>;
}

/// Struct that associate to a variant myth
pub struct Variant2Myth<'a> {
    annotations: &'a annotations_db::AnnotationsDataBase,
    sequences: &'a sequences_db::SequencesDataBase,
    annotators: [Vec<Box<dyn Annotator + std::marker::Send + std::marker::Sync + 'a>>; 5],
    annotators_choices: AnnotatorsChoices,
}

impl<'a> Variant2Myth<'a> {
    /// Create Variant2Myth struct
    pub fn new(
        annotations: &'a annotations_db::AnnotationsDataBase,
        translate: &'a translate::Translate,
        sequences: &'a sequences_db::SequencesDataBase,
        annotators_choices: AnnotatorsChoices,
    ) -> Self {
        let mut annotators: [Vec<Box<dyn Annotator + std::marker::Send + std::marker::Sync>>; 5] =
            [Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()];

        annotators[(AnnotatorsChoicesRaw::Gene as u8).ilog2() as usize].extend([]);
        annotators[(AnnotatorsChoicesRaw::Feature as u8).ilog2() as usize].extend([
            Box::new(feature_presence::FeaturePresence::new(
                b"upstream",
                effect::Effect::UpstreamGeneVariant,
            )) as Box<dyn Annotator + Send + Sync>,
            Box::new(feature_presence::FeaturePresence::new(
                b"downstream",
                effect::Effect::DownstreamGeneVariant,
            )) as Box<dyn Annotator + Send + Sync>,
            Box::new(feature_presence::FeaturePresence::new(
                b"five_prime_UTR",
                effect::Effect::FivePrimeUtrVariant,
            )) as Box<dyn Annotator + Send + Sync>,
            Box::new(feature_presence::FeaturePresence::new(
                b"three_prime_UTR",
                effect::Effect::ThreePrimeUtrVariant,
            )) as Box<dyn Annotator + Send + Sync>,
        ]);
        annotators[(AnnotatorsChoicesRaw::Effect as u8).ilog2() as usize].push(Box::new(
            sequence_analysis::SequenceAnalysis::new(translate, sequences),
        ));
        annotators[(AnnotatorsChoicesRaw::Hgvs as u8).ilog2() as usize].extend([]);

        Self {
            annotations,
            sequences,
            annotators,
            annotators_choices,
        }
    }

    /// Generate myth associate to variant
    pub fn myth(&self, variant: variant::Variant) -> myth::Myth {
        let mut myth = myth::Myth::from_variant(variant.clone());

        // Ignore not variant we could manage
        if !variant.valid() {
            myth.add_annotation(
                myth::AnnotationMyth::from_nowhere()
                    .effects(vec![effect::Effect::Ignore])
                    .build()
                    .unwrap(), // No possible error in build
            );

            return myth;
        }

        // Get annotation
        let not_coding_annotations = self
            .annotations
            .get_annotations(&variant.seqname, variant.get_interval());

        // Intergenic region
        if not_coding_annotations.is_empty() {
            myth.add_annotation(
                myth::AnnotationMyth::from_nowhere()
                    .effects(vec![effect::Effect::IntergenicRegion])
                    .build()
                    .unwrap(), // No possible error in build
            );

            return myth;
        }

        // Add myth with gene associate information
        if self.annotators_choices.contains(AnnotatorsChoicesRaw::Gene) {
            for gene in not_coding_annotations
                .iter()
                .filter(|a| a.get_feature() == b"gene")
            {
                myth.add_annotation(
                    myth::AnnotationMyth::from_annotation(gene)
                        .effects(vec![])
                        .build()
                        .unwrap(), // No possible error in build
                );
            }
        }

        // Get unique transcript
        let mut transcripts = ahash::AHashSet::new();
        for annotation in not_coding_annotations.iter() {
            if annotation.get_feature() == b"transcript" {
                transcripts.insert(annotation.get_attribute().get_id())
            } else {
                transcripts.insert(annotation.get_attribute().get_parent())
            };
        }

        for transcript_id in transcripts.iter() {
            let mut memoizor = memoizor::Memoizor::new(
                transcript_id,
                self.annotations,
                self.sequences,
                &not_coding_annotations,
            );

            // found transcript annotations
            let mut annotation_myth = if let Some(transcript_annot) = memoizor.transcript() {
                myth::AnnotationMyth::from_annotation(transcript_annot).effects(vec![])
            } else {
                continue;
            };

            for flag in self.annotators_choices.iter() {
                self.annotators[(flag as u8).ilog2() as usize]
                    .iter()
                    .for_each(|a| {
                        annotation_myth.extend_effect(&a.annotate(&variant, &mut memoizor));
                    })
            }

            myth.add_annotation(annotation_myth.build().unwrap()) // No possible error in build
        }

        myth
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */
    use super::*;
    use crate::error;
    use crate::test_data;

    #[test]
    fn annotator_choices() {
        let mut flag = AnnotatorsChoices::empty();

        assert!(!flag.contains(AnnotatorsChoicesRaw::Effect));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Feature));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Gene));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Hgvs));

        flag |= AnnotatorsChoicesRaw::Effect;

        assert!(flag.contains(AnnotatorsChoicesRaw::Effect));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Feature));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Gene));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Hgvs));

        flag |= AnnotatorsChoicesRaw::Feature;

        assert!(flag.contains(AnnotatorsChoicesRaw::Effect));
        assert!(flag.contains(AnnotatorsChoicesRaw::Feature));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Gene));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Hgvs));

        flag |= AnnotatorsChoicesRaw::Gene;

        assert!(flag.contains(AnnotatorsChoicesRaw::Effect));
        assert!(flag.contains(AnnotatorsChoicesRaw::Feature));
        assert!(flag.contains(AnnotatorsChoicesRaw::Gene));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Hgvs));

        flag |= AnnotatorsChoicesRaw::Hgvs;

        assert!(flag.contains(AnnotatorsChoicesRaw::Effect));
        assert!(flag.contains(AnnotatorsChoicesRaw::Feature));
        assert!(flag.contains(AnnotatorsChoicesRaw::Gene));
        assert!(flag.contains(AnnotatorsChoicesRaw::Hgvs));

        flag = AnnotatorsChoicesRaw::Effect | AnnotatorsChoicesRaw::Hgvs;

        assert!(flag.contains(AnnotatorsChoicesRaw::Effect));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Feature));
        assert!(!flag.contains(AnnotatorsChoicesRaw::Gene));
        assert!(flag.contains(AnnotatorsChoicesRaw::Hgvs));
    }

    #[test]
    pub fn unvalid_variant() -> error::Result<()> {
        let reader: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>> =
            std::io::BufReader::new(Box::new(test_data::GFF));
        let annotations_db = annotations_db::AnnotationsDataBase::from_reader(reader, 100)?;

        let translate = translate::Translate::default();

        let reader: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>> =
            std::io::BufReader::new(Box::new(test_data::SEQUENCE));
        let sequences_db = sequences_db::SequencesDataBase::from_reader(reader)?;

        let variant2myth = Variant2Myth::new(
            &annotations_db,
            &translate,
            &sequences_db,
            AnnotatorsChoices::empty(),
        );

        let variant = variant::Variant::test_variant(b"chrA", 73_602, b"A", b"", None)?;
        let mut truth = myth::Myth::from_variant(variant.clone());
        truth.add_annotation(
            myth::AnnotationMyth::from_nowhere()
                .effects(vec![effect::Effect::Ignore])
                .build()
                .unwrap(),
        );

        assert_eq!(variant2myth.myth(variant), truth);

        Ok(())
    }
}
