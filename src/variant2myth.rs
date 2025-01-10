//! A collection of Annotator apply for each (variant, transcript)

/* std use */

/* crate use */
use bstr::ByteSlice as _;

/* module declaration */
mod feature_presence;
mod sequence_analysis;

/* project use */
use crate::annotation;
use crate::annotations_db;
use crate::effect;
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
        annotations: &[&annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<effect::Effect>;
}

/// Struct that associate to a variant myth
pub struct Variant2Myth<'a> {
    annotations: &'a annotations_db::AnnotationsDataBase,
    annotators: Vec<Box<dyn Annotator + std::marker::Send + std::marker::Sync + 'a>>,
}

impl<'a> Variant2Myth<'a> {
    /// Create Variant2Myth struct
    pub fn new(
        annotations: &'a annotations_db::AnnotationsDataBase,
        translate: &'a translate::Translate,
        sequences: &'a sequences_db::SequencesDataBase,
        annotators_choices: AnnotatorsChoices,
    ) -> Self {
        let mut annotators: Vec<Box<dyn Annotator + std::marker::Send + std::marker::Sync>> =
            Vec::new();

        if annotators_choices.contains(AnnotatorsChoicesRaw::Feature) {
            annotators.push(Box::new(feature_presence::FeaturePresence::new(
                b"upstream",
                effect::Effect::UpstreamGeneVariant,
            )));

            annotators.push(Box::new(feature_presence::FeaturePresence::new(
                b"downstream",
                effect::Effect::DownstreamGeneVariant,
            ))
                as Box<dyn Annotator + std::marker::Send + std::marker::Sync>);
            annotators.push(Box::new(feature_presence::FeaturePresence::new(
                b"5UTR",
                effect::Effect::FivePrimeUtrVariant,
            ))
                as Box<dyn Annotator + std::marker::Send + std::marker::Sync>);
            annotators.push(Box::new(feature_presence::FeaturePresence::new(
                b"3UTR",
                effect::Effect::ThreePrimeUtrVariant,
            ))
                as Box<dyn Annotator + std::marker::Send + std::marker::Sync>);
            annotators.push(Box::new(sequence_analysis::SequenceAnalysis::new(
                translate, sequences,
            ))
                as Box<dyn Annotator + std::marker::Send + std::marker::Sync>);
        }

        Self {
            annotations,
            annotators,
        }
    }

    /// Generate myth associate to variant
    pub fn myth(&self, variant: variant::Variant) -> myth::Myth {
        let mut myth = myth::Myth::from_variant(variant.clone());

        if variant.alt_seq.contains(&b'*') {
            return myth;
        }

        let annotations = self
            .annotations
            .get_annotation(&variant.seqname, variant.get_interval());

        if annotations.is_empty() {
            myth.add_annotation(
                myth::AnnotationMyth::builder()
                    .source(b"variant_myth".to_vec())
                    .transcript_id(vec![])
                    .effects(vec![effect::Effect::IntergenicRegion])
                    .build()
                    .unwrap(), // No possible error in build
            );

            return myth;
        }

        let mut transcript2annotations: ahash::AHashMap<
            (Vec<u8>, Vec<u8>),
            Vec<&annotation::Annotation>,
        > = ahash::AHashMap::new();

        for transcript in annotations
            .iter()
            .filter(|&&x| x.get_feature() == b"transcript")
        {
            for annotation in self
                .annotations
                .get_annotation(&variant.seqname, transcript.get_interval())
            {
                let key = (
                    annotation.get_source().to_vec(),
                    annotation.get_parent().to_vec(),
                );

                transcript2annotations
                    .entry(key)
                    .and_modify(|x| x.push(annotation))
                    .or_insert(vec![annotation]);
            }
        }

        for (key, annotations) in transcript2annotations {
            let mut transcript_myth = myth::AnnotationMyth::builder()
                .source(key.0)
                .transcript_id(key.1)
                .effects(vec![]);

            let gene_name = annotations
                .iter()
                .filter(|a| a.get_feature().contains_str("gene"))
                .map(|a| a.get_attribute().get_id())
                .collect::<Vec<&[u8]>>()
                .join(&b';');

            transcript_myth = transcript_myth.gene_name(gene_name);

            for annotator in &self.annotators[..] {
                transcript_myth.extend_effect(&annotator.annotate(&annotations, &variant))
            }

            myth.add_annotation(transcript_myth.build().unwrap());
        }

        myth
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */

    use crate::variant2myth::AnnotatorsChoicesRaw;

    use super::AnnotatorsChoices;

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
}
