//! A collection of Annotator apply for each (variant, transcript)

/* std use */

/* crate use */

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
        annotations: &[annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<effect::Effect>;
}

/// Struct that associate to a variant myth
pub struct Variant2Myth<'a> {
    annotations: &'a mut annotations_db::AnnotationsDataBase,
    annotators: Vec<Box<dyn Annotator + std::marker::Send + std::marker::Sync + 'a>>,
    annotators_choices: AnnotatorsChoices,
}

impl<'a> Variant2Myth<'a> {
    /// Create Variant2Myth struct
    pub fn new(
        annotations: &'a mut annotations_db::AnnotationsDataBase,
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
        }

        if annotators_choices.contains(AnnotatorsChoicesRaw::Effect) {
            annotators.push(Box::new(sequence_analysis::SequenceAnalysis::new(
                translate, sequences,
            ))
                as Box<dyn Annotator + std::marker::Send + std::marker::Sync>);
        }

        Self {
            annotations,
            annotators,
            annotators_choices,
        }
    }

    /// Generate myth associate to variant
    pub fn myth(&mut self, variant: variant::Variant) -> myth::Myth {
        let mut myth = myth::Myth::from_variant(variant.clone());

        if variant.alt_seq.contains(&b'*') {
            return myth;
        }

        let annotations = self
            .annotations
            .get_annotations(&variant.seqname, variant.get_interval());

        if annotations.is_empty() {
            myth.add_annotation(
                myth::AnnotationMyth::builder()
                    .source(b"variant_myth".to_vec())
                    .feature(b"unknow".to_vec())
                    .id(b"".to_vec())
                    .name(b"".to_vec())
                    .effects(vec![effect::Effect::IntergenicRegion])
                    .build()
                    .unwrap(), // No possible error in build
            );

            return myth;
        }

        if self.annotators_choices.contains(AnnotatorsChoicesRaw::Gene) {
            for gene in annotations.iter().filter(|x| x.get_feature() == b"gene") {
                let gene_myth = myth::AnnotationMyth::builder()
                    .source(gene.get_source().to_vec())
                    .feature(gene.get_feature().to_vec())
                    .id(gene.get_attribute().get_id().to_vec())
                    .name(gene.get_attribute().get_name().to_vec())
                    .effects(vec![]);

                myth.add_annotation(gene_myth.build().unwrap()); // build can't failled
            }
        }

        if self
            .annotators_choices
            .contains(AnnotatorsChoicesRaw::Feature)
            || self
                .annotators_choices
                .contains(AnnotatorsChoicesRaw::Effect)
        {
            for transcript in annotations
                .iter()
                .filter(|x| x.get_feature() == b"transcript")
            {
                let annotations = if let Some(a) = self
                    .annotations
                    .get_subannotations(transcript.get_attribute().get_id())
                {
                    a
                } else {
                    continue;
                };

                let mut transcript_myth = myth::AnnotationMyth::builder()
                    .source(transcript.get_source().to_vec())
                    .feature(transcript.get_feature().to_vec())
                    .id(transcript.get_attribute().get_id().to_vec())
                    .name(transcript.get_attribute().get_name().to_vec())
                    .effects(vec![]);

                for annotator in &self.annotators[..] {
                    transcript_myth.extend_effect(&annotator.annotate(annotations, &variant))
                }

                myth.add_annotation(transcript_myth.build().unwrap()); // build can't failled
            }
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
