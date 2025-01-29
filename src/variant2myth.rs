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
        annotations: &[&annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<effect::Effect>;
}

/// Struct that associate to a variant myth
pub struct Variant2Myth<'a> {
    annotations: &'a annotations_db::AnnotationsDataBase,
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
                b"5UTR",
                effect::Effect::FivePrimeUtrVariant,
            )) as Box<dyn Annotator + Send + Sync>,
            Box::new(feature_presence::FeaturePresence::new(
                b"3UTR",
                effect::Effect::ThreePrimeUtrVariant,
            )) as Box<dyn Annotator + Send + Sync>,
        ]);
        annotators[(AnnotatorsChoicesRaw::Effect as u8).ilog2() as usize].push(Box::new(
            sequence_analysis::SequenceAnalysis::new(translate, sequences),
        ));
        annotators[(AnnotatorsChoicesRaw::Hgvs as u8).ilog2() as usize].extend([]);

        Self {
            annotations,
            annotators,
            annotators_choices,
        }
    }

    /// Generate myth associate to variant
    pub fn myth(&self, variant: variant::Variant) -> myth::Myth {
        let mut myth = myth::Myth::from_variant(variant.clone());

        // Detect not normalize variant
        if variant.alt_seq.contains(&b'*') {
            return myth;
        }

        // Get annotation
        let not_coding_annotations = self
            .annotations
            .get_annotations(&variant.seqname, variant.get_interval());

        // Intergenic region
        if not_coding_annotations.is_empty() {
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

        // Group annotation by transcript
        let mut transcript2annotations = ahash::AHashMap::new();
        for annotation in not_coding_annotations {
            let entry = if annotation.get_feature() == b"transcript" {
                transcript2annotations.entry(annotation.get_attribute().get_id())
            } else {
                transcript2annotations.entry(annotation.get_attribute().get_parent())
            };
            entry
                .and_modify(|x: &mut Vec<&annotation::Annotation>| x.push(annotation))
                .or_insert(vec![annotation]);
        }

        for (transcript_name, annotations) in transcript2annotations.iter() {
            // found transcript
            let mut annotation_myth = if let Some(transcript_annot) = annotations
                .iter()
                .find(|x| x.get_attribute().get_id() == *transcript_name)
            {
                myth::AnnotationMyth::from_annotation(transcript_annot)
            } else {
                continue;
            };

            // Not intergenic region
            for flag in self.annotators_choices.iter() {
                match flag {
                    AnnotatorsChoicesRaw::Hgvs | AnnotatorsChoicesRaw::Effect => {
                        if let Some(coding_annotation) =
                            self.annotations.get_coding_annotation(transcript_name)
                        {
                            let proxy = coding_annotation
                                .iter()
                                .collect::<Vec<&annotation::Annotation>>(); // TODO: found a solution to remove this

                            self.annotators[flag as usize].iter().for_each(|a| {
                                annotation_myth.extend_effect(&a.annotate(&proxy, &variant))
                            });
                        } else {
                            continue;
                        }
                    }
                    _ => {
                        self.annotators[flag as usize].iter().for_each(|a| {
                            annotation_myth.extend_effect(&a.annotate(annotations, &variant))
                        });
                    }
                }
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
