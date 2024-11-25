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
        no_annotators: bool,
    ) -> Self {
        let mut annotators: Vec<Box<dyn Annotator + std::marker::Send + std::marker::Sync>> = vec![
            Box::new(feature_presence::FeaturePresence::new(
                b"upstream",
                effect::Effect::UpstreamGeneVariant,
            )),
            Box::new(feature_presence::FeaturePresence::new(
                b"downstream",
                effect::Effect::DownstreamGeneVariant,
            )),
            Box::new(feature_presence::FeaturePresence::new(
                b"5UTR",
                effect::Effect::FivePrimeUtrVariant,
            )),
            Box::new(feature_presence::FeaturePresence::new(
                b"3UTR",
                effect::Effect::ThreePrimeUtrVariant,
            )),
            Box::new(sequence_analysis::SequenceAnalysis::new(
                translate, sequences,
            )),
        ];

        if no_annotators {
            annotators.clear()
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

        for annotation in annotations {
            let key = (
                annotation.get_source().to_vec(),
                annotation.get_parent().to_vec(),
            );

            transcript2annotations
                .entry(key)
                .and_modify(|x| x.push(annotation))
                .or_insert(vec![annotation]);
        }

        for (key, annotations) in transcript2annotations {
            let mut transcript_myth = myth::AnnotationMyth::builder()
                .source(key.0)
                .transcript_id(key.1)
                .effects(vec![]);

            transcript_myth = transcript_myth.gene_name(
                annotations
                    .iter()
                    .filter(|a| a.get_feature().contains_str("gene"))
                    .map(|a| a.get_attribute().get_name())
                    .collect::<Vec<&[u8]>>()
                    .join(&b';'),
            );

            for annotator in &self.annotators[..] {
                transcript_myth.extend_effect(&annotator.annotate(&annotations, &variant))
            }

            myth.add_annotation(transcript_myth.build().unwrap());
        }

        myth
    }
}
