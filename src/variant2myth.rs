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

trait Annotator {
    fn annotate(
        &self,
        annotations: &[&annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<effect::Effect>;
}

struct Variant2Myth<'a> {
    annotations: &'a annotations_db::AnnotationsDataBase,
    translate: &'a translate::Translate,
    sequences: &'a sequences_db::SequencesDataBase,
    annotators: [&'static dyn Annotator; 4],
}

static UPSTREAM: feature_presence::FeaturePresence =
    feature_presence::FeaturePresence::new(b"upstream", effect::Effect::UpstreamGeneVariant);
static DOWNSTREAM: feature_presence::FeaturePresence =
    feature_presence::FeaturePresence::new(b"downstream", effect::Effect::DownstreamGeneVariant);
static FIVE_PRIME_UTR: feature_presence::FeaturePresence =
    feature_presence::FeaturePresence::new(b"5UTR", effect::Effect::FivePrimeUtrVariant);
static THREE_PRIME_UTR: feature_presence::FeaturePresence =
    feature_presence::FeaturePresence::new(b"3UTR", effect::Effect::ThreePrimeUtrVariant);
static SEQ_ANALYSIS: sequence_analysis::SequenceAnalysis =
    sequence_analysis::SequenceAnalysis::new();

impl<'a> Variant2Myth<'a> {
    pub fn new(
        annotations: &'a annotations_db::AnnotationsDataBase,
        translate: &'a translate::Translate,
        sequences: &'a sequences_db::SequencesDataBase,
    ) -> Self {
        let annotators = [
            &UPSTREAM as &'static (dyn Annotator + 'static),
            &DOWNSTREAM as &'static (dyn Annotator + 'static),
            &FIVE_PRIME_UTR as &'static (dyn Annotator + 'static),
            &THREE_PRIME_UTR as &'static (dyn Annotator + 'static),
            //    &SEQ_ANALYSIS as &'static (dyn Annotator + 'static),
        ];

        Self {
            annotations,
            translate,
            sequences,
            annotators,
        }
    }

    pub fn myth(&mut self, variant: variant::Variant) -> myth::Myth {
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
                annotation.get_transcript_id().to_vec(),
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

            for annotator in &self.annotators[..] {
                transcript_myth.extend_effect(&annotator.annotate(&annotations, &variant))
            }

            myth.add_annotation(transcript_myth.build().unwrap());
        }

        myth
    }
}
