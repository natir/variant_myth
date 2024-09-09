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
    ) -> Vec<myth::AnnotationMyth>;
}

struct Variant2Myth<'a> {
    annotations: &'a annotations_db::AnnotationsDataBase,
    translate: &'a translate::Translate,
    sequences: &'a sequences_db::SequencesDataBase,
    all_annotators: [&'static dyn Annotator; 4],
    transcript_annotators: [&'static dyn Annotator; 1],
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
        let all_annotators = [
            &UPSTREAM as &'static (dyn Annotator + 'static),
            &DOWNSTREAM as &'static (dyn Annotator + 'static),
            &FIVE_PRIME_UTR as &'static (dyn Annotator + 'static),
            &THREE_PRIME_UTR as &'static (dyn Annotator + 'static),
        ];

        let transcript_annotators = [&SEQ_ANALYSIS as &'static (dyn Annotator + 'static)];

        Self {
            annotations,
            translate,
            sequences,
            all_annotators,
            transcript_annotators,
        }
    }

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

        for annotator in &self.all_annotators[..] {
            myth.extend_annotation(&annotator.annotate(&annotations, &variant))
        }

        for transcript in annotations
            .iter()
            .filter(|a| a.get_feature() == b"transcript")
        {
            let transcript_annotations = self
                .annotations
                .get_annotation(transcript.get_seqname(), transcript.get_interval())
                .iter()
                .filter(|a| {
                    a.get_attribute().get_transcript_id()
                        == transcript.get_attribute().get_transcript_id()
                })
                .cloned()
                .collect::<Vec<&annotation::Annotation>>();

            for annotator in &self.transcript_annotators[..] {
                myth.extend_annotation(&annotator.annotate(&transcript_annotations, &variant))
            }
        }

        myth
    }
}
