//! An annotator for upstream variant

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::effect;
use crate::myth;
use crate::variant;
use crate::variant2myth;

pub struct SequenceAnalysis {}

impl SequenceAnalysis {
    pub const fn new() -> Self {
        Self {}
    }
}

impl variant2myth::Annotator for SequenceAnalysis {
    fn annotate(
        &self,
        annotations: &[&annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<myth::AnnotationMyth> {
        let transcript = annotations
            .iter()
            .filter(|a| a.get_feature() == b"transcript")
            .next()
            .unwrap();

        let mut transcript_myth = myth::AnnotationMyth::builder()
            .source(transcript.get_source().to_vec())
            .transcript_id(transcript.get_transcript_id().to_vec())
            .gene_name(transcript.get_attribute().get_gene_name().to_vec())
            .effects(vec![]);

        vec![transcript_myth.build().unwrap()]
    }
}
