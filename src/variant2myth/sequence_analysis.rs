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
    ) -> Vec<effect::Effect> {
        let transcript = annotations
            .iter()
            .filter(|a| a.get_feature() == b"transcript")
            .next()
            .unwrap();

        vec![]
    }
}
