//! An annotator for upstream variant

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::effect;
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
        _variant: &variant::Variant,
    ) -> Vec<effect::Effect> {
        let _start_position = annotations
            .iter()
            .find(|&&x| x.get_feature() == b"start_codon")
            .map(|x| x.get_start());

        let _stop_position = annotations
            .iter()
            .find(|&&x| x.get_feature() == b"stop_codon")
            .map(|x| x.get_stop());

        let _exon_annot = annotations
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        vec![]
    }
}
