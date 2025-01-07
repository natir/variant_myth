//! An annotator for upstream variant

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::effect;
use crate::sequences_db;
use crate::translate;
use crate::variant;
use crate::variant2myth;

pub struct SequenceAnalysis<'a> {
    translate: &'a translate::Translate,
    sequences: &'a sequences_db::SequencesDataBase,
}

impl<'a> SequenceAnalysis<'a> {
    pub const fn new(
        translate: &'a translate::Translate,
        sequences: &'a sequences_db::SequencesDataBase,
    ) -> Self {
        Self {
            translate,
            sequences,
        }
    }
}

impl variant2myth::Annotator for SequenceAnalysis<'_> {
    fn annotate(
        &self,
        annotations: &[&annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<effect::Effect> {
        let start_position = annotations
            .iter()
            .find(|&&x| x.get_feature() == b"start_codon")
            .map(|x| x.get_start());

        let stop_position = annotations
            .iter()
            .find(|&&x| x.get_feature() == b"stop_codon")
            .map(|x| x.get_stop());

        let _exon_annot = annotations
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .cloned()
            .collect::<Vec<&annotation::Annotation>>();

        let strand = annotations
            .first()
            .map(|x| x.get_strand())
            .unwrap_or(&annotation::Strand::Forward);

        let _epissed = match self.sequences.epissed(annotations, *strand) {
            Ok(sequence) => sequence,
            Err(error) => {
                log::error!("{:?}", error);
                return vec![];
            }
        };

        let _epissed_var = match self.sequences.epissed_edit(annotations, *strand, variant) {
            Ok(sequence) => sequence,
            Err(error) => {
                log::error!("{:?}", error);
                return vec![];
            }
        };

        let coding =
            match self
                .sequences
                .coding(annotations, *strand, start_position, stop_position)
            {
                Ok(sequence) => sequence,
                Err(error) => {
                    log::error!("{:?}", error);
                    return vec![];
                }
            };

        let coding_var = match self.sequences.coding_edit(
            annotations,
            *strand,
            variant,
            start_position,
            stop_position,
        ) {
            Ok(sequence) => sequence,
            Err(error) => {
                log::error!("{:?}", error);
                return vec![];
            }
        };

        let translate = self.translate.translate(&coding);
        let _translate_var = self.translate.translate(&coding_var);

        log::trace!("{}", String::from_utf8(translate).unwrap());

        vec![]
    }
}
