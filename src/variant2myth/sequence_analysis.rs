//! An annotator for upstream variant

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::effect;
use crate::memoizor;
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
        variant: &variant::Variant,
        memoizor: &mut memoizor::Memoizor,
    ) -> Vec<effect::Effect> {
        let start_position = memoizor
            .coding_annotation()
            .unwrap()
            .iter()
            .find(|&x| x.get_feature() == b"start_codon")
            .map(|x| x.get_start());

        let stop_position = memoizor
            .coding_annotation()
            .unwrap()
            .iter()
            .find(|&x| x.get_feature() == b"stop_codon")
            .map(|x| x.get_stop());

        let strand = memoizor
            .transcript()
            .map(|x| x.get_strand())
            .cloned()
            .unwrap_or(annotation::Strand::Forward);

        let exon_annot = memoizor.exons_annotation();
        let exon_annot_proxy = exon_annot.iter().collect::<Vec<&annotation::Annotation>>();

        // No exon in associate annotation no sequence analysis
        if exon_annot_proxy.is_empty() {
            return vec![];
        }

        let coding =
            match self
                .sequences
                .coding(&exon_annot_proxy, strand, start_position, stop_position)
            {
                Ok(sequence) => sequence,
                Err(error) => {
                    log::error!("{:?}", error);
                    return vec![];
                }
            };

        let coding_var = match self.sequences.coding_edit(
            &exon_annot_proxy,
            strand,
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

        // Coding sequence not change variant haven't impact
        if coding == coding_var {
            return vec![];
        }

        let translate = self.translate.translate(&coding);
        let translate_var = self.translate.translate(&coding_var);

        if translate != translate_var {
            log::debug!("ORIGINAL: {}", String::from_utf8(translate).unwrap());
            log::debug!("EDIT    : {}", String::from_utf8(translate_var).unwrap());
            log::debug!("VARIANT : {}", variant);
            log::debug!(
                "ID      : {}",
                String::from_utf8(
                    memoizor.coding_annotation().unwrap()[0]
                        .get_attribute()
                        .get_id()
                        .to_vec()
                )
                .unwrap()
            );
        }

        vec![]
    }
}
