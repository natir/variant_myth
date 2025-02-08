//! Ann annotator to detect start stop lost

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::effect;
use crate::sequences_db;
use crate::translate;
use crate::variant;
use crate::variant2myth;

pub struct StartStopLost<'a, const E: usize> {
    translate: &'a translate::Translate,
    sequences: &'a sequences_db::SequencesDataBase,
}

impl<'a, const E: usize> StartStopLost<'a, E> {
    pub fn new(
        translate: &'a translate::Translate,
        sequences: &'a sequences_db::SequencesDataBase,
    ) -> Self {
        const {
            assert!(
                E == effect::Effect::StartLost as usize || E == effect::Effect::StopLost as usize,
                "StartStopLost struct can only manage Effect::StartLost and StopLost"
            )
        }

        Self {
            translate,
            sequences,
        }
    }
}

impl<const E: usize> variant2myth::Annotator for StartStopLost<'_, E> {
    fn annotate(
        &self,
        annotations: &[&annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<effect::Effect> {
        const STARTLOST: usize = effect::Effect::StartLost as usize;
        const STOPLOST: usize = effect::Effect::StopLost as usize;
        let feature_name: &[u8] = match E {
            STARTLOST => b"start_position",
            STOPLOST => b"stop_position",
            _ => panic!("StartStopLost struct can only manage Effect::StartLost and StopLost"),
        };

        // Get start range
        let (feature_range, seq_name, strand) = if let Some(start_annotation) = annotations
            .iter()
            .find(|annotation| annotation.get_feature() == feature_name)
        {
            (
                start_annotation.get_interval(),
                start_annotation.get_seqname(),
                start_annotation.get_strand(),
            )
        } else if let Some(exon) = annotations.iter().find(|a| a.get_feature() == b"exon") {
            (
                exon.get_start()..exon.get_start() + 3,
                exon.get_seqname(),
                exon.get_strand(),
            )
        } else {
            return vec![]; // no start stop position, no exon -> no annotation
        };

        let variant_range = variant.get_interval();
        if variant_range.start > feature_range.end && variant_range.end < feature_range.start {
            // Variant overlap feature
            let feature_seq = if let Ok(seq) = self.sequences.get_interval(seq_name, &feature_range)
            {
                if strand == &annotation::Strand::Forward {
                    seq.to_vec()
                } else {
                    let mut val = seq.to_vec();
                    sequences_db::rev_comp(&mut val);
                    val
                }
            } else {
                return vec![];
            };

            vec![]
        } else {
            // Variant not overlap feature
            vec![]
        }
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */
    use super::*;
}
