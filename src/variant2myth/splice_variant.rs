use crate::annotation;
use crate::effect;
use crate::utils::CDNAPosition;
use crate::variant;
use crate::variant2myth;

pub struct SpliceVariant {}

fn cdna_annotate(
    position: u64,
    annotations: &[&annotation::Annotation],
    found_effects: &mut Vec<effect::Effect>,
) {
    if let Some(cdna_pos) = CDNAPosition::from_genomic_pos(position, annotations) {
        match cdna_pos {
            CDNAPosition::FivePrimeIntronic {
                last_exon_position: _,
                distance_to_prev_exon,
            } => {
                if distance_to_prev_exon.abs() <= 2 {
                    found_effects.push(effect::Effect::SpliceDonorVariant);
                }
                found_effects.push(effect::Effect::IntronVariant);
            }
            CDNAPosition::ThreePrimeIntronic {
                next_exon_position: _,
                distance_to_next_exon,
            } => {
                if distance_to_next_exon.abs() <= 2 {
                    found_effects.push(effect::Effect::SpliceAcceptorVariant);
                }
                found_effects.push(effect::Effect::IntronVariant);
            }
            //Not intronic
            _ => {}
        }
    } else {
        log::debug!("No cDNA location found");
    }
}

impl variant2myth::Annotator for SpliceVariant {
    fn annotate(
        &self,
        annotations: &[&annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<effect::Effect> {
        let mut res: Vec<_> = vec![];
        cdna_annotate(variant.position, annotations, &mut res);
        res
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */
}
