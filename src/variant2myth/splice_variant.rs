use crate::annotation;
use crate::annotation::Strand;
use crate::effect;
use crate::variant;
use crate::variant2myth;
use std::ops::Not;

pub struct SpliceVariant {}

impl SpliceVariant {
    pub const fn new() -> Self {
        Self {}
    }
}

impl variant2myth::Annotator for SpliceVariant {
    fn annotate(
        &self,
        annotations: &[&annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<effect::Effect> {
        let is_intronic = annotations
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .map(|a| {
                (a.get_start()..a.get_stop())
                    .contains(&(variant.position + (variant.ref_seq.len() as u64) - 1))
            })
            .any(|x| x)
            .not();
        let (exonstarts, exon_ends): (Vec<_>, Vec<_>) = annotations
            .iter()
            .filter(|a| a.get_feature() == b"exon")
            .map(|a| match a.get_strand() {
                Strand::Forward => (a.get_start(), a.get_stop()),
                Strand::Reverse => (a.get_stop(), a.get_start()),
            })
            .unzip();
        for ann in annotations {
            println!("{}", ann);
        }
        vec![]
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */
}
