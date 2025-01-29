use crate::annotation;
use crate::annotation::Annotation;
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

fn genomic_to_cdna(position: u64, exons: &[&Annotation]) -> (u64, u64) {
    enum Location {
        Intronic,
        Exonic,
    }
    let (mut last_exon_pos, mut pos) = (0, 0);
    let mut exons = exons.to_vec();
    exons.sort_by(|a, b| {
        a.get_start()
            .partial_cmp(&b.get_start())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut location = Location::Exonic;
    pos = match exons[0].get_strand() {
        Strand::Forward => exons[0].get_start() - position,
        Strand::Reverse => exons[0].get_start() - position,
    };
    for exon in exons {
        if (exon.get_start()..exon.get_stop()).contains(&position) {
            last_exon_pos = pos;
            pos = match exon.get_strand() {
                Strand::Forward => exons[0].get_start() - position,
                Strand::Reverse => exons[0].get_start() - position,
            };
            location = Location::Exonic;
        } else {
            location = Location::Intronic;
        }
    }
    match exons[0].get_strand() {
        Strand::Forward => {}
        Strand::Reverse => {}
    }
    (last_exon_pos, pos)
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
