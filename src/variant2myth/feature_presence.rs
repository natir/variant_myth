//! An annotator for upstream variant

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::effect;
use crate::myth;
use crate::variant;
use crate::variant2myth;

pub struct FeaturePresence {
    name: &'static [u8],
    effect: effect::Effect,
}

impl FeaturePresence {
    pub const fn new(name: &'static [u8], effect: effect::Effect) -> Self {
        Self { name, effect }
    }
}

impl variant2myth::Annotator for FeaturePresence {
    fn annotate(
        &self,
        annotations: &[&annotation::Annotation],
        _variant: &variant::Variant,
    ) -> Vec<myth::AnnotationMyth> {
        let mut result = Vec::new();
        for up_annot in annotations.iter().filter(|a| a.get_feature() == self.name) {
            result.push(
                myth::AnnotationMyth::builder()
                    .source(up_annot.get_source().to_vec())
                    .transcript_id(up_annot.get_transcript_id().to_vec())
                    .effects(vec![self.effect.clone()])
                    .build()
                    .unwrap(), // No possible error in build
            )
        }

        result
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */
    use super::*;
}
