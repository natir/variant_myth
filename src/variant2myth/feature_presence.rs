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
    ) -> Vec<effect::Effect> {
        if annotations.iter().any(|a| a.get_feature() == self.name) {
            vec![self.effect.clone()]
        } else {
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
