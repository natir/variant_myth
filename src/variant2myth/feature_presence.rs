//! An annotator for upstream variant

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::effect;
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
    use bstr::ByteSlice as _;

    /* project use */
    use crate::annotation;
    use crate::effect;
    use crate::error;
    use crate::tests::GFF;
    use crate::variant;
    use crate::variant2myth::Annotator as _;

    use super::FeaturePresence;

    #[test]
    fn feature_presence() -> error::Result<()> {
        let file = GFF.replace(b"{0}", b"chr1");

        let obj = FeaturePresence::new(b"five_prime_UTR", effect::Effect::FivePrimeUtrVariant);
        let annotation = annotation::Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(file[292..374].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        assert_eq!(
            obj.annotate(
                &[&annotation],
                &variant::Variant::test_variant(b"chr1", 900, b"A", b"C")
            ),
            vec![effect::Effect::FivePrimeUtrVariant]
        );

        Ok(())
    }

    #[test]
    fn feature_absence() -> error::Result<()> {
        let file = GFF.replace(b"{0}", b"chr1");

        let obj = FeaturePresence::new(b"three_prime_UTR", effect::Effect::ThreePrimeUtrVariant);
        let annotation = annotation::Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(file[292..374].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        assert_eq!(
            obj.annotate(
                &[&annotation],
                &variant::Variant::test_variant(b"chr1", 900, b"A", b"C")
            ),
            vec![]
        );

        Ok(())
    }
}
