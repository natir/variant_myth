//! An annotator for upstream variant

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::effect;
use crate::variant;
use crate::variant2myth;

pub struct FeaturePresence<const E: usize> {}

impl<const E: usize> FeaturePresence<E> {
    pub const fn new() -> Self {
        Self {}
    }
}

const fn effect2feature_name<const E: usize>() -> &'static [u8] {
    match effect::usize2effect(E) {
        effect::Effect::UpstreamGeneVariant => b"upstream",
        effect::Effect::DownstreamGeneVariant => b"downstream",
        effect::Effect::FivePrimeUtrVariant => b"five_prime_UTR",
        effect::Effect::ThreePrimeUtrVariant => b"three_prime_UTR",
        _ => panic!("Effect isn't support by feature presence."),
    }
}

impl<const E: usize> variant2myth::Annotator for FeaturePresence<E> {
    fn annotate(
        &self,
        annotations: &[&annotation::Annotation],
        _variant: &variant::Variant,
    ) -> Vec<effect::Effect> {
        if annotations
            .iter()
            .any(|a| a.get_feature() == effect2feature_name::<E>())
        {
            vec![effect::usize2effect(E)]
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
    use crate::variant2myth::feature_presence::effect2feature_name;
    use crate::variant2myth::Annotator as _;

    use super::FeaturePresence;

    #[test]
    fn convert_effect2feature_name() {
        assert_eq!(
            effect2feature_name::<{ effect::Effect::UpstreamGeneVariant as usize }>(),
            b"upstream"
        );
        assert_eq!(
            effect2feature_name::<{ effect::Effect::DownstreamGeneVariant as usize }>(),
            b"downstream"
        );
        assert_eq!(
            effect2feature_name::<{ effect::Effect::FivePrimeUtrVariant as usize }>(),
            b"five_prime_UTR"
        );
        assert_eq!(
            effect2feature_name::<{ effect::Effect::ThreePrimeUtrVariant as usize }>(),
            b"three_prime_UTR"
        );
    }

    #[test]
    fn feature_presence() -> error::Result<()> {
        let file = GFF.replace(b"{0}", b"chr1");

        let obj = FeaturePresence::<{ effect::Effect::FivePrimeUtrVariant as usize }>::new();
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

        let obj = FeaturePresence::<{ effect::Effect::ThreePrimeUtrVariant as usize }>::new();
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
