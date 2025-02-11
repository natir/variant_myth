//! An annotator for upstream variant

/* std use */

/* crate use */

/* project use */
use crate::effect;
use crate::memoizor;
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
        _variant: &variant::Variant,
        memoizor: &mut memoizor::Memoizor,
    ) -> Vec<effect::Effect> {
        if memoizor
            .not_coding_annotation()
            .iter()
            .any(|a| a.get_feature() == self.name)
        {
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
    use crate::annotations_db;
    use crate::effect;
    use crate::error;
    use crate::memoizor;
    use crate::sequences_db;
    use crate::test_data::GFF;
    use crate::test_data::SEQUENCE;
    use crate::variant;
    use crate::variant2myth::Annotator as _;

    use super::FeaturePresence;

    #[test]
    fn feature_presence() -> error::Result<()> {
        let obj = FeaturePresence::new(b"five_prime_UTR", effect::Effect::FivePrimeUtrVariant);

        let reader: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>> =
            std::io::BufReader::new(Box::new(GFF));
        let annotations_db = annotations_db::AnnotationsDataBase::from_reader(reader, 100)?;
        let reader: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>> =
            std::io::BufReader::new(Box::new(SEQUENCE));
        let sequences_db = sequences_db::SequencesDataBase::from_reader(reader)?;

        let variant = variant::Variant::test_variant(b"chrA", 56, b"A", b"C", None)?;
        let not_coding_annotation =
            annotations_db.get_annotations(&variant.seqname, variant.get_interval());

        let mut memoizor = memoizor::Memoizor::new(
            b"ENST00000797271.1",
            &annotations_db,
            &sequences_db,
            &not_coding_annotation,
        );

        assert_eq!(
            obj.annotate(&variant, &mut memoizor),
            vec![effect::Effect::FivePrimeUtrVariant]
        );

        Ok(())
    }

    #[test]
    fn feature_absence() -> error::Result<()> {
        let obj = FeaturePresence::new(b"three_prime_UTR", effect::Effect::ThreePrimeUtrVariant);

        let reader: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>> =
            std::io::BufReader::new(Box::new(GFF));
        let annotations_db = annotations_db::AnnotationsDataBase::from_reader(reader, 100)?;
        let reader: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>> =
            std::io::BufReader::new(Box::new(SEQUENCE));
        let sequences_db = sequences_db::SequencesDataBase::from_reader(reader)?;

        let variant = variant::Variant::test_variant(b"chrA", 56, b"A", b"C", None)?;
        let not_coding_annotation =
            annotations_db.get_annotations(&variant.seqname, variant.get_interval());

        let mut memoizor = memoizor::Memoizor::new(
            b"ENST00000797271.1",
            &annotations_db,
            &sequences_db,
            &not_coding_annotation,
        );

        assert_eq!(obj.annotate(&variant, &mut memoizor), vec![]);

        Ok(())
    }
}
