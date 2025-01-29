//! A Myth store annotation about a variant.

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::effect;
use crate::variant;

/// Struct to store annotation information
#[derive(Debug, derive_builder::Builder, Clone, PartialEq)]
#[cfg_attr(feature = "out_json", derive(serde::Serialize))]
#[builder(pattern = "owned")]
pub struct AnnotationMyth {
    /// Source of annotation
    #[cfg_attr(feature = "out_json", serde(serialize_with = "crate::serialize_bstr"))]
    pub source: Vec<u8>,

    /// Feature type
    #[cfg_attr(feature = "out_json", serde(serialize_with = "crate::serialize_bstr"))]
    pub feature: Vec<u8>,

    /// Feature id
    #[cfg_attr(feature = "out_json", serde(serialize_with = "crate::serialize_bstr"))]
    pub id: Vec<u8>,

    #[builder(default)]
    #[cfg_attr(feature = "out_json", serde(serialize_with = "crate::serialize_bstr"))]
    /// Feature name
    pub name: Vec<u8>,

    /// Store effect of this variants
    pub effects: Vec<effect::Effect>,

    #[builder(private, default)]
    /// Store impact of effect
    pub impact: effect::Impact,
}

impl AnnotationMyth {
    /// Get builder of AnnotationMyth
    pub fn builder() -> AnnotationMythBuilder {
        AnnotationMythBuilder::default()
    }

    /// Create a preset builder from annotation of AnnotationMyth, effect must be set before build
    pub fn from_annotation(annotation: &annotation::Annotation) -> AnnotationMythBuilder {
        AnnotationMythBuilder::default()
            .source(annotation.get_source().to_vec())
            .feature(annotation.get_feature().to_vec())
            .id(annotation.get_attribute().get_id().to_vec())
            .name(annotation.get_attribute().get_name().to_vec())
    }

    /// Create a preset builder of AnnotationMyth without annotation
    pub fn from_nowhere() -> AnnotationMythBuilder {
        AnnotationMythBuilder::default()
            .source(b"variant_myth".to_vec())
            .feature(b"unknow".to_vec())
            .id(b"".to_vec())
            .name(b"".to_vec())
    }
}

impl AnnotationMythBuilder {
    /// Add Effect in AnnotationMyth
    pub fn add_effect(&mut self, e: effect::Effect) {
        if let Some(effects) = &mut self.effects {
            self.impact = Some(core::cmp::max(
                effect::Impact::from(&e),
                self.impact.clone().unwrap_or(effect::Impact::Other),
            ));
            effects.push(e);
        } else {
            self.impact = Some(effect::Impact::from(&e));
            self.effects = Some(vec![e]);
        }
    }

    /// Extend Effect in AnnotationMyth
    pub fn extend_effect(&mut self, e: &[effect::Effect]) {
        if let Some(effects) = &mut self.effects {
            effects.extend_from_slice(e);
            self.impact = effects.iter().map(effect::Impact::from).max();
        } else {
            self.effects = Some(e.to_vec());
            self.impact = e.iter().map(effect::Impact::from).max();
        }
    }
}

/// Store information around variant
#[derive(Debug, PartialEq)]
#[cfg_attr(feature = "out_json", derive(serde::Serialize))]
pub struct Myth {
    /// Variant associate to Myth
    pub variant: variant::Variant,
    /// Annotation Myth associate to variant
    pub annotations: Vec<AnnotationMyth>,
}

impl Myth {
    /// Build a Myth from variant
    pub fn from_variant(variant: variant::Variant) -> Self {
        Myth {
            variant,
            annotations: vec![],
        }
    }

    /// Add annotation to variant
    pub fn add_annotation(&mut self, source: AnnotationMyth) {
        self.annotations.push(source);
    }

    /// Extend annotation to variant
    pub fn extend_annotation(&mut self, source: &[AnnotationMyth]) {
        self.annotations.extend_from_slice(source)
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */
    use super::*;
    use crate::error;

    #[test]
    fn annotation_myth() {
        let annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .feature(b"gene".to_vec())
            .name(b"gene1".to_vec())
            .id(b"1111".to_vec())
            .effects(vec![])
            .build()
            .unwrap();

        assert_eq!(annotation.source, b"test".to_vec());
        assert_eq!(annotation.feature, b"gene".to_vec());
        assert_eq!(annotation.name, b"gene1".to_vec());
        assert_eq!(annotation.id, b"1111".to_vec());
        assert_eq!(annotation.effects, Vec::<effect::Effect>::new());

        let mut annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .feature(b"gene".to_vec())
            .name(b"gene1".to_vec())
            .id(b"11111".to_vec());

        annotation.add_effect(effect::Effect::GeneVariant);
        annotation.add_effect(effect::Effect::ExonRegion);

        assert_eq!(
            annotation.build().unwrap(),
            AnnotationMyth {
                source: b"test".to_vec(),
                feature: b"gene".to_vec(),
                name: b"gene1".to_vec(),
                id: b"11111".to_vec(),
                effects: vec![effect::Effect::GeneVariant, effect::Effect::ExonRegion],
                impact: effect::Impact::Modifier,
            }
        );

        let mut annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .feature(b"gene".to_vec())
            .name(b"gene1".to_vec())
            .id(b"1111".to_vec());

        annotation.extend_effect(&[effect::Effect::GeneVariant, effect::Effect::ExonRegion]);

        assert_eq!(
            annotation.build().unwrap(),
            AnnotationMyth {
                source: b"test".to_vec(),
                feature: b"gene".to_vec(),
                name: b"gene1".to_vec(),
                id: b"1111".to_vec(),
                effects: vec![effect::Effect::GeneVariant, effect::Effect::ExonRegion],
                impact: effect::Impact::Modifier,
            }
        )
    }

    #[test]
    fn annotation_myth_from_annotation() -> error::Result<()> {
        let gff_annotation =
            annotation::Annotation::from_byte_record(&csv::ByteRecord::from(vec![
                "chr1",
                "knownGene",
                "transcript",
                "29554",
                "31097",
                ".",
                "+",
                ".",
                "Parent=ENST00000473358.1;ID=11;Name=ENST00001",
            ]))?;

        let mut annotation = AnnotationMyth::from_annotation(&gff_annotation);
        annotation.add_effect(effect::Effect::ExonRegion);

        assert_eq!(
            annotation.build()?,
            AnnotationMyth {
                source: b"knownGene".to_vec(),
                feature: b"transcript".to_vec(),
                name: b"ENST00001".to_vec(),
                id: b"11".to_vec(),
                effects: vec![effect::Effect::ExonRegion],
                impact: effect::Impact::Modifier,
            }
        );

        Ok(())
    }

    #[test]
    fn annotation_myth_from_nowhere() -> error::Result<()> {
        let mut annotation = AnnotationMyth::from_nowhere();
        annotation.add_effect(effect::Effect::Ignore);

        assert_eq!(
            annotation.build()?,
            AnnotationMyth {
                source: b"variant_myth".to_vec(),
                feature: b"unknow".to_vec(),
                name: b"".to_vec(),
                id: b"".to_vec(),
                effects: vec![effect::Effect::Ignore],
                impact: effect::Impact::Other,
            }
        );

        Ok(())
    }

    #[test]
    fn myth() {
        let mut annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .feature(b"gene".to_vec())
            .name(b"gene1".to_vec())
            .id(b"1111".to_vec());

        annotation.add_effect(effect::Effect::GeneVariant);
        annotation.add_effect(effect::Effect::ExonRegion);

        let mut myth = Myth::from_variant(variant::Variant {
            seqname: b"93".to_vec(),
            position: 2036067340,
            ref_seq: b"T".to_vec(),
            alt_seq: b".".to_vec(),
        });

        myth.add_annotation(annotation.build().unwrap());

        assert_eq!(
            myth,
            Myth {
                variant: variant::Variant {
                    seqname: b"93".to_vec(),
                    position: 2036067340,
                    ref_seq: b"T".to_vec(),
                    alt_seq: b".".to_vec(),
                },
                annotations: vec![AnnotationMyth {
                    source: b"test".to_vec(),
                    feature: b"gene".to_vec(),
                    name: b"gene1".to_vec(),
                    id: b"1111".to_vec(),
                    effects: vec![effect::Effect::GeneVariant, effect::Effect::ExonRegion],
                    impact: effect::Impact::Modifier
                }]
            }
        );
    }
}
