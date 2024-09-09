//! A Myth store annotation about a variant.

/* std use */

/* crate use */

/* project use */
use crate::effect;
use crate::variant;

/// Struct to store annotation information
#[derive(Debug, serde::Serialize, derive_builder::Builder, Clone, PartialEq)]
#[builder(pattern = "owned")]
pub struct AnnotationMyth {
    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Source of annotation
    pub source: Vec<u8>,

    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Transcript id
    pub transcript_id: Vec<u8>,

    #[serde(serialize_with = "crate::serialize_bstr")]
    #[builder(default)]
    /// Gene name
    pub gene_name: Vec<u8>,

    /// Store effect of this variants
    effects: Vec<effect::Effect>,

    #[builder(private, default)]
    /// Store impact of effect
    impact: effect::Impact,
}

impl AnnotationMyth {
    /// Get builder of AnnotationMyth
    pub fn builder() -> AnnotationMythBuilder {
        AnnotationMythBuilder::default()
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
#[derive(Debug, serde::Serialize, PartialEq)]
pub struct Myth {
    variant: variant::Variant,
    annotations: Vec<AnnotationMyth>,
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

    #[test]
    fn annotation_myth() {
        let annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .transcript_id(b"gene1".to_vec())
            .effects(vec![])
            .build()
            .unwrap();

        assert_eq!(annotation.source, b"test".to_vec());
        assert_eq!(annotation.transcript_id, b"gene1".to_vec());
        assert_eq!(annotation.gene_name, Vec::<u8>::new());
        assert_eq!(annotation.effects, Vec::<effect::Effect>::new());

        let mut annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .transcript_id(b"gene1".to_vec());

        annotation.add_effect(effect::Effect::GeneVariant);
        annotation.add_effect(effect::Effect::ExonRegion);

        assert_eq!(
            annotation.build().unwrap(),
            AnnotationMyth {
                source: b"test".to_vec(),
                transcript_id: b"gene1".to_vec(),
                gene_name: b"".to_vec(),
                effects: vec![effect::Effect::GeneVariant, effect::Effect::ExonRegion],
                impact: effect::Impact::Modifier,
            }
        );

        let mut annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .transcript_id(b"gene1".to_vec());

        annotation.extend_effect(&[effect::Effect::GeneVariant, effect::Effect::ExonRegion]);

        assert_eq!(
            annotation.build().unwrap(),
            AnnotationMyth {
                source: b"test".to_vec(),
                transcript_id: b"gene1".to_vec(),
                gene_name: b"".to_vec(),
                effects: vec![effect::Effect::GeneVariant, effect::Effect::ExonRegion],
                impact: effect::Impact::Modifier,
            }
        )
    }

    #[test]
    fn myth() {
        let mut annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .transcript_id(b"gene1".to_vec());

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
                    transcript_id: b"gene1".to_vec(),
                    gene_name: b"".to_vec(),
                    effects: vec![effect::Effect::GeneVariant, effect::Effect::ExonRegion],
                    impact: effect::Impact::Modifier
                }]
            }
        );
    }

    #[test]
    fn json_serde() -> crate::error::Result<()> {
        let mut annotation = AnnotationMyth::builder()
            .source(b"test".to_vec())
            .transcript_id(b"gene1".to_vec());

        annotation.add_effect(effect::Effect::GeneVariant);
        annotation.add_effect(effect::Effect::ExonRegion);

        let mut myth = Myth::from_variant(variant::Variant {
            seqname: b"93".to_vec(),
            position: 2036067340,
            ref_seq: b"T".to_vec(),
            alt_seq: b".".to_vec(),
        });

        myth.add_annotation(annotation.build().unwrap());

        let mut json = vec![];
        serde_json::to_writer(&mut json, &myth)?;
        assert_eq!(json, b"{\"variant\":{\"seqname\":\"93\",\"position\":2036067340,\"ref_seq\":\"T\",\"alt_seq\":\".\"},\"annotations\":[{\"source\":\"test\",\"transcript_id\":\"gene1\",\"gene_name\":\"\",\"effects\":[\"GeneVariant\",\"ExonRegion\"],\"impact\":\"Modifier\"}]}");

        Ok(())
    }
}
