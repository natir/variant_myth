//! A Myth store annotation about a variant.

/* std use */

/* crate use */

use serde::Serialize;

/* project use */
use crate::variant;

/// Struct to store annotation information
#[derive(Debug, Serialize, derive_builder::Builder, Clone)]
#[builder(pattern = "owned")]
pub struct AnnotationMyth {
    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Source of annotation
    pub source: Vec<u8>,

    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Transcript id
    pub transcript_id: Vec<u8>,
}

impl AnnotationMyth {
    /// Get builder of AnnotationMyth
    pub fn builder() -> AnnotationMythBuilder {
        AnnotationMythBuilder::default()
    }
}

/// Store information around variant
#[derive(Debug, Serialize)]
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

    /// Add transcript from source
    pub fn add_annotation(&mut self, source: AnnotationMyth) {
        self.annotations.push(source);
    }
}
