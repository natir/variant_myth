//! Find out which genes cover a variant.

/* std use */

/* crate use */
use bstr::ByteSlice as _;

/* project use */
use crate::annotations_db;
use crate::variant;

/// Struct to store annotation information
#[derive(Debug, serde::Serialize, Clone, PartialEq)]
pub struct GeneMyth {
    /// Variant
    variant: variant::Variant,

    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Gene name
    gene_name: Vec<u8>,

    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Gene id
    gene_id: Vec<u8>,
}

impl GeneMyth {
    /// Create a new GeneMyth
    pub fn new(variant: variant::Variant) -> Self {
        Self {
            variant,
            gene_name: vec![],
            gene_id: vec![],
        }
    }

    /// Add a gene_name in myth
    pub fn add_gene_name(&mut self, gene_name: &[u8]) {
        if !self.gene_name.is_empty() {
            self.gene_name.push(b',');
        }
        self.gene_name.extend(gene_name);
    }

    /// Add a gene_id in myth
    pub fn add_gene_id(&mut self, gene_id: &[u8]) {
        if !self.gene_id.is_empty() {
            self.gene_id.push(b',');
        }
        self.gene_id.extend(gene_id);
    }
}

/// Struct that associate variant 2 gene
pub struct Variant2Gene<'a> {
    annotations: &'a annotations_db::AnnotationsDataBase,
}

impl<'a> Variant2Gene<'a> {
    /// Create a new Variant2Gene
    pub fn new(annotations: &'a annotations_db::AnnotationsDataBase) -> Self {
        Self { annotations }
    }

    /// Methode that generate GeneMyth from variant
    pub fn gene(&self, variant: variant::Variant) -> GeneMyth {
        let mut myth = GeneMyth::new(variant.clone());

        let annotations = self
            .annotations
            .get_annotation(&variant.seqname, variant.get_interval());

        if annotations.is_empty() {
            return myth;
        }

        for annotation in &annotations {
            if annotation.get_feature().contains_str("gene") {
                myth.add_gene_id(annotation.get_attribute().get_id());
                myth.add_gene_name(annotation.get_attribute().get_name());
            }
        }

        myth
    }
}
