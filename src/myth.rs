//! A Myth store annotation about a variant.

/* std use */

/* crate use */

/* project use */
use crate::variant;

/// Effect of variant
#[derive(Debug, Clone, serde::Serialize)]
pub enum Effect {
    /// The variant hits a CDS.
    Cds,
    /// One or many codons are changed
    CodonChange,
    /// One codon is changed and one or more codons are deleted
    CodonChangePlusCodonDeletion,
    /// One codon is changed and one or many codons are inserted
    CodonChangePlusCodonInsertion,
    /// One or many codons are deleted
    CodonDeletion,
    /// One or many codons are inserted
    CodonInsertion,
    /// Downstream of a gene (default length: 5K bases)
    Downstream,
    /// The vairant hist an exon.
    Exon,
    /// A deletion removes the whole exon.
    ExonDeleted,
    ///Insertion or deletion causes a frame shift
    FrameShift,
    /// The variant hits a gene.
    Gene,
    /// Variant hist and intron. Technically, hits no exon in the transcript.
    Intron,
    /// The variant is in an intergenic region
    Intergenic,
    /// The variant is in a highly conserved intergenic region
    IntergenicConserved,
    /// The variant is in a highly conserved intronic region
    IntronConserved,
    /// Variant causes a codon that produces a different amino acid
    NonSynonymousCoding,
    /// The variant hits a splice acceptor site (defined as two bases before exon start, except for the first exon).
    SpliceSiteAcceptor,
    /// The variant hits a Splice donor site (defined as two bases after coding exon end, except for the last exon).
    SplitceSiteDonor,
    /// A variant in 5′UTR region produces a three base sequence that can be a START codon.
    StartGained,
    /// Variant causes start codon to be mutated into a non-start codon.
    StartLost,
    /// Variant causes a STOP codon
    StopGained,
    /// Variant causes stop codon to be mutated into a non-stop codon
    StopLost,
    /// Variant causes a codon that produces the same amino acid
    SynonymousCoding,
    /// Variant causes start codon to be mutated into another start codon.
    SynonymousStart,
    /// Variant causes stop codon to be mutated into another stop codon.
    SynonymousStop,
    /// The variant hits a transcript.
    Transcript,
    /// Upstream of a gene (default length: 5K bases)
    Upstream,
    /// The variant deletes an exon which is in the 3′UTR of the transcript
    Utr3Deleted,
    /// Variant hits 3′UTR region
    Utr3Prime,
    /// The variant deletes an exon which is in the 5′UTR of the transcript
    Utr5Deleted,
    /// Variant hits 5′UTR region
    Utr5Prime,
}

/// Struct to store annotation information
#[derive(Debug, serde::Serialize, derive_builder::Builder, Clone)]
#[builder(pattern = "owned")]
pub struct AnnotationMyth {
    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Source of annotation
    pub source: Vec<u8>,

    #[serde(serialize_with = "crate::serialize_bstr")]
    /// Transcript id
    pub transcript_id: Vec<u8>,

    /// Store effect of this variants
    pub effects: Vec<Effect>,
}

impl AnnotationMyth {
    /// Get builder of AnnotationMyth
    pub fn builder() -> AnnotationMythBuilder {
        AnnotationMythBuilder::default()
    }
}

impl AnnotationMythBuilder {
    /// Add Effect in AnnotationMyth
    pub fn add_effect(&mut self, e: Effect) {
        if let Some(effects) = &mut self.effects {
            effects.push(e)
        } else {
            self.effects = Some(vec![e])
        }
    }

    /// Extend Effect in AnnotationMyth
    pub fn extend_effect(&mut self, e: &[Effect]) {
        if let Some(effects) = &mut self.effects {
            effects.extend_from_slice(e)
        } else {
            self.effects = Some(e.to_vec())
        }
    }
}

/// Store information around variant
#[derive(Debug, serde::Serialize)]
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
