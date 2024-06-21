//! A Myth store annotation about a variant.

/* std use */

/* crate use */

/* project use */
use crate::variant;

#[derive(Debug, Clone, serde::Serialize)]
pub enum Effect {
    /// A variant in 5′UTR region produces a three base sequence that can be a START codon.
    StartGaine,
    /// The variant hits a splice acceptor site (defined as two bases before exon start, except for the first exon).
    SpliceSiteAcceptor,
    /// The variant hits a Splice donor site (defined as two bases after coding exon end, except for the last exon).
    SplitceSiteDonor,
    /// Variant causes start codon to be mutated into a non-start codon.
    StartLost,
    /// Variant causes start codon to be mutated into another start codon.
    SYNONYMOUS_START,
    /// The variant hits a CDS.
    CDS,
    /// The variant hits a gene.
    GENE,
    /// The variant hits a transcript.
    TRANSCRIPT,
    /// The vairant hist an exon.
    EXON,
    /// A deletion removes the whole exon.
    EXON_DELETED,
    ///Variant causes a codon that produces a different amino acid
    NON_SYNONYMOUS_CODING,
    ///Variant causes a codon that produces the same amino acid
    SYNONYMOUS_CODING,
    ///Insertion or deletion causes a frame shift
    FRAME_SHIFT,
    ///One or many codons are changed
    CODON_CHANGE,
    /// One or many codons are inserted
    CODON_INSERTION,
    /// One codon is changed and one or many codons are inserted
    CODON_CHANGE_PLUS_CODON_INSERTION,
    /// One or many codons are deleted
    CODON_DELETION,
    ///One codon is changed and one or more codons are deleted
    CODON_CHANGE_PLUS_CODON_DELETION,
    /// Variant causes a STOP codon
    STOP_GAINED,
    /// Variant causes stop codon to be mutated into another stop codon.
    SYNONYMOUS_STOP,
    /// Variant causes stop codon to be mutated into a non-stop codon
    STOP_LOST,
    /// Variant hist and intron. Technically, hits no exon in the transcript.
    INTRON,
    /// Variant hits 3′UTR region
    UTR_3_PRIME,
    /// The variant deletes an exon which is in the 3′UTR of the transcript
    UTR_3_DELETED,
    /// Downstream of a gene (default length: 5K bases)
    DOWNSTREAM,
    /// The variant is in a highly conserved intronic region
    INTRON_CONSERVED,
    /// The variant is in a highly conserved intergenic region
    INTERGENIC_CONSERVED,
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

    pub effects: Vec<Effect>,
}

impl AnnotationMyth {
    /// Get builder of AnnotationMyth
    pub fn builder() -> AnnotationMythBuilder {
        AnnotationMythBuilder::default()
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
