//! Enumeration of effect and impact

/* std use */

/* crate use */

/* project use */

/// Impact of variant
#[derive(Default, Debug, Clone, serde::Serialize, PartialEq, Eq, PartialOrd, Ord)]
pub enum Impact {
    /// Variant have an High Impact
    High = 4,
    /// Varinat have an Moderate Impact
    Moderate = 3,
    /// Variant have a Low Impact
    Low = 2,
    /// Variant have an Modifier Impact
    Modifier = 1,
    #[default]
    /// Variant have an Undeterminate Impact
    Other = 0,
}

impl From<&Effect> for Impact {
    fn from(value: &Effect) -> Impact {
        match value {
            // High Impact
		Effect::BidirectionalGeneFusion
		| Effect::ChromosomeNumberVariation
		| Effect::ExonLossVariant
		| Effect::FeatureAblation
		| Effect::FrameshiftVariant
		| Effect::GeneFusion
		| Effect::ProteinProteinContact
		| Effect::RareAminoAcidVariant
		| Effect::RearrangedAtDnaLevel
		| Effect::SpliceAcceptorVariant
		| Effect::SpliceDonorVariant
		| Effect::StartLost
		| Effect::StopGained
		| Effect::StopLost
		| Effect::StructuralInteractionVariant
		| Effect::TranscriptAblation
		=> Impact::High,
            // Moderate Impact
		Effect::ConservativeInframeDeletion
		| Effect::ConservativeInframeInsertion
		| Effect::DisruptiveInframeDeletion
		| Effect::MissenseVariant
		| Effect::P3PrimeUtrTruncation
		| Effect::P5PrimeUtrTruncation
            | Effect::DisruptiveInframeInsertion
		=> Impact::Moderate,
            // Low Impact
		Effect::InitiatorCodonVariant
		| Effect::P5PrimeUtrPrematureStartCodonGainVariant
		| Effect::SpliceRegionVariant
		| Effect::SynonymousVariant
		| Effect::TfBindingSiteVariant
		| Effect::TfbsAblation
		| Effect::FeatureFusion
		=> Impact::Low,
            // Modifier Impact
		Effect::Chromosome
		| Effect::CodingSequenceVariant
		| Effect::ConservedIntergenicVariant
		| Effect::ConservedIntronVariant
		| Effect::DownstreamGeneVariant
		| Effect::ExonRegion
		| Effect::FeatureElongation
		| Effect::GeneVariant
		| Effect::IntergenicRegion
		| Effect::IntragenicVariant
		| Effect::IntronVariant
		| Effect::NonCodingTranscriptExonVariant
		| Effect::NonCodingTranscriptVariant
		| Effect::P3PrimeUtrVariant
		| Effect::P5PrimeUtrVariant
		| Effect::RegulatoryRegionVariant
		| Effect::SequenceFeature
		| Effect::UpstreamGeneVariant
		=> Impact::Modifier,
            // Other Impact
		Effect::Inversion // Large, Exon -> High | Gene, Transcript -> Moderate
		| Effect::MiRna //
		| Effect::StartRetainedVariant // FRAME_SHIFT_BEFORE_CDS_START -> Modifier | SYNONYMOUS_START -> Low
		| Effect::StopRetainedVariant // FrameShiftAfterCDS -> Modifer | NonSynonymousStop -> Low | SynonymousStop -> Low
		| Effect::Duplication // Large, Exon -> High | Gene, Transcript -> Moderate
		=> Impact::Other,

        }
    }
}

/// Effect of variant
#[derive(Debug, Clone, serde::Serialize, PartialEq)]
pub enum Effect {
    /// A sequence variant whereby two genes, on alternate strands have become joined.
    BidirectionalGeneFusion,
    /// Structural unit composed of a nucleic acid molecule which controls its own replication through the interaction of specific proteins at one or more origins of replication.
    Chromosome,
    /// A kind of chromosome variation where the chromosome complement is not an exact multiple of the haploid number.
    ChromosomeNumberVariation,
    /// A sequence variant that changes the coding sequence.
    CodingSequenceVariant,
    /// An inframe decrease in cds length that deletes one or more entire codons from the coding sequence but does not change any remaining codons.
    ConservativeInframeDeletion,
    /// An inframe increase in cds length that inserts one or more codons into the coding sequence between existing codons.
    ConservativeInframeInsertion,
    /// A sequence variant located in a conserved intergenic region, between genes.
    ConservedIntergenicVariant,
    /// A transcript variant occurring within a conserved region of an intron.
    ConservedIntronVariant,
    /// An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon.
    DisruptiveInframeDeletion,
    /// An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon.
    DisruptiveInframeInsertion,
    /// A sequence variant located 3' of a gene.
    DownstreamGeneVariant,
    /// An insertion which derives from, or is identical in sequence to, nucleotides present at a known location in the genome.
    Duplication,
    /// A sequence variant whereby an exon is lost from the transcript.
    ExonLossVariant,
    /// A region of an exon.
    ExonRegion,
    /// A sequence variant, caused by an alteration of the genomic sequence, where the deletion, is greater than the extent of the underlying genomic features.
    FeatureAblation,
    /// A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence.
    FeatureElongation,
    /// A sequence variant, caused by an alteration of the genomic sequence, where a deletion fuses genomic features.
    FeatureFusion,
    /// A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three.
    FrameshiftVariant,
    /// A sequence variant whereby a two genes have become joined.
    GeneFusion,
    /// A sequence variant where the structure of the gene is changed.
    GeneVariant,
    /// A codon variant that changes at least one base of the first codon of a transcript.
    InitiatorCodonVariant,
    /// A region containing or overlapping no genes that is bounded on either side by a gene, or bounded by a gene and the end of the chromosome.
    IntergenicRegion,
    /// A sequence variant located in the intergenic region, between genes.
    IntragenicVariant,
    /// A transcript variant occurring within an intron.
    IntronVariant,
    /// A continuous nucleotide sequence is inverted in the same position.
    Inversion,
    /// Small, ~22-nt, RNA molecule that is the endogenous transcript of a miRNA gene (or the product of other non coding RNA genes). Micro RNAs are produced from precursor molecules (SO:0001244) that can form local hairpin structures, which ordinarily are processed (usually via the Dicer pathway) such that a single miRNA molecule accumulates from one arm of a hairpin precursor molecule. Micro RNAs may trigger the cleavage of their target molecules or act as translational repressors.
    MiRna,
    /// A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved.
    MissenseVariant,
    /// A sequence variant that changes non-coding exon sequence in a non-coding transcript.
    NonCodingTranscriptExonVariant,
    /// A transcript variant of a non coding RNA gene.
    NonCodingTranscriptVariant,
    /// A sequence variant that causes the reduction of a the 3' UTR with regard to the reference sequence.
    P3PrimeUtrTruncation,
    /// A UTR variant of the 3' UTR.
    P3PrimeUtrVariant,
    /// A 5' UTR variant where a premature start codon is introduced, moved or lost.
    P5PrimeUtrPrematureStartCodonGainVariant,
    /// A sequence variant that causes the reduction of a the 5'UTR with regard to the reference sequence.
    P5PrimeUtrTruncation,
    /// A UTR variant of the 5' UTR.
    P5PrimeUtrVariant,
    /// A binding site that, in the protein molecule, interacts selectively and non-covalently with polypeptide residues.
    ProteinProteinContact,
    /// A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid.
    RareAminoAcidVariant,
    /// An attribute to describe the sequence of a feature, where the DNA is rearranged.
    RearrangedAtDnaLevel,
    /// A sequence variant located within a regulatory region.
    RegulatoryRegionVariant,
    /// Any extent of continuous biological sequence.
    SequenceFeature,
    /// A splice variant that changes the 2 base region at the 3' end of an intron.
    SpliceAcceptorVariant,
    /// A splice variant that changes the 2 base pair region at the 5' end of an intron.
    SpliceDonorVariant,
    /// A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron.
    SpliceRegionVariant,
    /// A codon variant that changes at least one base of the canonical start codon.
    StartLost,
    /// A sequence variant where at least one base in the start codon is changed, but the start remains.
    StartRetainedVariant,
    /// A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened polypeptide.
    StopGained,
    /// A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript.
    StopLost,
    /// A sequence variant where at least one base in the terminator codon is changed, but the terminator remains.
    StopRetainedVariant,
    /// A variant that impacts the internal interactions of the resulting polypeptide structure.
    StructuralInteractionVariant,
    /// A sequence variant where there is no resulting change to the encoded amino acid.
    SynonymousVariant,
    /// A sequence variant located within a transcription factor binding site.
    TfBindingSiteVariant,
    /// A feature ablation whereby the deleted region includes a transcription factor binding site.
    TfbsAblation,
    /// A feature ablation whereby the deleted region includes a transcript feature.
    TranscriptAblation,
    /// A sequence variant located 5' of a gene.
    UpstreamGeneVariant,
}

impl From<Effect> for Impact {
    fn from(value: Effect) -> Impact {
        match value {
            // High Impact
		Effect::BidirectionalGeneFusion
		| Effect::ChromosomeNumberVariation
		| Effect::ExonLossVariant
		| Effect::FeatureAblation
		| Effect::FrameshiftVariant
		| Effect::GeneFusion
		| Effect::ProteinProteinContact
		| Effect::RareAminoAcidVariant
		| Effect::RearrangedAtDnaLevel
		| Effect::SpliceAcceptorVariant
		| Effect::SpliceDonorVariant
		| Effect::StartLost
		| Effect::StopGained
		| Effect::StopLost
		| Effect::StructuralInteractionVariant
		| Effect::TranscriptAblation
		=> Impact::High,
            // Moderate Impact
		Effect::ConservativeInframeDeletion
		| Effect::ConservativeInframeInsertion
		| Effect::DisruptiveInframeDeletion
		| Effect::MissenseVariant
		| Effect::P3PrimeUtrTruncation
		| Effect::P5PrimeUtrTruncation
            | Effect::DisruptiveInframeInsertion
		=> Impact::Moderate,
            // Low Impact
		Effect::InitiatorCodonVariant
		| Effect::P5PrimeUtrPrematureStartCodonGainVariant
		| Effect::SpliceRegionVariant
		| Effect::SynonymousVariant
		| Effect::TfBindingSiteVariant
		| Effect::TfbsAblation
		| Effect::FeatureFusion
		=> Impact::Low,
            // Modifier Impact
		Effect::Chromosome
		| Effect::CodingSequenceVariant
		| Effect::ConservedIntergenicVariant
		| Effect::ConservedIntronVariant
		| Effect::DownstreamGeneVariant
		| Effect::ExonRegion
		| Effect::FeatureElongation
		| Effect::GeneVariant
		| Effect::IntergenicRegion
		| Effect::IntragenicVariant
		| Effect::IntronVariant
		| Effect::NonCodingTranscriptExonVariant
		| Effect::NonCodingTranscriptVariant
		| Effect::P3PrimeUtrVariant
		| Effect::P5PrimeUtrVariant
		| Effect::RegulatoryRegionVariant
		| Effect::SequenceFeature
		| Effect::UpstreamGeneVariant
		=> Impact::Modifier,
            // Other Impact
		Effect::Inversion // Large, Exon -> High | Gene, Transcript -> Moderate
		| Effect::MiRna //
		| Effect::StartRetainedVariant // FRAME_SHIFT_BEFORE_CDS_START -> Modifier | SYNONYMOUS_START -> Low
		| Effect::StopRetainedVariant // FrameShiftAfterCDS -> Modifer | NonSynonymousStop -> Low | SynonymousStop -> Low
		| Effect::Duplication // Large, Exon -> High | Gene, Transcript -> Moderate
		=> Impact::Other,

        }
    }
}

impl From<Effect> for Vec<u8> {
    fn from(value: Effect) -> Vec<u8> {
        match value {
            Effect::BidirectionalGeneFusion => b"bidirectional_gene_fusion".to_vec(),
            Effect::Chromosome => b"chromosome".to_vec(),
            Effect::ChromosomeNumberVariation => b"chromosome_number_variation".to_vec(),
            Effect::CodingSequenceVariant => b"coding_sequence_variant".to_vec(),
            Effect::ConservativeInframeDeletion => b"conservative_inframe_deletion".to_vec(),
            Effect::ConservativeInframeInsertion => b"conservative_inframe_insertion".to_vec(),
            Effect::ConservedIntergenicVariant => b"conserved_intergenic_variant".to_vec(),
            Effect::ConservedIntronVariant => b"conserved_intron_variant".to_vec(),
            Effect::DisruptiveInframeDeletion => b"disruptive_inframe_deletion".to_vec(),
            Effect::DisruptiveInframeInsertion => b"disruptive_inframe_insertion".to_vec(),
            Effect::DownstreamGeneVariant => b"downstream_gene_variant".to_vec(),
            Effect::Duplication => b"duplication".to_vec(),
            Effect::ExonLossVariant => b"exon_loss_variant".to_vec(),
            Effect::ExonRegion => b"exon_region".to_vec(),
            Effect::FeatureAblation => b"feature_ablation".to_vec(),
            Effect::FeatureElongation => b"feature_elongation".to_vec(),
            Effect::FeatureFusion => b"feature_fusion".to_vec(),
            Effect::FrameshiftVariant => b"frameshift_variant".to_vec(),
            Effect::GeneFusion => b"gene_fusion".to_vec(),
            Effect::GeneVariant => b"gene_variant".to_vec(),
            Effect::InitiatorCodonVariant => b"initiator_codon_variant".to_vec(),
            Effect::IntergenicRegion => b"intergenic_region".to_vec(),
            Effect::IntragenicVariant => b"intragenic_variant".to_vec(),
            Effect::IntronVariant => b"intron_variant".to_vec(),
            Effect::Inversion => b"inversion".to_vec(),
            Effect::MiRna => b"miRNA".to_vec(),
            Effect::MissenseVariant => b"missense_variant".to_vec(),
            Effect::NonCodingTranscriptExonVariant => {
                b"non_coding_transcript_exon_variant".to_vec()
            }
            Effect::NonCodingTranscriptVariant => b"non_coding_transcript_variant".to_vec(),
            Effect::P3PrimeUtrTruncation => b"3_prime_UTR_truncation".to_vec(),
            Effect::P3PrimeUtrVariant => b"3_prime_UTR_variant".to_vec(),
            Effect::P5PrimeUtrPrematureStartCodonGainVariant => {
                b"5_prime_UTR_premature_start_codon_gain_variant".to_vec()
            }
            Effect::P5PrimeUtrTruncation => b"5_prime_UTR_truncation".to_vec(),
            Effect::P5PrimeUtrVariant => b"5_prime_UTR_variant".to_vec(),
            Effect::ProteinProteinContact => b"protein_protein_contact".to_vec(),
            Effect::RareAminoAcidVariant => b"rare_amino_acid_variant".to_vec(),
            Effect::RearrangedAtDnaLevel => b"rearranged_at_DNA_level".to_vec(),
            Effect::RegulatoryRegionVariant => b"regulatory_region_variant".to_vec(),
            Effect::SequenceFeature => b"sequence_feature".to_vec(),
            Effect::SpliceAcceptorVariant => b"splice_acceptor_variant".to_vec(),
            Effect::SpliceDonorVariant => b"splice_donor_variant".to_vec(),
            Effect::SpliceRegionVariant => b"splice_region_variant".to_vec(),
            Effect::StartLost => b"start_lost".to_vec(),
            Effect::StartRetainedVariant => b"start_retained_variant".to_vec(),
            Effect::StopGained => b"stop_gained".to_vec(),
            Effect::StopLost => b"stop_lost".to_vec(),
            Effect::StopRetainedVariant => b"stop_retained_variant".to_vec(),
            Effect::StructuralInteractionVariant => b"structural_interaction_variant".to_vec(),
            Effect::SynonymousVariant => b"synonymous_variant".to_vec(),
            Effect::TfBindingSiteVariant => b"TF_binding_site_variant".to_vec(),
            Effect::TfbsAblation => b"TFBS_ablation".to_vec(),
            Effect::TranscriptAblation => b"transcript_ablation".to_vec(),
            Effect::UpstreamGeneVariant => b"upstream_gene_variant".to_vec(),
        }
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */
    use super::*;

    #[test]
    fn impact_order() {
        assert!(Impact::High == Impact::High);
        assert!(Impact::High > Impact::Moderate);
        assert!(Impact::High > Impact::Low);
        assert!(Impact::High > Impact::Modifier);
        assert!(Impact::High > Impact::Other);

        assert!(Impact::Moderate < Impact::High);
        assert!(Impact::Moderate == Impact::Moderate);
        assert!(Impact::Moderate > Impact::Low);
        assert!(Impact::Moderate > Impact::Modifier);
        assert!(Impact::Moderate > Impact::Other);

        assert!(Impact::Low < Impact::High);
        assert!(Impact::Low < Impact::Moderate);
        assert!(Impact::Low == Impact::Low);
        assert!(Impact::Low > Impact::Modifier);
        assert!(Impact::Low > Impact::Other);

        assert!(Impact::Modifier < Impact::High);
        assert!(Impact::Modifier < Impact::Moderate);
        assert!(Impact::Modifier < Impact::Low);
        assert!(Impact::Modifier == Impact::Modifier);
        assert!(Impact::Modifier > Impact::Other);

        assert!(Impact::Other < Impact::High);
        assert!(Impact::Other < Impact::Moderate);
        assert!(Impact::Other < Impact::Low);
        assert!(Impact::Other < Impact::Modifier);
        assert!(Impact::Other == Impact::Other);
    }

    #[test]
    fn impact_max() {
        assert_eq!(
            [
                Impact::High,
                Impact::Moderate,
                Impact::Low,
                Impact::Modifier,
                Impact::Other
            ]
            .into_iter()
            .max(),
            Some(Impact::High)
        );

        assert_eq!(
            [
                Impact::Moderate,
                Impact::Low,
                Impact::Modifier,
                Impact::Other
            ]
            .into_iter()
            .max(),
            Some(Impact::Moderate)
        );

        assert_eq!(
            [Impact::Low, Impact::Modifier, Impact::Other]
                .into_iter()
                .max(),
            Some(Impact::Low)
        );

        assert_eq!(
            [Impact::Modifier, Impact::Other].into_iter().max(),
            Some(Impact::Modifier)
        );

        assert_eq!(
            [Impact::Other, Impact::Other].into_iter().max(),
            Some(Impact::Other)
        );

        assert_eq!(
            [Impact::Modifier, Impact::High].into_iter().max(),
            Some(Impact::High)
        );
    }

    #[test]
    fn effect_to_impact() {
        assert_eq!(Impact::from(Effect::BidirectionalGeneFusion), Impact::High);
        assert_eq!(Impact::from(Effect::Chromosome), Impact::Modifier);
        assert_eq!(
            Impact::from(Effect::ChromosomeNumberVariation),
            Impact::High,
        );
        assert_eq!(
            Impact::from(Effect::CodingSequenceVariant),
            Impact::Modifier
        );
        assert_eq!(
            Impact::from(Effect::ConservativeInframeDeletion),
            Impact::Moderate
        );
        assert_eq!(
            Impact::from(Effect::ConservativeInframeInsertion),
            Impact::Moderate
        );
        assert_eq!(
            Impact::from(Effect::ConservedIntergenicVariant),
            Impact::Modifier
        );
        assert_eq!(
            Impact::from(Effect::ConservedIntronVariant),
            Impact::Modifier
        );
        assert_eq!(
            Impact::from(Effect::DisruptiveInframeDeletion),
            Impact::Moderate
        );
        assert_eq!(
            Impact::from(Effect::DisruptiveInframeInsertion),
            Impact::Moderate
        );
        assert_eq!(
            Impact::from(Effect::DownstreamGeneVariant),
            Impact::Modifier
        );
        assert_eq!(Impact::from(Effect::Duplication), Impact::Other);
        assert_eq!(Impact::from(Effect::ExonLossVariant), Impact::High);
        assert_eq!(Impact::from(Effect::ExonRegion), Impact::Modifier);
        assert_eq!(Impact::from(Effect::FeatureAblation), Impact::High);
        assert_eq!(Impact::from(Effect::FeatureElongation), Impact::Modifier);
        assert_eq!(Impact::from(Effect::FeatureFusion), Impact::Low);
        assert_eq!(Impact::from(Effect::FrameshiftVariant), Impact::High);
        assert_eq!(Impact::from(Effect::GeneFusion), Impact::High);
        assert_eq!(Impact::from(Effect::GeneVariant), Impact::Modifier);
        assert_eq!(Impact::from(Effect::InitiatorCodonVariant), Impact::Low);
        assert_eq!(Impact::from(Effect::IntergenicRegion), Impact::Modifier);
        assert_eq!(Impact::from(Effect::IntragenicVariant), Impact::Modifier);
        assert_eq!(Impact::from(Effect::IntronVariant), Impact::Modifier);
        assert_eq!(Impact::from(Effect::Inversion), Impact::Other);
        assert_eq!(Impact::from(Effect::MiRna), Impact::Other);
        assert_eq!(Impact::from(Effect::MissenseVariant), Impact::Moderate);
        assert_eq!(
            Impact::from(Effect::NonCodingTranscriptExonVariant),
            Impact::Modifier
        );
        assert_eq!(
            Impact::from(Effect::NonCodingTranscriptVariant),
            Impact::Modifier
        );
        assert_eq!(Impact::from(Effect::P3PrimeUtrTruncation), Impact::Moderate);
        assert_eq!(Impact::from(Effect::P3PrimeUtrVariant), Impact::Modifier);
        assert_eq!(
            Impact::from(Effect::P5PrimeUtrPrematureStartCodonGainVariant),
            Impact::Low
        );
        assert_eq!(Impact::from(Effect::P5PrimeUtrTruncation), Impact::Moderate);
        assert_eq!(Impact::from(Effect::P5PrimeUtrVariant), Impact::Modifier);
        assert_eq!(Impact::from(Effect::ProteinProteinContact), Impact::High);
        assert_eq!(Impact::from(Effect::RareAminoAcidVariant), Impact::High);
        assert_eq!(Impact::from(Effect::RearrangedAtDnaLevel), Impact::High);
        assert_eq!(
            Impact::from(Effect::RegulatoryRegionVariant),
            Impact::Modifier
        );
        assert_eq!(Impact::from(Effect::SequenceFeature), Impact::Modifier);
        assert_eq!(Impact::from(Effect::SpliceAcceptorVariant), Impact::High);
        assert_eq!(Impact::from(Effect::SpliceDonorVariant), Impact::High);
        assert_eq!(Impact::from(Effect::SpliceRegionVariant), Impact::Low);
        assert_eq!(Impact::from(Effect::StartLost), Impact::High);
        assert_eq!(Impact::from(Effect::StartRetainedVariant), Impact::Other);
        assert_eq!(Impact::from(Effect::StopGained), Impact::High);
        assert_eq!(Impact::from(Effect::StopLost), Impact::High);
        assert_eq!(Impact::from(Effect::StopRetainedVariant), Impact::Other);
        assert_eq!(
            Impact::from(Effect::StructuralInteractionVariant),
            Impact::High
        );
        assert_eq!(Impact::from(Effect::SynonymousVariant), Impact::Low);
        assert_eq!(Impact::from(Effect::TfBindingSiteVariant), Impact::Low);
        assert_eq!(Impact::from(Effect::TfbsAblation), Impact::Low);
        assert_eq!(Impact::from(Effect::TranscriptAblation), Impact::High);
        assert_eq!(Impact::from(Effect::UpstreamGeneVariant), Impact::Modifier);
    }

    #[test]
    fn effect_to_byte_string() {
        assert_eq!(
            Vec::<u8>::from(Effect::BidirectionalGeneFusion),
            b"bidirectional_gene_fusion".to_vec(),
        );
        assert_eq!(Vec::<u8>::from(Effect::Chromosome), b"chromosome".to_vec());
        assert_eq!(
            Vec::<u8>::from(Effect::ChromosomeNumberVariation),
            b"chromosome_number_variation".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::CodingSequenceVariant),
            b"coding_sequence_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::ConservativeInframeDeletion),
            b"conservative_inframe_deletion".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::ConservativeInframeInsertion),
            b"conservative_inframe_insertion".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::ConservedIntergenicVariant),
            b"conserved_intergenic_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::ConservedIntronVariant),
            b"conserved_intron_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::DisruptiveInframeDeletion),
            b"disruptive_inframe_deletion".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::DisruptiveInframeInsertion),
            b"disruptive_inframe_insertion".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::DownstreamGeneVariant),
            b"downstream_gene_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::Duplication),
            b"duplication".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::ExonLossVariant),
            b"exon_loss_variant".to_vec()
        );
        assert_eq!(Vec::<u8>::from(Effect::ExonRegion), b"exon_region".to_vec());
        assert_eq!(
            Vec::<u8>::from(Effect::FeatureAblation),
            b"feature_ablation".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::FeatureElongation),
            b"feature_elongation".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::FeatureFusion),
            b"feature_fusion".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::FrameshiftVariant),
            b"frameshift_variant".to_vec()
        );
        assert_eq!(Vec::<u8>::from(Effect::GeneFusion), b"gene_fusion".to_vec());
        assert_eq!(
            Vec::<u8>::from(Effect::GeneVariant),
            b"gene_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::InitiatorCodonVariant),
            b"initiator_codon_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::IntergenicRegion),
            b"intergenic_region".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::IntragenicVariant),
            b"intragenic_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::IntronVariant),
            b"intron_variant".to_vec()
        );
        assert_eq!(Vec::<u8>::from(Effect::Inversion), b"inversion".to_vec());
        assert_eq!(Vec::<u8>::from(Effect::MiRna), b"miRNA".to_vec());
        assert_eq!(
            Vec::<u8>::from(Effect::MissenseVariant),
            b"missense_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::NonCodingTranscriptExonVariant),
            b"non_coding_transcript_exon_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::NonCodingTranscriptVariant),
            b"non_coding_transcript_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::P3PrimeUtrTruncation),
            b"3_prime_UTR_truncation".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::P3PrimeUtrVariant),
            b"3_prime_UTR_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::P5PrimeUtrPrematureStartCodonGainVariant),
            b"5_prime_UTR_premature_start_codon_gain_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::P5PrimeUtrTruncation),
            b"5_prime_UTR_truncation".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::P5PrimeUtrVariant),
            b"5_prime_UTR_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::ProteinProteinContact),
            b"protein_protein_contact".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::RareAminoAcidVariant),
            b"rare_amino_acid_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::RearrangedAtDnaLevel),
            b"rearranged_at_DNA_level".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::RegulatoryRegionVariant),
            b"regulatory_region_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::SequenceFeature),
            b"sequence_feature".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::SpliceAcceptorVariant),
            b"splice_acceptor_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::SpliceDonorVariant),
            b"splice_donor_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::SpliceRegionVariant),
            b"splice_region_variant".to_vec()
        );
        assert_eq!(Vec::<u8>::from(Effect::StartLost), b"start_lost".to_vec());
        assert_eq!(
            Vec::<u8>::from(Effect::StartRetainedVariant),
            b"start_retained_variant".to_vec()
        );
        assert_eq!(Vec::<u8>::from(Effect::StopGained), b"stop_gained".to_vec());
        assert_eq!(Vec::<u8>::from(Effect::StopLost), b"stop_lost".to_vec());
        assert_eq!(
            Vec::<u8>::from(Effect::StopRetainedVariant),
            b"stop_retained_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::StructuralInteractionVariant),
            b"structural_interaction_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::SynonymousVariant),
            b"synonymous_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::TfBindingSiteVariant),
            b"TF_binding_site_variant".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::TfbsAblation),
            b"TFBS_ablation".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::TranscriptAblation),
            b"transcript_ablation".to_vec()
        );
        assert_eq!(
            Vec::<u8>::from(Effect::UpstreamGeneVariant),
            b"upstream_gene_variant".to_vec()
        );
    }
}
