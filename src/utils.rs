use std::{
    fmt::Display,
    ops::{Range, RangeBounds},
};

use crate::annotation::{self, Annotation};

/// Describes a locus position in cDNA coordinates
pub enum CDNAPosition {
    /// The locus is in an exonic region, in the 5'UTR region
    ExonicFivePrimeUTR {
        /// Distance to start codon, always negative
        distance_to_start_codon: i64,
    },

    /// The locus is in an exonic region, and within the coding sequence
    ExonicCoding {
        /// Distance to start codon (negative means 5', positive means 3', never zero!)
        distance_to_start_codon: i64,
    },

    /// The locus is in an intronic region, and the closest exon is 5' of this location
    FivePrimeIntronic {
        /// Distance from previous exon end to start codon
        last_exon_position: i64,
        /// Distance from previous exon end to this locus (always positive)
        distance_to_prev_exon: i64,
    },
    /// The locus is in an intronic region, and the closest exon is 3' of this location
    ThreePrimeIntronic {
        /// Distance from next exon start to start codon
        next_exon_position: i64,
        /// Distance from next exon start to this locus (always negative)
        distance_to_next_exon: i64,
    },

    /// The locus is in an exonic region, in the 3'UTR region
    ExonicThreePrimeUTR {
        /// Distance to stop codon (always positive)
        distance_to_stop_codon: i64,
    },
}

impl Display for CDNAPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ExonicFivePrimeUTR {
                distance_to_start_codon,
            } => {
                write!(f, "{}", distance_to_start_codon)?;
            }
            Self::ExonicCoding {
                distance_to_start_codon,
            } => {
                write!(f, "{}", distance_to_start_codon)?;
            }
            Self::FivePrimeIntronic {
                last_exon_position,
                distance_to_prev_exon,
            } => {
                write!(f, "{}+{}", last_exon_position, distance_to_prev_exon)?;
            }
            Self::ThreePrimeIntronic {
                next_exon_position,
                distance_to_next_exon,
            } => {
                write!(f, "{}{}", next_exon_position, distance_to_next_exon)?;
            }
            Self::ExonicThreePrimeUTR {
                distance_to_stop_codon,
            } => {
                write!(f, "*{}", distance_to_stop_codon)?;
            }
        }
        Ok(())
    }
}

impl CDNAPosition {
    pub fn from_genomic_pos(pos: u64, annotations: &[&Annotation]) -> Option<Self> {
        if annotations.is_empty() {
            return None;
        }
        //Get strandness once and for all
        let strandness = annotations[0].get_strand();

        // start_ann is the annotation feature describing TSS
        let start_ann = annotations
            .iter()
            .find(|a| a.get_feature() == b"five_prime_utr")?;

        // Depending on the orientation, the 5'UTR feature gives us the transcription start site either as first or last genomic coordinate
        let tss = match strandness {
            &annotation::Strand::Forward => start_ann.get_start(),
            &annotation::Strand::Reverse => start_ann.get_stop(),
        };

        // Get 0-based pre-mRNA coordinate, regardless of orientation
        let tx_pos = match strandness {
            &annotation::Strand::Forward => pos - tss,
            &annotation::Strand::Reverse => tss - pos,
        };
        let features_coords = {
            let mut f: Vec<_> = annotations
                .iter()
                .map(|ann| match strandness {
                    &annotation::Strand::Forward => ann.get_start() - tss..ann.get_stop() - tss,
                    &annotation::Strand::Reverse => tss - ann.get_stop()..tss - ann.get_start(),
                })
                .collect();
            f.sort_by_key(|k| k.start);
        };
        None
    }
}

#[cfg(test)]
mod tests {
    use std::default;

    use crate::annotation;

    #[test]
    fn test_genomic_to_cdna_pos() {
        let genomic = 270i64;
    }
}
