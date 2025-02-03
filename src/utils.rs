//! Utility functions for DNA location mapping
//!
//! For now, there is only a CDNAPosition struct, to get the HGVS.c locus description from genomic position.
use std::fmt::Display;

use crate::annotation::{self, Annotation};

/// Describes a locus position in cDNA coordinates
#[derive(Debug)]
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
    /// Converts a 0-based genomic position into a CDNAGenomicPosition, using the provided annotations slice.
    /// Returns None if either
    pub fn from_genomic_pos(g_pos: u64, annotations: &[&Annotation]) -> Option<Self> {
        if annotations.is_empty() {
            return None;
        }
        //Get strandness once and for all
        let strandness = annotations[0].get_strand();

        // Depending on the orientation, the transcript feature gives us the transcription start site either as first or last genomic coordinate
        let g_tx_start = match strandness {
            &annotation::Strand::Forward => annotations
                .iter()
                .min_by_key(|a| a.get_start())
                .map(|a| a.get_start())?,
            &annotation::Strand::Reverse => annotations
                .iter()
                .max_by_key(|a| a.get_stop())
                .map(|a| a.get_stop())?,
        };

        let tx_len = match strandness {
            &annotation::Strand::Forward => annotations
                .iter()
                .max_by_key(|a| a.get_stop())
                .map(|a| a.get_stop() - g_tx_start)?,
            &annotation::Strand::Reverse => annotations
                .iter()
                .min_by_key(|a| a.get_start())
                .map(|a| g_tx_start - a.get_start())?,
        };

        // Get 0-based pre-mRNA coordinate of pos, regardless of orientation. Will return early if negative (i.e. 5' upstream if forward, 3' downstream if reverse)
        let tx_pos = match strandness {
            &annotation::Strand::Forward => g_pos.checked_sub(g_tx_start)?,
            &annotation::Strand::Reverse => g_tx_start.checked_sub(g_pos)?,
        };

        // Check if tx_pos is in bounds. Note that the previous statement already eliminated out of bounds positions
        match strandness {
            &annotation::Strand::Forward if tx_pos > g_tx_start + tx_len => return None,
            &annotation::Strand::Reverse if g_tx_start - tx_len < tx_pos => return None,
            _ => {}
        }

        // 0-based coordinate, relative to transcript start, of the first start codon base
        let tx_start_codon = {
            let s = annotations
                .iter()
                .find(|a| a.get_feature() == b"start_codon")?;
            match strandness {
                &annotation::Strand::Forward => s.get_start() - g_tx_start,
                &annotation::Strand::Reverse => g_tx_start - s.get_stop(),
            }
        };

        // 0-based coordinate, relative to transcript start, of the last stop codon base
        let tx_stop_codon = {
            let s = annotations
                .iter()
                .find(|a| a.get_feature() == b"stop_codon")?;
            match strandness {
                &annotation::Strand::Forward => s.get_stop() - g_tx_start,
                &annotation::Strand::Reverse => g_tx_start - s.get_start(),
            }
        };

        // features coords are now sorted by exon
        let exons_coords = {
            let mut f: Vec<_> = annotations
                .iter()
                .filter(|ann| ann.get_feature() == b"exon")
                .map(|ann| match strandness {
                    &annotation::Strand::Forward => {
                        (ann.get_start() - g_tx_start, ann.get_stop() - g_tx_start)
                    }
                    &annotation::Strand::Reverse => {
                        (g_tx_start - ann.get_stop(), g_tx_start - ann.get_start())
                    }
                })
                .collect();
            f.sort_by_key(|k| k.0);
            f
        };

        let mut last_mrna_base = 0;
        let mut last_exon_end = 0;
        for exon in exons_coords.iter() {
            // tx_pos is intronic, we just missed it in the previous iteration
            if tx_pos < exon.0 {
                if tx_pos < (exon.0 - last_exon_end).div_euclid(2) {
                    return Some(CDNAPosition::FivePrimeIntronic {
                        last_exon_position: last_mrna_base as i64,
                        distance_to_prev_exon: (tx_pos - last_exon_end) as i64,
                    });
                } else {
                    return Some(CDNAPosition::ThreePrimeIntronic {
                        next_exon_position: (last_mrna_base + 1) as i64,
                        distance_to_next_exon: 0 - (exons_coords[1].0 - tx_pos) as i64,
                    });
                }
            }

            // tx_pos is exonic
            if tx_pos >= exon.0 && tx_pos < exon.1 {
                // Is it 5'UTR, coding, or 3'UTR ?

                // Coding
                if tx_pos >= tx_start_codon && tx_pos < tx_stop_codon {
                    return Some(CDNAPosition::ExonicCoding {
                        distance_to_start_codon: (tx_pos as i64 - tx_start_codon as i64),
                    });
                }
                // 5' UTR
                else if tx_pos < tx_start_codon {
                    return Some(CDNAPosition::ExonicFivePrimeUTR {
                        distance_to_start_codon: (tx_pos as i64 - tx_start_codon as i64),
                    });
                }
                // 3'UTR
                else if tx_pos > tx_stop_codon {
                    return Some(CDNAPosition::ExonicThreePrimeUTR {
                        distance_to_stop_codon: (tx_pos as i64 - tx_stop_codon as i64),
                    });
                }
            }
            last_mrna_base += exon.1 - exon.0;
            last_exon_end = exon.1;
        }
        None
    }
}

#[cfg(test)]
mod tests {
    //TODO: implement tests
}
