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
        /// Distance to start codon (always positive)
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
            } => write!(f, "{}", distance_to_start_codon),
            Self::ExonicCoding {
                distance_to_start_codon,
            } => write!(f, "{}", distance_to_start_codon),
            Self::FivePrimeIntronic {
                last_exon_position,
                distance_to_prev_exon,
            } => write!(f, "{}+{}", last_exon_position, distance_to_prev_exon),
            Self::ThreePrimeIntronic {
                next_exon_position,
                distance_to_next_exon,
            } => write!(f, "{}{}", next_exon_position, distance_to_next_exon),
            Self::ExonicThreePrimeUTR {
                distance_to_stop_codon,
            } => write!(f, "*{}", distance_to_stop_codon),
        }
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
        } + 1;

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

        let mut running_intron_len = 0;
        let mut last_intron_len = 0;
        for (exon_id, exon) in exons_coords.iter().enumerate() {
            if exon_id > 0 {
                last_intron_len = exon.0.abs_diff(exons_coords[exon_id - 1].1);
                running_intron_len += last_intron_len;
            }
            // tx_pos is exonic
            if tx_pos >= exon.0 && tx_pos < exon.1 {
                // Is it 5'UTR, coding, or 3'UTR ?

                // Coding
                if tx_pos >= tx_start_codon && tx_pos < tx_stop_codon {
                    return Some(CDNAPosition::ExonicCoding {
                        distance_to_start_codon: tx_pos as i64
                            - tx_start_codon as i64
                            - running_intron_len as i64,
                    });
                }
                // 5' UTR
                else if tx_pos < tx_start_codon {
                    return Some(CDNAPosition::ExonicFivePrimeUTR {
                        distance_to_start_codon: tx_start_codon as i64
                            - tx_pos as i64
                            - running_intron_len as i64,
                    });
                }
                // 3'UTR
                else if tx_pos > tx_stop_codon {
                    return Some(CDNAPosition::ExonicThreePrimeUTR {
                        distance_to_stop_codon: tx_stop_codon as i64
                            - tx_pos as i64
                            - running_intron_len as i64,
                    });
                }
            }
            // tx_pos is intronic, we just missed it in the previous iteration
            if tx_pos < exon.0 {
                // Check distance to previous exon, to see which one is closer
                if (exons_coords[exon_id - 1].1.abs_diff(tx_pos))
                    < (exons_coords[exon_id].0.abs_diff(tx_pos))
                {
                    return Some(CDNAPosition::FivePrimeIntronic {
                        last_exon_position: exons_coords[exon_id - 1].1 as i64
                            - tx_start_codon as i64
                            - (running_intron_len - last_intron_len) as i64,
                        distance_to_prev_exon: tx_pos as i64 - exons_coords[exon_id - 1].1 as i64,
                    });
                } else {
                    return Some(CDNAPosition::ThreePrimeIntronic {
                        next_exon_position: exon.1 as i64
                            - tx_start_codon as i64
                            - (running_intron_len - last_intron_len) as i64
                            + 1,
                        distance_to_next_exon: tx_pos as i64 - exon.0 as i64,
                    });
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::annotations_db;

    use super::CDNAPosition;

    // This GFF is an actual part of ncbiRefSeq, with hg19 coordinates.
    const	GFF:	&'static	[u8]	=	b"chr11	ncbiRefSeq.2021-05-17	transcript	5246694	5248301	.	-	.	ID=NM_000518.5;Parent=HBB;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	exon	5246694	5246956	.	-	.	ID=NM_000518.5.3;Parent=NM_000518.5;exon_id=NM_000518.5.3;exon_number=3;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	exon	5247807	5248029	.	-	.	ID=NM_000518.5.2;Parent=NM_000518.5;exon_id=NM_000518.5.2;exon_number=2;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	exon	5248160	5248301	.	-	.	ID=NM_000518.5.1;Parent=NM_000518.5;exon_id=NM_000518.5.1;exon_number=1;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	CDS	5246828	5246956	.	-	0	ID=agat-cds-343287;Parent=NM_000518.5;exon_id=NM_000518.5.3;exon_number=3;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	CDS	5247807	5248029	.	-	1	ID=agat-cds-343288;Parent=NM_000518.5;exon_id=NM_000518.5.2;exon_number=2;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	CDS	5248160	5248251	.	-	0	ID=agat-cds-343289;Parent=NM_000518.5;exon_id=NM_000518.5.1;exon_number=1;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	3UTR	5246694	5246827	.	-	.	ID=agat-3utr-35972;Parent=NM_000518.5;exon_id=NM_000518.5.3;exon_number=3;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	5UTR	5248252	5248301	.	-	.	ID=agat-5utr-66412;Parent=NM_000518.5;exon_id=NM_000518.5.1;exon_number=1;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	start_codon	5248249	5248251	.	-	0	ID=agat-start_codon-34632;Parent=NM_000518.5;exon_id=NM_000518.5.1;exon_number=1;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
chr11	ncbiRefSeq.2021-05-17	stop_codon	5246828	5246830	.	-	0	ID=agat-stop_codon-34241;Parent=NM_000518.5;exon_id=NM_000518.5.3;exon_number=3;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5";

    #[test]
    fn test_cdna_from_genomic_position() {
        // chr11:5248158A>G should be a splice donor variant with HGVS.c = c.92+2T>C

        let reader: Box<dyn std::io::Read + Send> = Box::new(std::io::Cursor::new(GFF));
        let annotation_db =
            annotations_db::AnnotationsDataBase::from_reader(std::io::BufReader::new(reader), 100)
                .unwrap();

        let annotations = annotation_db.get_coding_annotation(b"NM_000518.5").unwrap();

        let proxy = annotations
            .iter()
            .collect::<Vec<&crate::annotation::Annotation>>();
        let cdna_pos = CDNAPosition::from_genomic_pos(5247141u64, &proxy).unwrap();
        eprintln!("{:?}", &cdna_pos);

        // chr11:5246958T>C should be a splice acceptor variant
    }
}
