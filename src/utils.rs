//! Utility functions for DNA location mapping
//!
//! For now, there is only a CDNAPosition struct, to get the HGVS.c locus description from genomic position.
use std::fmt::Display;

use crate::annotation::{self, Annotation};

use std::slice::Windows;

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

fn intron_length_before_position(exons_coords: &[(u64, u64)], position: u64) -> u64 {
    exons_coords
        .windows(2)
        .filter_map(|x| {
            // Iterating through:
            // |=============EXON=============|--------INTRON---------|=============EXON=============|
            // ^x[0].0                        ^x[0].1                 ^x[1].0                        ^x[1].1

            // Position is beyond this region. So we should count intron size
            if position > x[1].1 {
                Some(x[1].0 - x[0].1)
            } else {
                // Position is within this region

                // Position is within the first exon: do not count
                if position < x[0].1 {
                    None
                } else if position < x[1].0 {
                    Some(position - x[0].1)
                } else {
                    Some(x[1].0 - x[0].1)
                }
            }
        })
        .sum()
}

fn intron_length_after_position(exons_coords: &[(u64, u64)], position: u64) -> u64 {
    exons_coords
        .windows(2)
        .filter_map(|x| {
            // Iterating through:
            // |=============EXON=============|--------INTRON---------|=============EXON=============|
            // ^x[0].0                        ^x[0].1                 ^x[1].0                        ^x[1].1

            // Position is before this region. So we should count intron size
            if position < x[0].0 {
                Some(x[1].0 - x[0].1)
            } else {
                // Position is within this region

                // Position is within the first exon: do not count
                if position < x[0].1 {
                    None
                } else if position < x[1].0 {
                    Some(position - x[0].1)
                } else {
                    Some(x[1].0 - x[0].1)
                }
            }
        })
        .sum()
}

impl CDNAPosition {
    /// Converts a 0-based genomic position into a CDNAGenomicPosition, using the provided annotations slice.
    /// Returns None if either
    pub fn from_genomic_pos(g_pos: u64, annotations: &[&Annotation]) -> Option<Self> {
        if annotations.is_empty() {
            log::debug!("No annotations for genomic position {}", g_pos);
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

        // Get 0-based **pre-mRNA** coordinate of pos, regardless of orientation. Will return early if negative (i.e. 5' upstream if forward, 3' downstream if reverse)
        let tx_pos = match strandness {
            &annotation::Strand::Forward => g_pos as i64 - g_tx_start as i64,
            &annotation::Strand::Reverse => g_tx_start as i64 - g_pos as i64,
        };

        // Check if tx_pos is in bounds. Note that the previous statement already eliminated out of bounds positions
        match strandness {
            &annotation::Strand::Forward if tx_pos > (g_tx_start + tx_len) as i64 => return None,
            &annotation::Strand::Reverse if ((g_tx_start - tx_len) as i64) < tx_pos => return None,
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

        // features coords are now sorted by exon. They are 0-based coordinates, relative to pre-mRNA first nucleotide
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

        // Weird variable to keep track of the distance between current exon end and start codon.
        let mut last_cds_position: i64 = 0;
        let intron_length_before_pos = intron_length_before_position(&exons_coords, tx_pos as u64);

        for (exon_id, exon) in exons_coords.iter().enumerate() {
            // tx_pos is exonic
            if tx_pos >= exon.0 as i64 && tx_pos < exon.1 as i64 {
                // Is it 5'UTR, coding, or 3'UTR ?

                let distance_to_start_codon =
                    exon.1 as i64 - tx_pos as i64 - intron_length_before_pos as i64;

                // Coding
                if tx_pos >= tx_start_codon as i64 && tx_pos < tx_stop_codon as i64 {
                    return Some(CDNAPosition::ExonicCoding {
                        distance_to_start_codon,
                    });
                }
                // 5' UTR
                else if tx_pos < tx_start_codon as i64 {
                    return Some(CDNAPosition::ExonicFivePrimeUTR {
                        distance_to_start_codon,
                    });
                }
                // 3'UTR
                else if tx_pos > tx_stop_codon as i64 {
                    return Some(CDNAPosition::ExonicThreePrimeUTR {
                        distance_to_stop_codon: tx_pos as i64 - tx_stop_codon as i64,
                    });
                }
            }
            // tx_pos is intronic, we just missed it in the previous iteration
            if tx_pos < exon.0 as i64 {
                // Check distance to previous exon, to see which one is closer
                if (exons_coords[exon_id - 1].1.abs_diff(tx_pos as u64))
                    < (exons_coords[exon_id].0.abs_diff(tx_pos as u64))
                {
                    return Some(CDNAPosition::FivePrimeIntronic {
                        last_exon_position: last_cds_position,
                        distance_to_prev_exon: tx_pos as i64 - exons_coords[exon_id - 1].1 as i64
                            + 1,
                    });
                } else {
                    return Some(CDNAPosition::ThreePrimeIntronic {
                        next_exon_position: last_cds_position + 1,
                        distance_to_next_exon: tx_pos as i64 - exon.0 as i64,
                    });
                }
            }
            if tx_pos > tx_len as i64 {
                eprintln!("g.{} is outside transcript!", g_pos);
            }
            // First exon
            last_cds_position = if exon_id == 0 {
                intron_length_before_position(&exons_coords, tx_pos as u64) as i64
                    - tx_start_codon as i64
            } else {
                last_cds_position + (exon.1 - exon.0) as i64
            };
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::annotations_db;

    use super::*;

    use super::CDNAPosition;

    // This GFF is an actual part of ncbiRefSeq, with hg19 coordinates.
    const	GFF_HBB:	&'static	[u8]	=	b"chr11	ncbiRefSeq.2021-05-17	transcript	5246694	5248301	.	-	.	ID=NM_000518.5;Parent=HBB;gene_id=HBB;gene_name=HBB;transcript_id=NM_000518.5
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
    fn test_intron_len_before_position() {
        assert_eq!(
            intron_length_before_position(&vec![(0, 100), (200, 300), (400, 500)], 150),
            50
        );
        assert_eq!(
            intron_length_before_position(&vec![(0, 100), (200, 300), (400, 500)], 250),
            100
        );
        assert_eq!(
            intron_length_before_position(&vec![(0, 100), (200, 300), (400, 500)], 350),
            150
        );
        assert_eq!(
            intron_length_before_position(&vec![(0, 100), (200, 300), (400, 500)], 450),
            200
        );
        assert_eq!(
            intron_length_before_position(&vec![(0, 100), (200, 300), (400, 500)], 101),
            1
        );
        assert_eq!(
            intron_length_before_position(&vec![(0, 100), (200, 300), (400, 500)], 199),
            99
        );
        assert_eq!(
            intron_length_before_position(&vec![(0, 100), (200, 300), (400, 500)], 200),
            100
        );
        assert_eq!(
            intron_length_before_position(&vec![(0, 100), (200, 300), (400, 500)], 501),
            200
        );
    }

    fn test_intron_len_after_position() {
        assert_eq!(
            intron_length_after_position(&vec![(0, 100), (200, 300), (400, 500), (600, 700)], 250),
            200
        );
    }

    #[test]
    fn test_cdna_from_genomic_position() {
        // chr11:5248158A>G should be a splice donor variant with HGVS.c = c.92+2T>C

        let reader: Box<dyn std::io::Read + Send> = Box::new(std::io::Cursor::new(GFF_HBB));
        let annotation_db =
            annotations_db::AnnotationsDataBase::from_reader(std::io::BufReader::new(reader), 100)
                .unwrap();

        let annotations = annotation_db.get_coding_annotation(b"NM_000518.5").unwrap();

        let proxy = annotations
            .iter()
            .collect::<Vec<&crate::annotation::Annotation>>();
        let genomic_positions = [
            5246620, 5246634, 5246642, 5246646, 5246653, 5246662, 5246684, 5246690, 5246717,
            5246718, 5246732, 5246735, 5246737, 5246757, 5246772, 5246772, 5246795, 5246801,
            5246822, 5246822, 5246827, 5246827, 5246829, 5246836, 5246837, 5246840, 5246845,
            5246850, 5246862, 5246863, 5246868, 5246870, 5246873, 5246878, 5246880, 5246883,
            5246886, 5246889, 5246892, 5246892, 5246892, 5246898, 5246901, 5246908, 5246908,
            5246913, 5246929, 5246931, 5246944, 5246945, 5246948, 5246957, 5246958, 5246958,
            5246959, 5246959, 5246971, 5246975, 5246984, 5246989, 5246992, 5246993, 5246998,
            5247001, 5247026, 5247036, 5247052, 5247058, 5247062, 5247070, 5247081, 5247134,
            5247135, 5247141, 5247153, 5247155, 5247194, 5247195, 5247220, 5247234, 5247250,
            5247251, 5247274, 5247720, 5247724, 5247726, 5247733, 5247736, 5247737, 5247781,
            5247791, 5247806, 5247806, 5247808, 5247812, 5247827, 5247828, 5247836, 5247836,
            5247871, 5247876, 5247902, 5247904, 5247914, 5247936, 5247942, 5247985, 5247992,
            5248000, 5248004, 5248008, 5248009, 5248014, 5248024, 5248049, 5248050, 5248052,
            5248055, 5248097, 5248106, 5248107, 5248113, 5248121, 5248143, 5248149, 5248154,
            5248155, 5248155, 5248155, 5248158, 5248158, 5248159, 5248159, 5248160, 5248163,
            5248164, 5248168, 5248170, 5248173, 5248177, 5248182, 5248184, 5248185, 5248191,
            5248192, 5248200, 5248202, 5248205, 5248214, 5248218, 5248219, 5248219, 5248223,
            5248224, 5248232, 5248233, 5248233, 5248243, 5248250, 5248250, 5248259, 5248282,
            5248282, 5248302, 5248313, 5248329, 5248330, 5248331, 5248333, 5248333, 5248343,
            5248352, 5248357, 5248383, 5248384, 5248387, 5248388, 5248388, 5248389, 5248391,
            5248393, 5248402, 5248427, 5248440, 5248440, 5248456, 5248480, 5248487, 5248491,
            5248499, 5248500, 5248502, 5248507, 5248510,
        ];
        let cdna_positions = [
            "*208", "*194", "*186", "*182", "*175", "*166", "*144", "*138", "*111", "*110", "*96",
            "*93", "*91", "*71", "*56", "*56", "*33", "*27", "*6", "*6", "*1", "*1", "443", "436",
            "435", "432", "427", "422", "410", "409", "404", "402", "399", "394", "392", "389",
            "386", "383", "380", "380", "380", "374", "371", "364", "364", "359", "343", "341",
            "328", "327", "324", "316-1", "316-2", "316-2", "316-3", "316-3", "316-15", "316-19",
            "316-28", "316-33", "316-36", "316-37", "316-42", "316-45", "316-70", "316-80",
            "316-96", "316-102", "316-106", "316-114", "316-125", "316-178", "316-179", "316-185",
            "316-197", "316-199", "316-238", "316-239", "316-264", "316-278", "316-294", "316-295",
            "316-318", "315+87", "315+83", "315+81", "315+74", "315+71", "315+70", "315+26",
            "315+16", "315+1", "315+1", "314", "310", "295", "294", "286", "286", "251", "246",
            "220", "218", "208", "186", "180", "137", "130", "122", "118", "114", "113", "108",
            "98", "93-20", "93-21", "93-23", "93-26", "92+63", "92+54", "92+53", "92+47", "92+39",
            "92+17", "92+11", "92+6", "92+5", "92+5", "92+5", "92+2", "92+2", "92+1", "92+1", "92",
            "89", "88", "84", "82", "79", "75", "70", "68", "67", "61", "60", "52", "50", "47",
            "38", "34", "33", "33", "29", "28", "20", "19", "19", "9", "2", "2", "-8", "-31",
            "-31", "-51", "-62", "-78", "-79", "-80", "-82", "-82", "-92", "-101", "-106", "-132",
            "-133", "-136", "-137", "-137", "-138", "-140", "-142", "-151", "-176", "-189", "-189",
            "-205", "-229", "-236", "-240", "-248", "-249", "-251", "-256", "-259", "",
        ];
        for (g_pos, cdna_pos) in genomic_positions.into_iter().zip(cdna_positions.iter()) {
            let cdna_pos_pred = CDNAPosition::from_genomic_pos(g_pos, &proxy).unwrap();
            assert_eq!(
                cdna_pos,
                &cdna_pos_pred.to_string(),
                "Found {:?}",
                &cdna_pos_pred
            );
        }
    }
}
