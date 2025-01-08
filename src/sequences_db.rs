//! Sequence database

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::error;
use crate::variant;

/// Perform reverse complement
pub fn rev_comp(seq: &mut [u8]) {
    // Reverse the sequence
    seq.reverse();

    // Complement the sequence
    seq.iter_mut().for_each(|c| {
        if *c & 2 == 0 {
            *c ^= 21;
        } else {
            *c ^= 4;
        }
    });
}

/// Store sequence data
pub struct SequencesDataBase(ahash::AHashMap<Vec<u8>, Vec<u8>>);

impl SequencesDataBase {
    /// Build a SequenceDataBase from a reader
    pub fn from_reader(
        input: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>,
    ) -> error::Result<Self> {
        let mut inner = ahash::AHashMap::new();

        let mut reader = noodles::fasta::io::Reader::new(input);
        let records = reader.records();

        for result in records {
            let record = result?;
            inner.insert(record.name().to_vec(), record.sequence().as_ref().to_vec());
        }

        Ok(Self(inner))
    }

    /// Get interval
    pub fn get_interval(
        &self,
        seqname: &[u8],
        interval: &core::ops::Range<u64>,
    ) -> error::Result<&[u8]> {
        Ok(self
            .0
            .get(seqname)
            .ok_or(error::Error::SeqNotInReferences(unsafe {
                String::from_utf8_unchecked(seqname.to_vec())
            }))?
            .get((interval.start as usize)..(interval.end as usize))
            .ok_or(error::Error::IntervalNotInSeq {
                interval: interval.clone(),
                name: unsafe { String::from_utf8_unchecked(seqname.to_vec()) },
            })?)
    }

    /// Get concatenation of sequence covered by annotations
    pub fn epissed(
        &self,
        annotations: &[&annotation::Annotation],
        strand: annotation::Strand,
    ) -> error::Result<Vec<u8>> {
        if let Some(first) = annotations.first() {
            let seqname = first.get_seqname();

            let mut result = Vec::new();
            for annotation in annotations {
                result.extend(self.get_interval(seqname, &annotation.get_interval())?);
            }

            if strand == annotation::Strand::Reverse {
                rev_comp(&mut result)
            }

            Ok(result)
        } else {
            Ok(vec![])
        }
    }

    /// Get concatenation of sequence covered by annotations edited by variant
    pub fn epissed_edit(
        &self,
        annotations: &[&annotation::Annotation],
        strand: annotation::Strand,
        variant: &variant::Variant,
    ) -> error::Result<Vec<u8>> {
        let epissed = self.epissed(annotations, annotation::Strand::Forward)?;
        let mut edit = Vec::new();

        let mut pos_in_epissed = 0;
        for annotation in annotations {
            if annotation.get_stop() < variant.position {
                pos_in_epissed += annotation.get_stop() - annotation.get_start();
            } else if annotation.get_start() < variant.position
                && annotation.get_stop() > variant.position
            {
                pos_in_epissed += variant.position - annotation.get_start();

                edit.extend(&epissed[..pos_in_epissed as usize]);
                edit.extend(variant.alt_seq.to_vec());
                edit.extend(
                    &epissed[std::cmp::min(
                        pos_in_epissed as usize + variant.ref_seq.len(),
                        epissed.len(),
                    )..],
                );
                break;
            }
        }

        if strand == annotation::Strand::Reverse {
            rev_comp(&mut edit);
        }

        Ok(edit)
    }

    /// Get coding sequence covered by annotations
    pub fn coding(
        &self,
        annotations: &[&annotation::Annotation],
        strand: annotation::Strand,
        mut start_position: std::option::Option<u64>,
        mut stop_position: std::option::Option<u64>,
    ) -> error::Result<Vec<u8>> {
        if let Some(first) = annotations.first() {
            let seqname = first.get_seqname();

            if strand == annotation::Strand::Reverse {
                (start_position, stop_position) = (stop_position, start_position)
            }

            let start = if let Some(start) = start_position {
                start - 1
            } else {
                0
            };
            let stop = if let Some(stop) = stop_position {
                stop - 1
            } else {
                u64::MAX
            };

            let (coding, _, _) =
                self.coding_internal(seqname, annotations, start, stop, u64::MAX)?;

            Self::coding_end(coding, strand)
        } else {
            Ok(vec![])
        }
    }

    /// Get coding sequence covered by annotations edited by variant
    pub fn coding_edit(
        &self,
        annotations: &[&annotation::Annotation],
        strand: annotation::Strand,
        variant: &variant::Variant,
        mut start_position: std::option::Option<u64>,
        mut stop_position: std::option::Option<u64>,
    ) -> error::Result<Vec<u8>> {
        if let Some(first) = annotations.first() {
            let seqname = first.get_seqname();

            if strand == annotation::Strand::Reverse {
                (start_position, stop_position) = (stop_position, start_position)
            }

            let start = if let Some(start) = start_position {
                start - 1
            } else {
                0
            };
            let stop = if let Some(stop) = stop_position {
                stop - 1
            } else {
                u64::MAX
            };

            let (coding, variant_pos, change) =
                self.coding_internal(seqname, annotations, start, stop, variant.position)?;
            if change {
                let mut edit = coding[..variant_pos as usize].to_vec();
                edit.extend(&variant.alt_seq);
                edit.extend(
                    &coding[std::cmp::min(
                        variant_pos as usize + variant.ref_seq.len(),
                        coding.len(),
                    )..],
                );

                Self::coding_end(edit, strand)
            } else {
                Self::coding_end(coding, strand)
            }
        } else {
            Ok(vec![])
        }
    }

    fn coding_internal(
        &self,
        seqname: &[u8],
        annotations: &[&annotation::Annotation],
        start: u64,
        stop: u64,
        mut variant_pos: u64,
    ) -> error::Result<(Vec<u8>, u64, bool)> {
        let mut result = Vec::new();
        let mut in_coding = false;
        for annotation in annotations {
            if annotation.get_stop() < start {
                // before start skip
                continue;
            } else if start < annotation.get_stop()
                && start > annotation.get_start()
                && stop < annotation.get_stop()
                && stop > annotation.get_start()
            {
                // start stop in same exon
                result.extend(self.get_interval(seqname, &(start..stop))?);
                if start < variant_pos && stop > variant_pos {
                    variant_pos -= start;
                    in_coding = true;
                }
                break;
            } else if start < annotation.get_stop() && start > annotation.get_start() {
                // start exon
                result.extend(self.get_interval(seqname, &(start..annotation.get_stop()))?);
                if start < variant_pos && annotation.get_stop() > variant_pos {
                    variant_pos -= start;
                    in_coding = true;
                }
            } else if stop < annotation.get_stop() && stop > annotation.get_start() {
                // stop exon
                result.extend(self.get_interval(seqname, &(annotation.get_start()..stop))?);
                if annotation.get_start() < variant_pos && stop > variant_pos {
                    variant_pos -= annotation.get_start();
                    in_coding = true;
                }
                break;
            } else {
                // all other case
                result.extend(self.get_interval(seqname, &annotation.get_interval())?);
                if annotation.get_start() < variant_pos && stop > variant_pos {
                    variant_pos -= annotation.get_start();
                    in_coding = true;
                }
            }
        }

        Ok((result, variant_pos, in_coding))
    }

    fn coding_end(mut sequence: Vec<u8>, strand: annotation::Strand) -> error::Result<Vec<u8>> {
        if strand == annotation::Strand::Reverse {
            rev_comp(&mut sequence);
        }

        Ok(sequence)
    }
}

#[cfg(test)]
mod tests {
    /* std use */
    use std::io::{Seek, Write as _};

    /* crate use */
    use biotest::Format as _;

    /* project use */
    use super::*;

    #[test]
    fn _rev_comp() -> error::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Sequence::default();

        let mut seq = Vec::new();
        generator.record(&mut seq, &mut rng)?;

        rev_comp(&mut seq);
        assert_eq!(seq, b"tgTCGcAagTTcCAgagGAtGAataCCtgTTAaCGgtAaTTGCaGcaTgtCTCcccttttgcAcaGTACCAGcagaCatGaGCaaccAtCtAtAaTTcgAtaTcCTgcGtacaAgCatTaccgTggctTAACTAACacGCGaTTcATAta".to_vec());

        Ok(())
    }

    fn seqdb_setup() -> error::Result<SequencesDataBase> {
        let mut seq_file = std::io::Cursor::new(Vec::new());

        seq_file.write_all(b">sequence\n")?;
        seq_file.write_all(b"ttcgtctag")?;
        seq_file.write_all(b"ATGACCGCCATGCAAAGGCTCACTGGG")?;
        seq_file.write_all(b"ctctcttcaccc")?;
        seq_file.write_all(b"CTTAAGCATCTACGTATGCGG")?;
        seq_file.write_all(b"gatcgcaggcctctctcggtg")?;
        seq_file.write_all(b"TGTCGTCGGTCGAGGGTTTAACAT")?;
        seq_file.write_all(b"atcctgcttggccaa")?;

        seq_file.seek(std::io::SeekFrom::Start(0))?;

        SequencesDataBase::from_reader(std::io::BufReader::new(
            Box::new(seq_file) as Box<dyn std::io::Read + Send>
        ))
    }

    fn annotation_setup() -> Vec<annotation::Annotation> {
        vec![
            annotation::Annotation::test_annotation(b"sequence".to_vec(), 10, 37),
            annotation::Annotation::test_annotation(b"sequence".to_vec(), 49, 70),
            annotation::Annotation::test_annotation(b"sequence".to_vec(), 91, 115),
        ]
    }

    #[test]
    fn get_interval() -> error::Result<()> {
        let seqdb = seqdb_setup()?;

        assert_eq!(
            b"TGACCGCCATGCAAAGGCTCACTGGGctctcttcacccCT",
            seqdb.get_interval(b"sequence", &(10..50))?
        );

        Ok(())
    }

    #[test]
    fn epissed() -> error::Result<()> {
        let seqdb = seqdb_setup()?;

        let annotations = annotation_setup();

        assert_eq!(
            b"".to_vec(),
            seqdb.epissed(&[], annotation::Strand::Forward)?
        );

        assert_eq!(
            b"ATGACCGCCATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAACAT".to_vec(),
            seqdb.epissed(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward
            )?
        );

        assert_eq!(
            b"ATGTTAAACCCTCGACCGACGACACCGCATACGTAGATGCTTAAGCCCAGTGAGCCTTTGCATGGCGGTCAT".to_vec(),
            seqdb.epissed(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Reverse
            )?
        );

        Ok(())
    }

    #[test]
    fn epissed_edit() -> error::Result<()> {
        let seqdb = seqdb_setup()?;

        let annotations = annotation_setup();

        let variant = variant::Variant::test_variant(b"sequence", 15, b"G", b"ggg");

        assert_eq!(
            b"ATGACCgggCCATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAACAT".to_vec(),
            seqdb.epissed_edit(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward,
                &variant,
            )?
        );

        assert_eq!(
            b"ATGTTAAACCCTCGACCGACGACACCGCATACGTAGATGCTTAAGCCCAGTGAGCCTTTGCATGGcccGGTCAT".to_vec(),
            seqdb.epissed_edit(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Reverse,
                &variant,
            )?
        );

        Ok(())
    }

    #[test]
    fn coding() -> error::Result<()> {
        let seqdb = seqdb_setup()?;
        let annotations = annotation_setup();

        let start_forward = Some(19);
        let stop_forward = Some(112);

        let start_reverse = Some(115);
        let stop_reverse = Some(50);

        assert_eq!(
            b"ATGACCGCCATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAACAT".to_vec(),
            seqdb.coding(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward,
                None,
                None,
            )?
        );

        assert_eq!(
            b"ATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAACAT".to_vec(),
            seqdb.coding(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward,
                start_forward,
                None,
            )?
        );

        assert_eq!(
            b"ATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAA".to_vec(),
            seqdb.coding(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward,
                start_forward,
                stop_forward,
            )?
        );

        assert_eq!(
            b"ATGTTAAACCCTCGACCGACGACACCGCATACGTAGATGCTTAAGCCCAGTGAGCCTTTGCATGGCGGTCAT".to_vec(),
            seqdb.coding(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Reverse,
                None,
                None,
            )?
        );

        assert_eq!(
            b"ATGTTAAACCCTCGACCGACGACACCGCATACGTAGATGCTTAAGCCCAGTGAGCCTTTGCATGGCGGTCAT".to_vec(),
            seqdb.coding(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Reverse,
                start_reverse,
                None,
            )?
        );

        assert_eq!(
            b"ATGTTAAACCCTCGACCGACGACACCGCATACGTAGATGCTTAA".to_vec(),
            seqdb.coding(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Reverse,
                start_reverse,
                stop_reverse,
            )?
        );

        Ok(())
    }

    #[test]
    fn coding_edit() -> error::Result<()> {
        let seqdb = seqdb_setup()?;
        let annotations = annotation_setup();

        let start_forward = Some(19);
        let stop_forward = Some(112);

        let start_reverse = Some(115);
        let stop_reverse = Some(50);

        let forward_before = variant::Variant::test_variant(b"sequence", 15, b"G", b"ggg");
        assert_eq!(
            b"ATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAA".to_vec(),
            seqdb.coding_edit(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward,
                &forward_before,
                start_forward,
                stop_forward,
            )?
        );

        let forward_in = variant::Variant::test_variant(b"sequence", 23, b"A", b"t");
        assert_eq!(
            b"ATGCAtAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAA".to_vec(),
            seqdb.coding_edit(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward,
                &forward_in,
                start_forward,
                stop_forward,
            )?
        );

        let forward_after = variant::Variant::test_variant(b"sequence", 113, b"A", b"g");
        assert_eq!(
            b"ATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAA".to_vec(),
            seqdb.coding_edit(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward,
                &forward_after,
                start_forward,
                stop_forward,
            )?
        );

        let reverse_before = variant::Variant::test_variant(b"sequence", 114, b"A", b"g");
        assert_eq!(
            b"ATGTTAAACCCTCGACCGACGACACCGCATACGTAGATGCTTAA".to_vec(),
            seqdb.coding_edit(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Reverse,
                &reverse_before,
                start_reverse,
                stop_reverse,
            )?
        );

        let reverse_in = variant::Variant::test_variant(b"sequence", 61, b"A", b"g");
        assert_eq!(
            b"ATGTTAAACCCTCGACCGACGACACCGCATAcGTAGATGCTTAA".to_vec(),
            seqdb.coding_edit(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Reverse,
                &reverse_in,
                start_reverse,
                stop_reverse,
            )?
        );

        let reverse_before = variant::Variant::test_variant(b"sequence", 15, b"G", b"ggg");
        assert_eq!(
            b"ATGTTAAACCCTCGACCGACGACACCGCATACGTAGATGCTTAA".to_vec(),
            seqdb.coding_edit(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Reverse,
                &reverse_before,
                start_reverse,
                stop_reverse,
            )?
        );

        Ok(())
    }

    #[test]
    fn single_exon() -> error::Result<()> {
        let seqdb = seqdb_setup()?;

        let annotations = annotation_setup();

        assert_eq!(
            b"ATGACCGCCATGCAAAGGCTCACTGGG".to_vec(),
            seqdb.coding(
                &[&annotations[0]],
                annotation::Strand::Forward,
                //&forward_in,
                None,
                None,
            )?
        );

        let variant = variant::Variant::test_variant(b"sequence", 23, b"A", b"t");
        assert_eq!(
            b"ATGACCGCCATGCAtAGGCTCACTGGG".to_vec(),
            seqdb.coding_edit(
                &[&annotations[0]],
                annotation::Strand::Forward,
                &variant,
                None,
                None,
            )?
        );

        Ok(())
    }
}
