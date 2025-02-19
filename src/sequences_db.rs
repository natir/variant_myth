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
                if annotation.get_start() < variant_pos && annotation.get_stop() > variant_pos {
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

    /* crate use */

    /* project use */
    use super::*;
    use crate::test_data;

    #[test]
    fn _rev_comp() -> error::Result<()> {
        let mut seq = test_data::SEQUENCE_DB
            .get_interval(b"chrA", &(231..291))?
            .to_vec();

        assert_eq!(
            seq,
            b"ACCCACATAGTTCAATTTCAATATACGAAGgtaggcactgagatttcaatatcagaagaa".to_vec()
        );
        rev_comp(&mut seq);
        assert_eq!(
            seq,
            b"ttcttctgatattgaaatctcagtgcctacCTTCGTATATTGAAATTGAACTATGTGGGT".to_vec()
        );

        Ok(())
    }

    #[test]
    fn get_interval() -> error::Result<()> {
        assert_eq!(
            b"TTCAATTTCAATATACGAAGgtaggcactgagatttcaat",
            test_data::SEQUENCE_DB.get_interval(b"chrA", &(241..281))?
        );

        Ok(())
    }

    #[test]
    fn epissed() -> error::Result<()> {
        let annotations = test_data::GFF_ANNOTATION[3..5]
            .iter()
            .collect::<Vec<&annotation::Annotation>>();

        assert_eq!(
            b"".to_vec(),
            test_data::SEQUENCE_DB.epissed(&[], annotation::Strand::Forward)?
        );

        assert_eq!(
            b"AGCTGACTTAAGAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAGTTCTTCAGACGACGGTCCCTGAGTTACTGAAGCCACTTCACCTGTTTGGGCAGACAGCTGGGAGTGCCCAGAGCTGACACCCTCCAGGTGACCCACAGGTAACGGCTGACCCACGCTGGAGTGTAGGAGCCTTGCTTCAAGACCACACAGACTTTGAG".to_vec(),
            test_data::SEQUENCE_DB.epissed(&annotations, annotation::Strand::Forward)?
        );

        assert_eq!(
            b"CTCAAAGTCTGTGTGGTCTTGAAGCAAGGCTCCTACACTCCAGCGTGGGTCAGCCGTTACCTGTGGGTCACCTGGAGGGTGTCAGCTCTGGGCACTCCCAGCTGTCTGCCCAAACAGGTGAAGTGGCTTCAGTAACTCAGGGACCGTCGTCTGAAGAACTTCGTATATTGAAATTGAACTATGTGGGTATTCATCAGGCCTCTTTCATTTGCAAGAGGCAGAAATTCAACTCAGGAGGGCTTCAGTAATGCAAGTTCAGCCATGGTGGGATTCAGGGGCTTAATGTCACTGGAACTTGGCGTCTCCCTCCATGTCTCCATTTTGCTTTCCTCTGCGTTGAGTTCCTTCTTAAGTCAGCT".to_vec(),
            test_data::SEQUENCE_DB.epissed(&annotations, annotation::Strand::Reverse)?
        );

        Ok(())
    }

    #[test]
    fn epissed_edit() -> error::Result<()> {
        let annotations = test_data::GFF_ANNOTATION[3..5]
            .iter()
            .collect::<Vec<&annotation::Annotation>>();

        let variant = variant::Variant::test_variant(b"chrA", 61, b"G", b"ggg", None)?;

        assert_eq!(
            b"AgggCTGACTTAAGAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAGTTCTTCAGACGACGGTCCCTGAGTTACTGAAGCCACTTCACCTGTTTGGGCAGACAGCTGGGAGTGCCCAGAGCTGACACCCTCCAGGTGACCCACAGGTAACGGCTGACCCACGCTGGAGTGTAGGAGCCTTGCTTCAAGACCACACAGACTTTGAG".to_vec(),
            test_data::SEQUENCE_DB.epissed_edit(&annotations, annotation::Strand::Forward, &variant)?
        );

        assert_eq!(
            b"CTCAAAGTCTGTGTGGTCTTGAAGCAAGGCTCCTACACTCCAGCGTGGGTCAGCCGTTACCTGTGGGTCACCTGGAGGGTGTCAGCTCTGGGCACTCCCAGCTGTCTGCCCAAACAGGTGAAGTGGCTTCAGTAACTCAGGGACCGTCGTCTGAAGAACTTCGTATATTGAAATTGAACTATGTGGGTATTCATCAGGCCTCTTTCATTTGCAAGAGGCAGAAATTCAACTCAGGAGGGCTTCAGTAATGCAAGTTCAGCCATGGTGGGATTCAGGGGCTTAATGTCACTGGAACTTGGCGTCTCCCTCCATGTCTCCATTTTGCTTTCCTCTGCGTTGAGTTCCTTCTTAAGTCAGcccT".to_vec(),
            test_data::SEQUENCE_DB.epissed_edit(&annotations, annotation::Strand::Reverse, &variant)?
        );

        Ok(())
    }

    #[test]
    fn coding() -> error::Result<()> {
        let annotations = test_data::GFF_ANNOTATION[3..5]
            .iter()
            .collect::<Vec<&annotation::Annotation>>();

        let start_forward = Some(63);
        let stop_forward = Some(13324);

        let start_reverse = Some(13340);
        let stop_reverse = Some(220);

        assert_eq!(
            b"AGCTGACTTAAGAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAGTTCTTCAGACGACGGTCCCTGAGTTACTGAAGCCACTTCACCTGTTTGGGCAGACAGCTGGGAGTGCCCAGAGCTGACACCCTCCAGGTGACCCACAGGTAACGGCTGACCCACGCTGGAGTGTAGGAGCCTTGCTTCAAGACCACACAGACTTTGAG".to_vec(),
            test_data::SEQUENCE_DB.coding(&annotations, annotation::Strand::Forward, None, None,)?
        );

        assert_eq!(
            b"CTGACTTAAGAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAGTTCTTCAGACGACGGTCCCTGAGTTACTGAAGCCACTTCACCTGTTTGGGCAGACAGCTGGGAGTGCCCAGAGCTGACACCCTCCAGGTGACCCACAGGTAACGGCTGACCCACGCTGGAGTGTAGGAGCCTTGCTTCAAGACCACACAGACTTTGAG".to_vec(),
            test_data::SEQUENCE_DB.coding(
                &annotations,
                annotation::Strand::Forward,
                start_forward,
                None,
            )?
        );

        assert_eq!(
            b"CTGACTTAAGAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAGTTCTTCAGACGACGGTCCCTGAGTTACTGAAGCCACTTCACCTGTTTGGGCAGACAGCTGGGAGTGCCCAGAGCTGACACCCTCCAGGTGACCCACAGGTAACGGCTGACCCACGCTGGAGT".to_vec(),
            test_data::SEQUENCE_DB.coding(
                &annotations,
                annotation::Strand::Forward,
                start_forward,
                stop_forward,
            )?
        );

        assert_eq!(
            b"CTCAAAGTCTGTGTGGTCTTGAAGCAAGGCTCCTACACTCCAGCGTGGGTCAGCCGTTACCTGTGGGTCACCTGGAGGGTGTCAGCTCTGGGCACTCCCAGCTGTCTGCCCAAACAGGTGAAGTGGCTTCAGTAACTCAGGGACCGTCGTCTGAAGAACTTCGTATATTGAAATTGAACTATGTGGGTATTCATCAGGCCTCTTTCATTTGCAAGAGGCAGAAATTCAACTCAGGAGGGCTTCAGTAATGCAAGTTCAGCCATGGTGGGATTCAGGGGCTTAATGTCACTGGAACTTGGCGTCTCCCTCCATGTCTCCATTTTGCTTTCCTCTGCGTTGAGTTCCTTCTTAAGTCAGCT".to_vec(),
            test_data::SEQUENCE_DB.coding(&annotations, annotation::Strand::Reverse, None, None,)?
        );

        assert_eq!(
            b"GAAGCAAGGCTCCTACACTCCAGCGTGGGTCAGCCGTTACCTGTGGGTCACCTGGAGGGTGTCAGCTCTGGGCACTCCCAGCTGTCTGCCCAAACAGGTGAAGTGGCTTCAGTAACTCAGGGACCGTCGTCTGAAGAACTTCGTATATTGAAATTGAACTATGTGGGTATTCATCAGGCCTCTTTCATTTGCAAGAGGCAGAAATTCAACTCAGGAGGGCTTCAGTAATGCAAGTTCAGCCATGGTGGGATTCAGGGGCTTAATGTCACTGGAACTTGGCGTCTCCCTCCATGTCTCCATTTTGCTTTCCTCTGCGTTGAGTTCCTTCTTAAGTCAGCT".to_vec(),
            test_data::SEQUENCE_DB.coding(
                &annotations,
                annotation::Strand::Reverse,
                start_reverse,
                None,
            )?
        );

        assert_eq!(
            b"GAAGCAAGGCTCCTACACTCCAGCGTGGGTCAGCCGTTACCTGTGGGTCACCTGGAGGGTGTCAGCTCTGGGCACTCCCAGCTGTCTGCCCAAACAGGTGAAGTGGCTTCAGTAACTCAGGGACCGTCGTCTGAAGAACTTCGTATATTGAAATTGAACTATGTGGGTATTCATCAGGCC".to_vec(),
            test_data::SEQUENCE_DB.coding(
                &annotations,
                annotation::Strand::Reverse,
                start_reverse,
                stop_reverse,
            )?
        );

        Ok(())
    }

    #[test]
    fn coding_edit() -> error::Result<()> {
        let annotations = test_data::GFF_ANNOTATION[3..5]
            .iter()
            .collect::<Vec<&annotation::Annotation>>();

        let start_forward = Some(63);
        let stop_forward = Some(13324);

        let start_reverse = Some(13340);
        let stop_reverse = Some(220);

        let forward_before = variant::Variant::test_variant(b"sequence", 62, b"G", b"ggg", None)?;
        assert_eq!(
            b"CTGACTTAAGAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAGTTCTTCAGACGACGGTCCCTGAGTTACTGAAGCCACTTCACCTGTTTGGGCAGACAGCTGGGAGTGCCCAGAGCTGACACCCTCCAGGTGACCCACAGGTAACGGCTGACCCACGCTGGAGT".to_vec(),
            test_data::SEQUENCE_DB.coding_edit(
                &annotations,
                annotation::Strand::Forward,
                &forward_before,
                start_forward,
                stop_forward,
            )?
        );

        let forward_in = variant::Variant::test_variant(b"sequence", 71, b"G", b"ggg", None)?;
        assert_eq!(
            b"CTGACTTAAgggAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAGTTCTTCAGACGACGGTCCCTGAGTTACTGAAGCCACTTCACCTGTTTGGGCAGACAGCTGGGAGTGCCCAGAGCTGACACCCTCCAGGTGACCCACAGGTAACGGCTGACCCACGCTGGAGT".to_vec(),
            test_data::SEQUENCE_DB.coding_edit(
                &annotations,
                annotation::Strand::Forward,
                &forward_in,
                start_forward,
                stop_forward,
            )?
        );

        let forward_after = variant::Variant::test_variant(b"sequence", 13341, b"A", b"g", None)?;
        assert_eq!(
            b"CTGACTTAAGAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAGTTCTTCAGACGACGGTCCCTGAGTTACTGAAGCCACTTCACCTGTTTGGGCAGACAGCTGGGAGTGCCCAGAGCTGACACCCTCCAGGTGACCCACAGGTAACGGCTGACCCACGCTGGAGT".to_vec(),
            test_data::SEQUENCE_DB.coding_edit(
                &annotations,
                annotation::Strand::Forward,
                &forward_after,
                start_forward,
                stop_forward,
            )?
        );

        let reverse_before = variant::Variant::test_variant(b"sequence", 13350, b"A", b"g", None)?;
        assert_eq!(
            b"GAAGCAAGGCTCCTACACTCCAGCGTGGGTCAGCCGTTACCTGTGGGTCACCTGGAGGGTGTCAGCTCTGGGCACTCCCAGCTGTCTGCCCAAACAGGTGAAGTGGCTTCAGTAACTCAGGGACCGTCGTCTGAAGAACTTCGTATATTGAAATTGAACTATGTGGGTATTCATCAGGCC".to_vec(),
            test_data::SEQUENCE_DB.coding_edit(
                &annotations,
                annotation::Strand::Reverse,
                &reverse_before,
                start_reverse,
                stop_reverse,
            )?
        );

        let reverse_in = variant::Variant::test_variant(b"sequence", 13292, b"C", b"aaa", None)?;
        assert_eq!(
            b"GAAGCAAGGCTCCTACACTCCAGCGTGGGTCAGCCGTTACCTGTGGGTCACCTGGAGGGTGTCAGCTCTGGGCACTCCCAGCTGTCTGtttCCAAACAGGTGAAGTGGCTTCAGTAACTCAGGGACCGTCGTCTGAAGAACTTCGTATATTGAAATTGAACTATGTGGGTATTCATCAGGCC".to_vec(),
            test_data::SEQUENCE_DB.coding_edit(
                &annotations,
                annotation::Strand::Reverse,
                &reverse_in,
                start_reverse,
                stop_reverse,
            )?
        );

        let reverse_before = variant::Variant::test_variant(b"sequence", 70, b"G", b"ggg", None)?;
        assert_eq!(
            b"GAAGCAAGGCTCCTACACTCCAGCGTGGGTCAGCCGTTACCTGTGGGTCACCTGGAGGGTGTCAGCTCTGGGCACTCCCAGCTGTCTGCCCAAACAGGTGAAGTGGCTTCAGTAACTCAGGGACCGTCGTCTGAAGAACTTCGTATATTGAAATTGAACTATGTGGGTATTCATCAGGCC".to_vec(),
            test_data::SEQUENCE_DB.coding_edit(
                &annotations,
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
        let annotations = &test_data::GFF_ANNOTATION[3];

        assert_eq!(
            b"AGCTGACTTAAGAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAG".to_vec(),
            test_data::SEQUENCE_DB.coding(&[annotations], annotation::Strand::Forward, None, None,)?
        );

        let variant = variant::Variant::test_variant(b"chrA", 61, b"G", b"ggg", None)?;
        assert_eq!(
            b"AgggCTGACTTAAGAAGGAACTCAACGCAGAGGAAAGCAAAATGGAGACATGGAGGGAGACGCCAAGTTCCAGTGACATTAAGCCCCTGAATCCCACCATGGCTGAACTTGCATTACTGAAGCCCTCCTGAGTTGAATTTCTGCCTCTTGCAAATGAAAGAGGCCTGATGAATACCCACATAGTTCAATTTCAATATACGAAG".to_vec(),
            test_data::SEQUENCE_DB.coding_edit(
                &[annotations],
                annotation::Strand::Forward,
                &variant,
                None,
                None,
            )?
        );

        Ok(())
    }
}
