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
                edit.extend(&epissed[(pos_in_epissed as usize + variant.ref_seq.len())..]);
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
        start_position: std::option::Option<u64>,
        stop_position: std::option::Option<u64>,
    ) -> error::Result<Vec<u8>> {
        if let Some(first) = annotations.first() {
            let seqname = first.get_seqname();

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

            let mut result = Vec::new();
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
                    break;
                } else if start < annotation.get_stop() && start > annotation.get_start() {
                    // start exon
                    result.extend(self.get_interval(seqname, &(start..annotation.get_stop()))?);
                } else if stop < annotation.get_stop() && stop > annotation.get_start() {
                    // stop exon
                    result.extend(self.get_interval(seqname, &(annotation.get_start()..stop))?);
                    break;
                } else {
                    // all other case
                    result.extend(self.get_interval(seqname, &annotation.get_interval())?)
                }
            }

            if strand == annotation::Strand::Reverse {
                rev_comp(&mut result);
            }

            Ok(result)
        } else {
            Ok(vec![])
        }
    }

    /// Get coding sequence covered by annotations edited by variant
    pub fn coding_edit(
        &self,
        annotations: &[&annotation::Annotation],
        strand: annotation::Strand,
        start_position: std::option::Option<u64>,
        stop_position: std::option::Option<u64>,
        variant: &variant::Variant,
    ) -> error::Result<Vec<u8>> {
        let mut coding = self.coding(
            annotations,
            annotation::Strand::Forward,
            start_position,
            stop_position,
        )?;
        let mut edit = Vec::new();

        let start = if let Some(start) = start_position {
            start - 1
        } else {
            1
        };
        let stop = if let Some(stop) = stop_position {
            stop - 1
        } else {
            u64::MAX
        };

        if variant.position < start || variant.position > stop {
            if strand == annotation::Strand::Reverse {
                rev_comp(&mut coding);
            }

            return Ok(coding);
        }

        let mut pos_in_coding = -(start as i64) + 1;
        for annotation in annotations {
            if annotation.get_stop() < variant.position {
                pos_in_coding += (annotation.get_stop() - annotation.get_start()) as i64;
            } else if annotation.get_start() < variant.position
                && annotation.get_stop() > variant.position
            {
                pos_in_coding += (variant.position - annotation.get_start()) as i64;

                edit.extend(&coding[..pos_in_coding as usize]);
                edit.extend(variant.alt_seq.to_vec());
                edit.extend(&coding[(pos_in_coding as usize + variant.ref_seq.len())..]);
                break;
            }
        }

        if strand == annotation::Strand::Reverse {
            rev_comp(&mut edit);
        }

        Ok(edit)
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */
    use biotest::Format as _;

    use crate::translate::Translate;

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

    fn setup() -> error::Result<SequencesDataBase> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::default();

        let mut temp_input = vec![];
        generator.records(&mut temp_input, &mut rng, 5)?;
        let input: std::io::BufReader<Box<dyn std::io::Read + Send>> =
            std::io::BufReader::new(Box::new(std::io::Cursor::new(temp_input.to_vec())));

        // seqname:
        // - GSWNPZYBHL
        // - RGKPRHQMHK
        // - CVZGQYSRGI
        // - ELUFGTSRIU
        // - RMIVZUTDJN

        SequencesDataBase::from_reader(input)
    }

    #[test]
    fn get_interval() -> error::Result<()> {
        let seqdb = setup()?;

        assert_eq!(
            b"gcGCaCcCGtCtATgTTgTATcaTTCGaCCttcAaGCGCA",
            seqdb.get_interval(b"CVZGQYSRGI", &(10..50))?
        );

        Ok(())
    }

    #[test]
    fn epissed() -> error::Result<()> {
        let seqdb = setup()?;

        let annotations = vec![];

        assert_eq!(
            b"".to_vec(),
            seqdb.epissed(&annotations, annotation::Strand::Forward)?
        );

        let annotations = vec![
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 2, 12),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 15, 32),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 41, 59),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 69, 80),
        ];

        // CttAacGtTtAtGTgACAGCCaCGctGagattTGtgCttaAGggTcCTGcGTAGCTGTCCACgTTTGagtGaGCatAGGACAAaacTaTTagagGtatAGCcTatTtaaaaCGgcttGGTtgaCtgACTacgtCTaTgTCAGgCtaGTtc
        assert_eq!(
            b"ttAacGtTtAgACAGCCaCGctGagatAGggTcCTGcGTAGCTGTgtGaGCatAGG".to_vec(),
            seqdb.epissed(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward
            )?
        );

        assert_eq!(
            b"CCTatGCtCacACAGCTACgCAGgAccCTatctCagCGtGGCTGTcTaAaCgtTaa".to_vec(),
            seqdb.epissed(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Reverse
            )?
        );

        Ok(())
    }

    #[test]
    fn epissed_edit() -> error::Result<()> {
        let seqdb = setup()?;

        let annotations = vec![
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 2, 12),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 15, 32),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 41, 59),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 69, 80),
        ];

        let variant = variant::Variant::test_variant(b"RGKPRHQMHK", 22, b"A", b"Agg");

        assert_eq!(
            b"ttAacGtTtAgACAGCCaAggGctGagatAGggTcCTGcGTAGCTGTgtGaGCatAGG".to_vec(),
            seqdb.epissed_edit(
                &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                annotation::Strand::Forward,
                &variant,
            )?
        );

        assert_eq!(
            b"CCTatGCtCacACAGCTACgCAGgAccCTatctCagCccTtGGCTGTcTaAaCgtTaa".to_vec(),
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
        let seqdb = setup()?;

        let annotations = vec![];

        assert_eq!(
            b"".to_vec(),
            seqdb.epissed(&annotations, annotation::Strand::Forward)?
        );

        let annotations = vec![
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 2, 12),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 15, 32),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 41, 59),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 69, 80),
        ];

        let annotations_ref = annotations.iter().collect::<Vec<&annotation::Annotation>>();

        // no start no stop
        assert_eq!(
            b"ttAacGtTtAgACAGCCaCGctGagatAGggTcCTGcGTAGCTGTgtGaGCatAGG".to_vec(),
            seqdb.coding(&annotations_ref, annotation::Strand::Forward, None, None,)?
        );

        assert_eq!(
            b"CCTatGCtCacACAGCTACgCAGgAccCTatctCagCGtGGCTGTcTaAaCgtTaa".to_vec(),
            seqdb.coding(&annotations_ref, annotation::Strand::Reverse, None, None,)?
        );

        // start in first exon
        assert_eq!(
            b"AacGtTtAgACAGCCaCGctGagatAGggTcCTGcGTAGCTGTgtGaGCatAGG".to_vec(),
            seqdb.coding(&annotations_ref, annotation::Strand::Forward, Some(4), None,)?
        );

        assert_eq!(
            b"CCTatGCtCacACAGCTACgCAGgAccCTatctCagCGtGGCTGTcTaAaCgtT".to_vec(),
            seqdb.coding(&annotations_ref, annotation::Strand::Reverse, Some(4), None,)?
        );

        // start in second exon
        assert_eq!(
            b"AGCCaCGctGagatAGggTcCTGcGTAGCTGTgtGaGCatAGG".to_vec(),
            seqdb.coding(
                &annotations_ref,
                annotation::Strand::Forward,
                Some(18),
                None,
            )?
        );

        assert_eq!(
            b"CCTatGCtCacACAGCTACgCAGgAccCTatctCagCGtGGCT".to_vec(),
            seqdb.coding(
                &annotations_ref,
                annotation::Strand::Reverse,
                Some(18),
                None,
            )?
        );

        // start and stop in second exon
        assert_eq!(
            b"CCaCGctGag".to_vec(),
            seqdb.coding(
                &annotations_ref,
                annotation::Strand::Forward,
                Some(20),
                Some(30),
            )?
        );

        assert_eq!(
            b"ctCagCGtGG".to_vec(),
            seqdb.coding(
                &annotations_ref,
                annotation::Strand::Reverse,
                Some(20),
                Some(30),
            )?
        );

        // stop in last exon
        assert_eq!(
            b"ttAacGtTtAgACAGCCaCGctGagatAGggTcCTGcGTAGCTGTgtGaGC".to_vec(),
            seqdb.coding(
                &annotations_ref,
                annotation::Strand::Forward,
                None,
                Some(75),
            )?
        );

        assert_eq!(
            b"GCtCacACAGCTACgCAGgAccCTatctCagCGtGGCTGTcTaAaCgtTaa".to_vec(),
            seqdb.coding(
                &annotations_ref,
                annotation::Strand::Reverse,
                None,
                Some(75),
            )?
        );

        // stop in penultimate exon
        assert_eq!(
            b"ttAacGtTtAgACAGCCaCGctGagatAGggTcCTG".to_vec(),
            seqdb.coding(
                &annotations_ref,
                annotation::Strand::Forward,
                None,
                Some(50),
            )?
        );

        assert_eq!(
            b"CAGgAccCTatctCagCGtGGCTGTcTaAaCgtTaa".to_vec(),
            seqdb.coding(
                &annotations_ref,
                annotation::Strand::Reverse,
                None,
                Some(50),
            )?
        );

        Ok(())
    }

    #[test]
    fn coding_edit() -> error::Result<()> {
        let seqdb = setup()?;

        let annotations = vec![];

        assert_eq!(
            b"".to_vec(),
            seqdb.epissed(&annotations, annotation::Strand::Forward)?
        );

        let annotations = vec![
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 2, 12),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 15, 32),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 41, 59),
            annotation::Annotation::test_annotation(b"RGKPRHQMHK".to_vec(), 69, 80),
        ];

        let variant = variant::Variant::test_variant(b"RGKPRHQMHK", 22, b"A", b"Agg");

        let annotations_ref = annotations.iter().collect::<Vec<&annotation::Annotation>>();

        // no start no stop
        assert_eq!(
            b"ttAacGtTtAgACAGCCaAggGctGagatAGggTcCTGcGTAGCTGTgtGaGCatAGG".to_vec(),
            seqdb.coding_edit(
                &annotations_ref,
                annotation::Strand::Forward,
                None,
                None,
                &variant
            )?
        );

        assert_eq!(
            b"CCTatGCtCacACAGCTACgCAGgAccCTatctCagCccTtGGCTGTcTaAaCgtTaa".to_vec(),
            seqdb.coding_edit(
                &annotations_ref,
                annotation::Strand::Reverse,
                None,
                None,
                &variant
            )?
        );

        // variant between start stop
        assert_eq!(
            b"AacGtTtAgACAGCCaAggGctGagatAGggTcCTGcGTAGCTGTgtGaGCatAGG".to_vec(),
            seqdb.coding_edit(
                &annotations_ref,
                annotation::Strand::Forward,
                Some(4),
                None,
                &variant
            )?
        );

        assert_eq!(
            b"CCTatGCtCacACAGCTACgCAGgAccCTatctCagCccTtGGCTGTcTaAaCgtT".to_vec(),
            seqdb.coding_edit(
                &annotations_ref,
                annotation::Strand::Reverse,
                Some(4),
                None,
                &variant
            )?
        );

        // variant before start
        assert_eq!(
            b"atAGggTcCTGcGTAGCTGTgtGaGCatAGG".to_vec(),
            seqdb.coding_edit(
                &annotations_ref,
                annotation::Strand::Forward,
                Some(30),
                None,
                &variant
            )?
        );

        assert_eq!(
            b"CCTatGCtCacACAGCTACgCAGgAccCTat".to_vec(),
            seqdb.coding_edit(
                &annotations_ref,
                annotation::Strand::Reverse,
                Some(30),
                None,
                &variant
            )?
        );

        // variant after stop
        assert_eq!(
            b"ttAacGtTtAgACAG".to_vec(),
            seqdb.coding_edit(
                &annotations_ref,
                annotation::Strand::Forward,
                None,
                Some(20),
                &variant
            )?
        );

        assert_eq!(
            b"CTGTcTaAaCgtTaa".to_vec(),
            seqdb.coding_edit(
                &annotations_ref,
                annotation::Strand::Reverse,
                None,
                Some(20),
                &variant
            )?
        );

        Ok(())
    }

    #[test]
    fn not_random() {
        let mut seq_file: Vec<u8> = Vec::new();

        seq_file.extend(b">hbb\n");
        seq_file.extend(b"
ttcgtctagATGACCGCCATGCAAAGGCTCACTGGGctctcttcacccCTTAAGCATCTACGTATGCGGgatcgcaggcctctctcggtgTGTCGTCGGTCGAGGGTTTAAATGatcctgcttggccaa");

        let seqdb = SequencesDataBase::from_reader(std::io::BufReader::new(Box::new(
            std::io::Cursor::new(seq_file),
        )
            as Box<dyn std::io::Read + Send>))
        .unwrap();

        let annotations = vec![
            annotation::Annotation::test_annotation(b"hbb".to_vec(), 10, 37),
            annotation::Annotation::test_annotation(b"hbb".to_vec(), 49, 70),
            annotation::Annotation::test_annotation(b"hbb".to_vec(), 91, 115),
        ];

        assert_eq!(
            b"ATGACCGCCATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAAATG".to_vec(),
            seqdb
                .epissed(
                    &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                    annotation::Strand::Forward
                )
                .unwrap()
        );

        assert_eq!(
            b"ATGACCGCCATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAAATG".to_vec(),
            seqdb
                .coding(
                    &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                    annotation::Strand::Forward,
                    None,
                    None,
                )
                .unwrap()
        );

        assert_eq!(
            b"ATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAAATG".to_vec(),
            seqdb
                .coding(
                    &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                    annotation::Strand::Forward,
                    Some(19),
                    None,
                )
                .unwrap()
        );

        assert_eq!(
            b"ATGCAAAGGCTCACTGGGCTTAAGCATCTACGTATGCGGTGTCGTCGGTCGAGGGTTTAA".to_vec(),
            seqdb
                .coding(
                    &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                    annotation::Strand::Forward,
                    Some(19),
                    Some(112),
                )
                .unwrap()
        );

        let translate = Translate::default();
        assert_eq!(
            b"MQRLTGLKHLRMRCRRSRV*".to_vec(),
            translate.translate(
                &seqdb
                    .coding(
                        &annotations.iter().collect::<Vec<&annotation::Annotation>>(),
                        annotation::Strand::Forward,
                        Some(19),
                        Some(112),
                    )
                    .unwrap()
            )
        );
    }
}
