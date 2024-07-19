//! Sequence database

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::error;

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
    pub fn get_interval(&self, seqname: &[u8], interval: &core::ops::Range<u64>) -> Option<&[u8]> {
        self.0
            .get(seqname)
            .map(|seq| &seq[interval.start as usize..interval.end as usize])
    }

    /// Get transcript sequence
    pub fn get_transcript(
        &self,
        seqname: &[u8],
        intervals: &[(core::ops::Range<u64>, annotation::Strand)],
    ) -> Vec<u8> {
        let mut transcript = Vec::with_capacity(
            intervals
                .iter()
                .map(|i| (i.0.end - i.0.start) as usize)
                .sum(),
        );

        for (i, s) in intervals {
            match self.get_interval(seqname, i) {
                Some(seq) => {
                    let len_before = transcript.len();
                    transcript.extend(seq);
                    match s {
                        annotation::Strand::Forward => (),
                        annotation::Strand::Reverse => rev_comp(&mut transcript[len_before..]),
                    }
                }
                None => log::error!(
                    "Can't get sequence between position {}..{} {} in {}",
                    i.start,
                    i.end,
                    s,
                    String::from_utf8(seqname.to_vec()).unwrap(),
                ),
            }
        }

        transcript
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */
    use biotest::Format as _;

    /* project use */
    use super::*;

    #[test]
    fn sequence_db() -> error::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::default();

        let mut temp_input = vec![];
        generator.records(&mut temp_input, &mut rng, 5)?;
        let input: std::io::BufReader<Box<dyn std::io::Read + Send>> =
            std::io::BufReader::new(Box::new(std::io::Cursor::new(temp_input.to_vec())));

        let sequences = SequencesDataBase::from_reader(input)?;

        assert_eq!(
            sequences.get_interval(b"RGKPRHQMHK", &(50..60)).unwrap(),
            b"GTAGCTGTCC".to_vec()
        );

        assert_eq!(
            sequences.get_transcript(
                b"ELUFGTSRIU",
                &vec![
                    (1..5u64, annotation::Strand::Forward),
                    (20..25u64, annotation::Strand::Reverse),
                    (32..35u64, annotation::Strand::Forward)
                ]
            ),
            b"taTgTTgccTct".to_vec()
        );

        assert_eq!(
            sequences.get_transcript(
                b"chr1",
                &vec![
                    (1..5u64, annotation::Strand::Forward),
                    (20..25u64, annotation::Strand::Reverse),
                    (32..35u64, annotation::Strand::Forward)
                ]
            ),
            b"".to_vec()
        );

        Ok(())
    }
}
