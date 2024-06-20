//! Sequence database

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::interval_tree;

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
        interval: &interval_tree::Interval<u64>,
    ) -> Option<&[u8]> {
        self.0
            .get(seqname)
            .map(|seq| &seq[interval.start as usize..interval.end as usize])
    }

    /// Get transcript sequence
    pub fn get_transcript(
        &self,
        seqname: &[u8],
        intervals: &[interval_tree::Interval<u64>],
    ) -> Vec<u8> {
        let mut transcript =
            Vec::with_capacity(intervals.iter().map(|i| (i.end - i.start) as usize).sum());

        for interval in intervals {
            match self.get_interval(seqname, interval) {
                Some(seq) => transcript.extend(seq),
                None => log::error!(
                    "Can't get sequence between position {}..{} in {}",
                    interval.start,
                    interval.end,
                    String::from_utf8(seqname.to_vec()).unwrap(),
                ),
            }
        }

        transcript
    }
}
