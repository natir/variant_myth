//! Annotations database

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::error;
use crate::interval_tree;

/// Store annotations information associate to intervals
pub struct AnnotationsDataBase {
    intervals: ahash::AHashMap<Vec<u8>, interval_tree::IntervalTree<u64, annotation::Annotation>>,
}

impl AnnotationsDataBase {
    /// Build a AnnotationsDataBase from a reader
    pub fn from_reader(input: std::io::BufReader<Box<dyn std::io::Read>>) -> error::Result<Self> {
        let mut intervals: ahash::AHashMap<
            Vec<u8>,
            interval_tree::IntervalTree<u64, annotation::Annotation>,
        > = ahash::AHashMap::new();

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(input);

        for result in reader.byte_records() {
            let annotation = annotation::Annotation::from_byte_record(&result?)?;
            let (seqname, interval) = annotation.get_interval();

            intervals
                .entry(seqname.to_vec())
                .and_modify(
                    |e: &mut interval_tree::IntervalTree<u64, annotation::Annotation>| {
                        e.insert(interval.clone(), annotation.clone())
                    },
                )
                .or_insert({
                    let mut i = interval_tree::IntervalTree::new();

                    i.insert(interval, annotation);
                    i
                });
        }

        intervals.values_mut().for_each(|i| i.index());

        Ok(Self { intervals })
    }

    /// Get annotation match with seqname and interval
    pub fn get_annotation(
        &self,
        seqname: &[u8],
        interval: interval_tree::Interval<u64>,
    ) -> Vec<&annotation::Annotation> {
        if let Some(chr) = self.intervals.get(seqname) {
            chr.find(interval).into_iter().map(|e| e.data()).collect()
        } else {
            vec![]
        }
    }
}
