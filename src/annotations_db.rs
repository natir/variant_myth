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
    pub fn from_reader(
        input: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>,
    ) -> error::Result<Self> {
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
            let seqname = annotation.get_seqname();
            let interval = annotation.get_interval();

            intervals
                .entry(seqname.to_vec())
                .and_modify(
                    |tree: &mut interval_tree::IntervalTree<u64, annotation::Annotation>| {
                        Self::add_annotion(tree, interval.clone(), annotation.clone());
                    },
                )
                .or_insert({
                    let mut tree = interval_tree::IntervalTree::new();

                    Self::add_annotion(&mut tree, interval, annotation);
                    tree
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

    /// Add annotation
    fn add_annotion(
        tree: &mut interval_tree::IntervalTree<u64, annotation::Annotation>,
        interval: interval_tree::Interval<u64>,
        annotation: annotation::Annotation,
    ) {
        if annotation.get_feature() == b"transcript" {
            tree.insert(
                interval.start - 5000..interval.start,
                annotation::Annotation::from_annotation(&annotation, b"upstream"),
            );
            tree.insert(
                interval.end..interval.end + 5000,
                annotation::Annotation::from_annotation(&annotation, b"downstream"),
            );
        }

        tree.insert(interval.clone(), annotation.clone())
    }
}
