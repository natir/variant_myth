//! Annotations database

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::error;

const DOMAIN_NUMBER: usize = 8092;

/// Store annotations information associate to intervals
pub struct AnnotationsDataBase {
    intervals: ahash::AHashMap<Vec<u8>, iitiiri::Iitii<u64, annotation::Annotation, DOMAIN_NUMBER>>,
}

impl AnnotationsDataBase {
    /// Build a AnnotationsDataBase from a reader
    pub fn from_reader(
        input: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>,
        updown_distance: u64,
    ) -> error::Result<Self> {
        let mut intervals_builder: ahash::AHashMap<
            Vec<u8>,
            Vec<iitiiri::Node<u64, annotation::Annotation>>,
        > = ahash::AHashMap::new();

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(input);

        for result in reader.byte_records() {
            let annotation = match annotation::Annotation::from_byte_record(&result?) {
                Ok(annot) => annot,
                Err(error) => {
                    log::error!("{}", error);
                    continue;
                }
            };

            let seqname = annotation.get_seqname();
            let interval = annotation.get_interval();

            intervals_builder
                .entry(seqname.to_vec())
                .and_modify(
                    |tree: &mut Vec<iitiiri::Node<u64, annotation::Annotation>>| {
                        Self::add_annotion(
                            tree,
                            interval.clone(),
                            annotation.clone(),
                            updown_distance,
                        );
                    },
                )
                .or_insert({
                    let mut tree = Vec::new();

                    Self::add_annotion(&mut tree, interval, annotation, updown_distance);
                    tree
                });
        }

        let mut intervals: ahash::AHashMap<
            Vec<u8>,
            iitiiri::Iitii<u64, annotation::Annotation, DOMAIN_NUMBER>,
        > = ahash::AHashMap::with_capacity(intervals_builder.len());
        for (key, values) in intervals_builder.drain() {
            intervals.insert(key, iitiiri::Iitii::new(values.clone()));
        }

        Ok(Self { intervals })
    }

    /// Get annotation match with seqname and interval
    pub fn get_annotation(
        &self,
        seqname: &[u8],
        interval: core::ops::Range<u64>,
    ) -> Vec<&annotation::Annotation> {
        if let Some(chr) = self.intervals.get(seqname) {
            chr.overlap(interval.start, interval.end)
        } else {
            vec![]
        }
    }

    /// Add annotation
    fn add_annotion(
        tree: &mut Vec<iitiiri::Node<u64, annotation::Annotation>>,
        interval: core::ops::Range<u64>,
        annotation: annotation::Annotation,
        updown_distance: u64,
    ) {
        if annotation.get_feature() == b"transcript" {
            let upstream = if interval.start < updown_distance {
                0
            } else {
                interval.start - updown_distance
            };

            tree.push(iitiiri::Node::new(
                upstream,
                interval.start,
                annotation::Annotation::from_annotation(&annotation, b"upstream"),
            ));
            tree.push(iitiiri::Node::new(
                interval.end,
                interval.end + updown_distance,
                annotation::Annotation::from_annotation(&annotation, b"downstream"),
            ));
        }

        tree.push(iitiiri::Node::new(interval.start, interval.end, annotation))
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */
    use bstr::ByteSlice as _;

    use crate::annotation::Annotation;

    /* project use */
    use super::*;
    use crate::tests::GFF;

    #[test]
    fn annotations() -> error::Result<()> {
        let file = GFF.replace(b"{0}", b"chr1");

        let b1 = Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(file[..76].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        let b2 = Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(file[77..191].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        let b3 = Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(file[192..285].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        let b4 = Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(file[286..368].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        let mut truth = vec![&b1, &b2, &b3, &b4];
        truth.sort_by_key(|a| (a.get_start(), a.get_stop()));
        for annot in &truth {
            println!("truth: {}", annot);
        }

        let reader: Box<dyn std::io::Read + Send> = Box::new(std::io::Cursor::new(file));
        let annotations = AnnotationsDataBase::from_reader(std::io::BufReader::new(reader), 100)?;

        let mut result = annotations.get_annotation(b"chr1", 840..841);
        result.sort_by_key(|a| (a.get_start(), a.get_stop()));
        for annot in &result {
            println!("result: {}", annot);
        }
        assert_eq!(result, truth);

        // seqname not present
        assert_eq!(
            annotations.get_annotation(b"chrX", 2300..2301),
            Vec::<&annotation::Annotation>::new()
        );

        Ok(())
    }
}
