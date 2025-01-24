//! Annotations database

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::error;

use superintervals::SuperIntervalsEytz;

/// Store annotations information associate to intervals
pub struct AnnotationsDataBase {
    transcripts_intervals: ahash::AHashMap<Vec<u8>, SuperIntervalsEytz<annotation::Annotation>>,
    transcripts2other: ahash::AHashMap<Vec<u8>, Vec<annotation::Annotation>>,
}

impl AnnotationsDataBase {
    /// Build a AnnotationsDataBase from a reader
    pub fn from_reader(
        input: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>,
        updown_distance: u64,
    ) -> error::Result<Self> {
        let mut transcripts_intervals: ahash::AHashMap<
            Vec<u8>,
            superintervals::SuperIntervalsEytz<annotation::Annotation>,
        > = ahash::AHashMap::new();

        let mut transcripts2other: ahash::AHashMap<Vec<u8>, Vec<annotation::Annotation>> =
            ahash::AHashMap::new();

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

            match annotation.get_feature() {
                b"transcript" | b"gene" => {
                    transcripts_intervals
                        .entry(seqname.to_vec())
                        .and_modify(|tree| {
                            Self::add_annotation(
                                tree,
                                interval.clone(),
                                annotation.clone(),
                                updown_distance,
                            );
                        })
                        .or_insert({
                            let mut tree = SuperIntervalsEytz::new();

                            Self::add_annotation(&mut tree, interval, annotation, updown_distance);
                            tree
                        });
                }
                _ => {
                    if !transcripts2other.contains_key(annotation.get_parent()) {
                        transcripts2other.insert(annotation.get_parent().to_vec(), Vec::new());
                    }
                    transcripts2other
                        .get_mut(annotation.get_parent())
                        .unwrap() // we check previoulsy
                        .push(annotation);
                }
            }
        }
        for (_, value) in transcripts_intervals.iter_mut() {
            value.index();
        }

        Ok(Self {
            transcripts_intervals,
            transcripts2other,
        })
    }

    /// Get gene and transcript match with seqname and interval
    pub fn get_annotations(
        &mut self,
        seqname: &[u8],
        interval: core::ops::Range<u64>,
    ) -> Vec<annotation::Annotation> {
        let mut res: Vec<_> = Vec::new();
        // if let Some(chr) = self.transcripts_intervals.get(seqname) {
        //     chr.find_overlaps(interval.start as i32, interval.end as i32, &mut res);
        // }
        // res
        match self.transcripts_intervals.get_mut(seqname) {
            Some(chr) => {
                chr.find_overlaps(interval.start as i32, interval.end as i32, &mut res);
            }
            None => {}
        }
        res
    }

    /// Get sub annotation present in gene transcript
    pub fn get_subannotations(&self, id: &[u8]) -> Option<&Vec<annotation::Annotation>> {
        self.transcripts2other.get(id)
    }

    /// Add annotation
    fn add_annotation(
        tree: &mut superintervals::SuperIntervalsEytz<annotation::Annotation>,
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

            tree.add(
                upstream as i32,
                interval.start as i32,
                annotation::Annotation::from_annotation(&annotation, b"upstream"),
            );
            tree.add(
                interval.end as i32,
                (interval.end + updown_distance) as i32,
                annotation::Annotation::from_annotation(&annotation, b"downstream"),
            );
        }

        tree.add(interval.start as i32, interval.end as i32, annotation);
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
            String::from_utf8(file[77..197].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        let b3 = Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(file[198..291].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        let b4 = Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(file[292..374].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        let mut truth = vec![b1, b2];
        truth.sort_by_key(|a| (a.get_start(), a.get_stop()));

        let reader: Box<dyn std::io::Read + Send> = Box::new(std::io::Cursor::new(file));
        let mut annotations =
            AnnotationsDataBase::from_reader(std::io::BufReader::new(reader), 100)?;

        let mut result = annotations.get_annotations(b"chr1", 840..841);
        result.sort_by_key(|a| (a.get_start(), a.get_stop()));
        assert_eq!(result, truth);

        let mut truth = vec![b3, b4];
        truth.sort_by_key(|a| (a.get_start(), a.get_stop()));

        let result = annotations
            .get_subannotations(result[1].get_attribute().get_id())
            .unwrap();
        let mut value = result.clone();
        value.sort_by_key(|a| (a.get_start(), a.get_stop()));
        value = value[..2].to_vec();

        assert_eq!(value, truth);

        // seqname not present
        assert_eq!(
            annotations.get_annotations(b"chrX", 2300..2301),
            Vec::<annotation::Annotation>::new()
        );

        Ok(())
    }
}
