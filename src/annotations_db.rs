//! Annotations database

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::error;

const DOMAIN_NUMBER: usize = 128;

/// Store annotations information associate to intervals
pub struct AnnotationsDataBase {
    transcripts_intervals: ahash::AHashMap<
        Vec<u8>,
        clairiere::InterpolateTree<u64, annotation::Annotation, DOMAIN_NUMBER>,
    >,
    transcripts2codings: ahash::AHashMap<Vec<u8>, Vec<annotation::Annotation>>,
    transcripts_id2annotation: ahash::AHashMap<Vec<u8>, annotation::Annotation>,
}

impl AnnotationsDataBase {
    /// Build a AnnotationsDataBase from a reader
    pub fn from_reader(
        input: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>,
        updown_distance: u64,
    ) -> error::Result<Self> {
        let mut intervals_builder: ahash::AHashMap<
            Vec<u8>,
            Vec<clairiere::Node<u64, annotation::Annotation>>,
        > = ahash::AHashMap::new();

        let mut transcripts2codings: ahash::AHashMap<Vec<u8>, Vec<annotation::Annotation>> =
            ahash::AHashMap::new();
        let mut transcripts_id2annotation: ahash::AHashMap<Vec<u8>, annotation::Annotation> =
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
                b"exon" | b"start_codon" | b"stop_codon" => {
                    transcripts2codings
                        .entry(annotation.get_parent().to_vec())
                        .or_default()
                        .push(annotation);
                }
                _ => {
                    if annotation.get_feature() == b"transcript" {
                        transcripts_id2annotation.insert(
                            annotation.get_attribute().get_id().to_vec(),
                            annotation.clone(),
                        );
                    }
                    intervals_builder
                        .entry(seqname.to_vec())
                        .and_modify(
                            |tree: &mut Vec<clairiere::Node<u64, annotation::Annotation>>| {
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
            }
        }

        let mut transcripts_intervals: ahash::AHashMap<
            Vec<u8>,
            clairiere::InterpolateTree<u64, annotation::Annotation, DOMAIN_NUMBER>,
        > = ahash::AHashMap::with_capacity(intervals_builder.len());
        for (key, values) in intervals_builder.drain() {
            transcripts_intervals.insert(key, clairiere::InterpolateTree::new(values.clone()));
        }

        Ok(Self {
            transcripts_intervals,
            transcripts2codings,
            transcripts_id2annotation,
        })
    }

    /// Get gene and transcript match with seqname and interval
    pub fn get_annotations(
        &self,
        seqname: &[u8],
        interval: core::ops::Range<u64>,
    ) -> Vec<&annotation::Annotation> {
        if let Some(chr) = self.transcripts_intervals.get(seqname) {
            chr.overlap(interval.start, interval.end)
        } else {
            vec![]
        }
    }

    /// Get coding annotation present in transcript
    pub fn get_coding_annotation(
        &self,
        transcript_id: &[u8],
    ) -> Option<&Vec<annotation::Annotation>> {
        self.transcripts2codings.get(transcript_id)
    }

    /// Get transcript annotation from id
    pub fn get_transcript(&self, transcript_id: &[u8]) -> Option<&annotation::Annotation> {
        self.transcripts_id2annotation.get(transcript_id)
    }

    /// Add annotation
    fn add_annotion(
        tree: &mut Vec<clairiere::Node<u64, annotation::Annotation>>,
        interval: core::ops::Range<u64>,
        annotation: annotation::Annotation,
        updown_distance: u64,
    ) {
        if annotation.get_feature() == b"transcript" {
            let upstream = interval.start.saturating_sub(updown_distance);

            tree.push(clairiere::Node::new(
                upstream,
                interval.start,
                annotation::Annotation::create_child(
                    &annotation,
                    b"upstream",
                    interval.end,
                    interval.end + updown_distance,
                ),
            ));

            tree.push(clairiere::Node::new(
                interval.end,
                interval.end + updown_distance,
                annotation::Annotation::create_child(
                    &annotation,
                    b"downstream",
                    interval.end,
                    interval.end + updown_distance,
                ),
            ));
        }

        tree.push(clairiere::Node::new(
            interval.start,
            interval.end,
            annotation,
        ))
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
    fn annotations() -> error::Result<()> {
        let mut truth = vec![&test_data::GFF_ANNOTATION[0], &test_data::GFF_ANNOTATION[1]];
        truth.sort_by_key(|a| (a.get_start(), a.get_stop()));

        let reader: Box<dyn std::io::Read + Send> = Box::new(test_data::GFF);
        let annotations = AnnotationsDataBase::from_reader(std::io::BufReader::new(reader), 100)?;

        let mut result = annotations.get_annotations(b"chrA", 13250..13251);
        result.sort_by_key(|a| (a.get_start(), a.get_stop()));
        assert_eq!(result, truth);

        let mut truth = test_data::GFF_ANNOTATION[3..8].to_vec();
        truth.sort_by_key(|a| (a.get_start(), a.get_stop()));

        let result = annotations
            .get_coding_annotation(test_data::GFF_ANNOTATION[1].get_attribute().get_id())
            .unwrap();
        let mut value = result.clone();
        value.sort_by_key(|a| (a.get_start(), a.get_stop()));
        assert_eq!(value, truth);

        // seqname not present
        assert_eq!(
            annotations.get_annotations(b"chrX", 2300..2301),
            Vec::<&annotation::Annotation>::new()
        );

        Ok(())
    }
}
