//! Annotations database

/* std use */

/* crate use */
use rust_lapper;

/* project use */
use crate::annotation;
use crate::error;

/// Store annotations information associate to intervals
pub struct AnnotationsDataBase {
    intervals: ahash::AHashMap<Vec<u8>, rust_lapper::Lapper<u64, annotation::Annotation>>,
}

impl AnnotationsDataBase {
    /// Build a AnnotationsDataBase from a reader
    pub fn from_reader(
        input: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>,
        updown_distance: u64,
    ) -> error::Result<Self> {
        let mut intervals_builder: ahash::AHashMap<
            Vec<u8>,
            Vec<rust_lapper::Interval<u64, annotation::Annotation>>,
        > = ahash::AHashMap::new();

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(input);

        for result in reader.byte_records() {
            let annotation = annotation::Annotation::from_byte_record(&result?)?;
            let seqname = annotation.get_seqname();
            let interval = annotation.get_interval();

            intervals_builder
                .entry(seqname.to_vec())
                .and_modify(
                    |tree: &mut Vec<rust_lapper::Interval<u64, annotation::Annotation>>| {
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
            rust_lapper::Lapper<u64, annotation::Annotation>,
        > = ahash::AHashMap::with_capacity(intervals_builder.len());
        for (key, values) in intervals_builder.drain() {
            intervals.insert(key, rust_lapper::Lapper::new(values));
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
            chr.find(interval.start, interval.end)
                .map(|e| &e.val)
                .collect()
        } else {
            vec![]
        }
    }

    /// Add annotation
    fn add_annotion(
        tree: &mut Vec<rust_lapper::Interval<u64, annotation::Annotation>>,
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

            tree.push(rust_lapper::Interval {
                start: upstream,
                stop: interval.start,
                val: annotation::Annotation::from_annotation(&annotation, b"upstream"),
            });
            tree.push(rust_lapper::Interval {
                start: interval.end,
                stop: interval.end + updown_distance,
                val: annotation::Annotation::from_annotation(&annotation, b"downstream"),
            });
        }

        tree.push(rust_lapper::Interval {
            start: interval.start,
            stop: interval.end,
            val: annotation,
        })
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    use crate::annotation::Annotation;

    /* project use */
    use super::*;

    const DATA: &[u8] = b"chr1\tknownGene\ttranscript\t11869\t14409\t.\t+\t.\tgene_id=gene1
chr1\tknownGene\texon\t11869\t12227\t.\t+\t.\tgene_id=gene1;transcript_id=gene1;exon_number=1;exon_id=gene1.1
chr1\tknownGene\texon\t12613\t12721\t.\t+\t.\tgene_id=gene1;transcript_id=gene1;exon_number=2;exon_id=gene1.2
chr1\tknownGene\texon\t13221\t14409\t.\t+\t.\tgene_id=gene1;transcript_id=gene1;exon_number=3;exon_id=gene1.3
chr1\tknownGene\ttranscript\t17369\t17436\t.\t-\t.\tgene_id=gene2;transcript_id=gene2
chr1\tknownGene\texon\t17369\t17436\t.\t-\t.\tgene_id=gene2;transcript_id=gene2;exon_number=1;exon_id=gene2.1
chr1\tknownGene\ttranscript\t29554\t31097\t.\t+\t.\tgene_id=gene3;transcript_id=gene3
chr1\tknownGene\texon\t29554\t30039\t.\t+\t.\tgene_id=gene3;transcript_id=gene3;exon_number=1;exon_id=gene3.1
chr1\tknownGene\texon\t30564\t30667\t.\t+\t.\tgene_id=gene3;transcript_id=gene3;exon_number=2;exon_id=gene3.2
chr1\tknownGene\texon\t30976\t31097\t.\t+\t.\tgene_id=gene3;transcript_id=gene3;exon_number=3;exon_id=gene3.3
chr2\ttest\ttranscript\t50\t200\t.\t-\t.\tgene_id=gene4";

    #[test]
    fn annotations() -> error::Result<()> {
        let annotations =
            AnnotationsDataBase::from_reader(std::io::BufReader::new(Box::new(DATA)), 100)?;

        let b1 = Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(DATA[..57].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;
        println!("{}", String::from_utf8(DATA[262..363].to_vec()).unwrap());
        let b2 = Annotation::from_byte_record(&csv::ByteRecord::from(
            String::from_utf8(DATA[262..363].to_vec())
                .unwrap()
                .split('\t')
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        ))?;

        let truth = vec![&b1, &b2];
        assert_eq!(annotations.get_annotation(b"chr1", 14000..14001), truth);

        // seqname not present
        assert_eq!(
            annotations.get_annotation(b"chrX", 14000..14001),
            Vec::<&annotation::Annotation>::new()
        );

        Ok(())
    }
}
