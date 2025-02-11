//! An object to cache call of many database use to annotate variant

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::annotations_db;
use crate::error;
use crate::sequences_db;
use crate::variant;

/// Object cache call of Annotations and Sequences DataBase
pub struct Memoizor<'a> {
    transcript_id: &'a [u8],
    annotations_db: &'a annotations_db::AnnotationsDataBase,
    sequences: &'a sequences_db::SequencesDataBase,
    not_coding_annotations: &'a [&'a annotation::Annotation],

    option_coding_annotation: Option<Option<&'a Vec<annotation::Annotation>>>,
    option_exons_annotation: Option<Vec<annotation::Annotation>>,
    option_transcript: Option<Option<&'a annotation::Annotation>>,
    option_epissed: Option<Vec<u8>>,
    option_epissed_edit: Option<Vec<u8>>,
    option_coding: Option<Vec<u8>>,
    option_coding_edit: Option<Vec<u8>>,
}

impl<'a> Memoizor<'a> {
    /// Create a new Memoizor
    pub fn new(
        transcript_id: &'a [u8],
        annotations_db: &'a annotations_db::AnnotationsDataBase,
        sequences: &'a sequences_db::SequencesDataBase,
        not_coding_annotations: &'a [&'a annotation::Annotation],
    ) -> Self {
        Self {
            transcript_id,
            annotations_db,
            sequences,
            not_coding_annotations,
            option_coding_annotation: None,
            option_exons_annotation: None,
            option_transcript: None,
            option_epissed: None,
            option_epissed_edit: None,
            option_coding: None,
            option_coding_edit: None,
        }
    }

    /// Get not coding annotation associate with variant
    #[inline(always)]
    pub fn not_coding_annotation(&self) -> &'a [&'a annotation::Annotation] {
        self.not_coding_annotations
    }

    /// Get coding annotation present in transcript
    #[inline(always)]
    pub fn coding_annotation(&mut self) -> Option<&Vec<annotation::Annotation>> {
        if self.option_coding_annotation.is_none() {
            self.option_coding_annotation = Some(
                self.annotations_db
                    .get_coding_annotation(self.transcript_id),
            );
        }

        self.option_coding_annotation.unwrap() // value isn't none we check it
    }

    /// Get exon extract from annotation present in transcript
    #[inline(always)]
    pub fn exons_annotation(&mut self) -> &[annotation::Annotation] {
        if self.option_exons_annotation.is_none() {
            self.option_exons_annotation = Some(
                self.coding_annotation()
                    .unwrap()
                    .iter()
                    .filter(|a| a.get_feature() == b"exon")
                    .cloned()
                    .collect::<Vec<annotation::Annotation>>(),
            );
        }

        self.option_exons_annotation.as_ref().unwrap()
    }

    /// Get transcript annotation from id
    #[inline(always)]
    pub fn transcript(&mut self) -> Option<&annotation::Annotation> {
        if self.option_transcript.is_none() {
            self.option_transcript = Some(self.annotations_db.get_transcript(self.transcript_id));
        }

        self.option_transcript.unwrap() // value isn't none we check it
    }

    /// Get concatenation of sequence covered by annotations
    pub fn epissed(
        &mut self,
        annotations: &[&annotation::Annotation],
        strand: annotation::Strand,
    ) -> error::Result<&[u8]> {
        if self.option_epissed.is_none() {
            self.option_epissed = Some(self.sequences.epissed(annotations, strand)?);
        }

        Ok(self.option_epissed.as_ref().unwrap())
    }

    /// Get concatenation of sequence covered by annotations edited by variant
    pub fn epissed_edit(
        &mut self,
        annotations: &[&annotation::Annotation],
        strand: annotation::Strand,
        variant: &variant::Variant,
    ) -> error::Result<&[u8]> {
        if self.option_transcript.is_none() {
            self.option_epissed_edit =
                Some(self.sequences.epissed_edit(annotations, strand, variant)?);
        }

        Ok(self.option_epissed.as_ref().unwrap())
    }

    /// Get coding sequence covered by annotations
    pub fn coding(
        &mut self,
        annotations: &[&annotation::Annotation],
        strand: annotation::Strand,
        start_position: std::option::Option<u64>,
        stop_position: std::option::Option<u64>,
    ) -> error::Result<&[u8]> {
        if self.option_coding.is_none() {
            self.option_coding =
                Some(
                    self.sequences
                        .coding(annotations, strand, start_position, stop_position)?,
                );
        }

        Ok(self.option_coding.as_ref().unwrap())
    }

    /// Get coding sequence covered by annotations edited by variant
    pub fn coding_edit(
        &mut self,
        annotations: &[&annotation::Annotation],
        strand: annotation::Strand,
        variant: &variant::Variant,
        start_position: std::option::Option<u64>,
        stop_position: std::option::Option<u64>,
    ) -> error::Result<&[u8]> {
        if self.option_transcript.is_none() {
            self.option_coding_edit = Some(self.sequences.coding_edit(
                annotations,
                strand,
                variant,
                start_position,
                stop_position,
            )?);
        }

        Ok(self.option_coding.as_ref().unwrap())
    }
}
