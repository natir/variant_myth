//! Struct to convert exon an reference sequence in different type of sequence

/* std use */

/* crate use */

/* project use */
use crate::annotation;
use crate::sequences_db;
use crate::translate;
use crate::variant;

struct Exons2Sequences<'a> {
    sequences: &'a sequences_db::SequencesDataBase,
    translate: &'a translate::Translate,
}

impl<'a> Exons2Sequences<'a> {
    /// Create a new Exons2Sequences from sequences database and translate table
    pub fn new(
        sequences: &'a sequences_db::SequencesDataBase,
        translate: &'a translate::Translate,
    ) -> Self {
        Self {
            sequences,
            translate,
        }
    }

    /// Get
    pub fn epissed(&self, exon_annot: &[&annotation::Annotation]) -> Vec<u8> {
        let intervals: Vec<(core::ops::Range<u64>, annotation::Strand)> = exon_annot
            .iter()
            .map(|obj| (obj.get_start()..obj.get_stop(), *obj.get_strand()))
            .collect();

        self.sequences
            .get_transcript(exon_annot[0].get_seqname(), &intervals)
    }

    pub fn epissed_with_variant(
        &self,
        exon_target: usize,
        exon_annot: &[&annotation::Annotation],
        variant: &variant::Variant,
    ) -> Vec<u8> {
        let epissed = self.epissed(exon_annot);

        let var_pos_in_trans = exon_annot[..exon_target]
            .iter()
            .map(|e| e.get_stop() - e.get_start())
            .sum::<u64>()
            + variant.position
            - exon_annot[exon_target].get_start();

        epissed[..=var_pos_in_trans as usize]
            .iter()
            .chain(&variant.alt_seq)
            .chain(&epissed[(var_pos_in_trans as usize + variant.ref_seq.len())..])
            .cloned()
            .collect::<Vec<u8>>()
    }

    pub fn coding(
        &self,
        exon_annot: &[&annotation::Annotation],
        start_position: std::option::Option<u64>,
        stop_position: std::option::Option<u64>,
    ) -> Vec<u8> {
        vec![]
    }

    pub fn coding_with_variant(
        &self,
        exon_annot: &[&annotation::Annotation],
        start_position: std::option::Option<u64>,
        stop_position: std::option::Option<u64>,
        variant: &variant::Variant,
    ) -> Vec<u8> {
        vec![]
    }

    pub fn aa(
        &self,
        exon_annot: &[&annotation::Annotation],
        start_position: std::option::Option<u64>,
        stop_position: std::option::Option<u64>,
    ) -> Vec<u8> {
        vec![]
    }

    pub fn aa_with_variant(
        &self,
        exon_annot: &[&annotation::Annotation],
        start_position: std::option::Option<u64>,
        stop_position: std::option::Option<u64>,
        variant: &variant::Variant,
    ) -> Vec<u8> {
        vec![]
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    use annotation::Annotation;
    /* crate use */
    use biotest::Format as _;

    /* project use */
    use super::*;
    use crate::error;

    #[test]
    fn epissed() -> error::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(100_000).build()?;

        let mut temp_input = vec![];
        generator.record(&mut temp_input, &mut rng)?;
        let input: std::io::BufReader<Box<dyn std::io::Read + Send>> =
            std::io::BufReader::new(Box::new(std::io::Cursor::new(temp_input.to_vec())));

        let sequences = sequences_db::SequencesDataBase::from_reader(input)?;
        let translate = translate::Translate::default();

        let exons2seq = Exons2Sequences::new(&sequences, &translate);

        let annotations = (0..10)
            .map(|i| {
                Annotation::from_byte_record(&csv::ByteRecord::from(vec![
                    "GSWNPZYBHL",
                    "knownGene",
                    "transcript",
                    &(i * 100 + 34).to_string(),
                    &(i * 100 + 134).to_string(),
                    ".",
                    "+",
                    ".",
                    "",
                ]))
            })
            .collect::<error::Result<Vec<annotation::Annotation>>>()?;

        let sub_annot = annotations[2..8]
            .iter()
            .collect::<Vec<&annotation::Annotation>>();

        assert_eq!(exons2seq.epissed(&sub_annot), b"GaGCatAGGACAAaacTaTTagagGtatAGCcTatTtaaaaCGgcttGGTtgaCtgACTacgtCTaTgTCAGgCtaGTtcCCTcgcTgAgGgAtCAAatTCTATTGTaggcGCaCcCGtCtATgTTgTATcaTTCGaCCttcAaGCGCAatgaTGAtaatcaCtGcTAGCCAgaTTgcAaTtaTGgACTTagGgtATACCtcTctCAtgCGCagTCTcaacCATAtGtGgtAtacAagtTGgAtgcGtTCtctTgctTtcGggATtcGAgtaTgacgtCCTAtActaGAggcAAGGACGaATctgCaaatgctgTcCaAgttcGtGAtcAttaTtGgCACgCcgcCgATtcGCaTatTGGGCTacgtgACCGttTCAttTacAGCaTctttAAgAcCGgACTctgTGTtaAGCAgcagAcGttCagTgCTAtccTGAAccCaaAcacagCATCTaTCgGcgcaGCaCaTATTacCGaTtgttCGTTaGccGaCAaGCGGATCgGGGATCaAaGcaACCGaTcGGCCGgGacTcATCTcaGcCgTGAgtTTACatTagaCtTTtccCCcAgagtctAGCCtCTgATTtTGCcGcGgCgTcG".to_vec());

        Ok(())
    }
}
