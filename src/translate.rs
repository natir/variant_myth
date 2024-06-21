//! Struct to perform translate

/* std use */
use std::io::BufRead as _;

/* crate use */

/* project use */
use crate::error;

/// Convert a sequence in 2 bit representation if suseq is larger than 32 only the last 32 nuc is store
#[inline(always)]
pub fn seq2bit(subseq: &[u8]) -> u8 {
    let mut kmer: u8 = 0;

    for n in subseq {
        kmer <<= 2;
        kmer |= nuc2bit(*n);
    }

    kmer
}

/// Convert a nucleotide in 2bit representation
#[inline(always)]
pub fn nuc2bit(nuc: u8) -> u8 {
    (nuc >> 1) & 0b11
}

/// Perform a sequence translation
pub struct Translate {
    aa: [u8; 64],
    start: [bool; 64],
    end: [bool; 64],
}

impl Translate {
    /// Build a Translate table from a reader
    pub fn from_reader(
        mut input: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>,
    ) -> error::Result<Self> {
        let mut aa = [b'-'; 64];
        let mut start = [false; 64];
        let mut end = [false; 64];

        let mut buffer = Vec::new();
        input.read_until(b'\n', &mut buffer)?;
        let peptide = buffer[11..75].to_ascii_uppercase();
        buffer.clear();

        input.read_until(b'\n', &mut buffer)?;
        let start_stop = buffer[11..75].to_ascii_uppercase();
        buffer.clear();

        input.read_until(b'\n', &mut buffer)?;
        let base1 = buffer[11..75].to_ascii_uppercase();
        buffer.clear();

        input.read_until(b'\n', &mut buffer)?;
        let base2 = buffer[11..75].to_ascii_uppercase();
        buffer.clear();

        input.read_until(b'\n', &mut buffer)?;
        let base3 = buffer[11..75].to_ascii_uppercase();
        buffer.clear();

        for (p, (s, (b1, (b2, b3)))) in peptide.into_iter().zip(
            start_stop.into_iter().zip(
                base1
                    .into_iter()
                    .zip(base2.into_iter().zip(base3.into_iter())),
            ),
        ) {
            let codon_hash = seq2bit(&[b1, b2, b3]);

            aa[codon_hash as usize] = p;
            start[codon_hash as usize] = s == b'M';
            end[codon_hash as usize] = s == b'*';
        }

        Ok(Translate { aa, start, end })
    }

    /// Get acide amine associate with codon
    pub fn get_aa(&self, codon: &[u8]) -> u8 {
        self.aa[seq2bit(codon) as usize]
    }

    /// Codon is start
    pub fn is_start(&self, codon: &[u8]) -> bool {
        self.start[seq2bit(codon) as usize]
    }

    /// Codon is stop
    pub fn is_stop(&self, codon: &[u8]) -> bool {
        self.end[seq2bit(codon) as usize]
    }

    /// Translate sequence
    pub fn translate(&self, seq: &[u8]) -> Vec<u8> {
        seq.to_ascii_uppercase()
            .chunks_exact(3)
            .map(|codon| self.get_aa(codon))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */
    use biotest::Format as _;

    /* project use */
    use super::*;

    const STANDARD: &[u8] =
        b"    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M------**--*----M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
";

    #[test]
    fn check_table() -> error::Result<()> {
        let trans = Translate::from_reader(std::io::BufReader::new(Box::new(STANDARD)))?;

        assert_eq!(trans.get_aa(&[b'T', b'T', b'T']), b'F');
        assert_eq!(trans.get_aa(&[b'A', b'T', b'G']), b'M');
        assert_eq!(trans.get_aa(&[b'G', b'G', b'G']), b'G');

        assert!(trans.is_start(&[b'T', b'T', b'G']));
        assert!(trans.is_start(&[b'C', b'T', b'G']));
        assert!(trans.is_start(&[b'A', b'T', b'G']));
        assert!(!trans.is_start(&[b'A', b'T', b'C']));

        assert!(trans.is_stop(&[b'T', b'A', b'A']));
        assert!(trans.is_stop(&[b'T', b'A', b'G']));
        assert!(trans.is_stop(&[b'T', b'G', b'A']));
        assert!(!trans.is_stop(&[b'A', b'T', b'C']));

        Ok(())
    }

    #[test]
    fn translate() -> error::Result<()> {
        let trans = Translate::from_reader(std::io::BufReader::new(Box::new(STANDARD)))?;

        let mut rng = biotest::rand();
        let generate = biotest::Sequence::default();
        let mut seq = Vec::new();
        generate.record(&mut seq, &mut rng)?;

        assert_eq!(
            trans.translate(&seq),
            b"YMNRVLVKPR*CLYAGYRIIDGCSCLLVLCKRGDMLQLPLTGIHPLELAT"
        );

        Ok(())
    }
}
