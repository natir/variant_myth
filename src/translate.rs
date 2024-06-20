//! Struct to perform translate

/* std use */
use std::io::BufRead as _;

/* crate use */

/* project use */
use crate::error;

/// Perform a sequence translation
pub struct Translate {
    table: ahash::AHashMap<[u8; 3], u8>,
}

impl Translate {
    /// Build a Translate table from a reader
    pub fn from_reader(
        mut input: std::io::BufReader<Box<dyn std::io::Read + std::marker::Send>>,
    ) -> error::Result<Self> {
        let mut table = ahash::AHashMap::new();

        let mut buffer = Vec::new();
        input.read_until(b'\n', &mut buffer)?;
        let aa = buffer[10..].to_ascii_uppercase();
        buffer.clear();

        input.read_until(b'\n', &mut buffer)?;
        let _start_stop = buffer[10..].to_ascii_uppercase();
        buffer.clear();

        input.read_until(b'\n', &mut buffer)?;
        let base1 = buffer[10..].to_ascii_uppercase();
        buffer.clear();

        input.read_until(b'\n', &mut buffer)?;
        let base2 = buffer[10..].to_ascii_uppercase();
        buffer.clear();

        input.read_until(b'\n', &mut buffer)?;
        let base3 = buffer[10..].to_ascii_uppercase();
        buffer.clear();

        for (aa, (b1, (b2, b3))) in aa.into_iter().zip(
            base1
                .into_iter()
                .zip(base2.into_iter().zip(base3.into_iter())),
        ) {
            table.insert([b1, b2, b3], aa);
        }

        Ok(Translate { table })
    }

    /// Get table
    pub fn table(&self) -> &ahash::AHashMap<[u8; 3], u8> {
        &self.table
    }

    /// Translate sequence
    pub fn translate(&self, seq: &[u8]) -> Vec<u8> {
        seq.to_ascii_uppercase()
            .chunks_exact(3)
            .map(|codon| self.table.get(codon))
            .map(core::option::Option::unwrap)
            .cloned()
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

        assert_eq!(trans.table().get(&[b'T', b'T', b'T']), Some(&b'F'));
        assert_eq!(trans.table().get(&[b'A', b'T', b'G']), Some(&b'M'));
        assert_eq!(trans.table().get(&[b'G', b'G', b'G']), Some(&b'G'));

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
