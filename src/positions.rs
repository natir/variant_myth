//! Define object related to position in different genomic point of view

/* std use */

/* crate use */
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */
use crate::annotation;

#[derive(Debug, Clone, PartialEq, Eq)]
/// Object that store a Genomic position
struct Genomic {
    position: u64,
}

impl Genomic {
    /// Create a new Genomic position
    pub fn new(position: u64) -> Self {
        Genomic { position }
    }

    /// Get the position
    pub fn position(&self) -> &u64 {
        &self.position
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// Object that store a Coding position
struct Coding {
    position: u64,
    not_coding: i64,
}

impl Coding {
    /// Create a new Coding position
    pub fn new(position: u64, not_coding: i64) -> Self {
        Coding {
            position,
            not_coding,
        }
    }

    /// Get position in coding
    pub fn position(&self) -> &u64 {
        &self.position
    }

    /// Get position in not coding
    ///
    /// if return value is 0 position are a coding position
    pub fn not_coding(&self) -> &i64 {
        &self.not_coding
    }

    /// Create a Coding position from:
    /// - Genomic value
    /// - exons limits
    /// - position of start codon
    /// - position of stop codon
    /// - strand of coding region
    ///
    fn try_from_genomic(
        genomic: Genomic,
        exons: Vec<u64>,
        start: u64,
        stop: u64,
        strand: annotation::Strand,
    ) -> Option<Self> {
        let position = if strand == annotation::Strand::Reverse {
            -(*genomic.position() as i64)
        } else {
            *genomic.position() as i64
        };

        let mut tmp: Vec<i64> = exons.iter().map(|p| *p as i64).collect();
        tmp.extend([start as i64, stop as i64]);

        if strand == annotation::Strand::Reverse {
            tmp.iter_mut().for_each(|v| *v *= -1);
        }

        #[cfg(not(feature = "parallel"))]
        tmp.sort();
        #[cfg(feature = "parallel")]
        tmp.par_sort();

        let not_coding_limit: Vec<std::ops::Range<i64>> =
            tmp.chunks(2).map(|x| (x[0])..(x[1])).collect();

        let mut tmp: Vec<i64> = exons[1..].iter().map(|p| *p as i64).collect();
        tmp.pop();
        tmp.extend([start as i64, stop as i64]);
        if strand == annotation::Strand::Reverse {
            tmp.iter_mut().for_each(|v| *v *= -1);
        }

        #[cfg(not(feature = "parallel"))]
        tmp.sort();
        #[cfg(feature = "parallel")]
        tmp.par_sort();
        let coding_limit: Vec<std::ops::Range<i64>> = tmp.chunks(2).map(|x| x[0]..x[1]).collect();

        dbg!(position);
        dbg!(&not_coding_limit);
        dbg!(&coding_limit);

        match not_coding_limit.iter().position(|x| x.contains(&position)) {
            Some(index) => {
                dbg!(index);

                if index == 0 {
                    Some(Coding::new(
                        1,
                        -(position.abs_diff(coding_limit[0].start) as i64),
                    ))
                } else {
                    Some(Coding::new(
                        (coding_limit[..index]
                            .iter()
                            .map(|x| x.end - x.start)
                            .sum::<i64>()) as u64,
                        position.abs_diff(coding_limit[index - 1].end) as i64,
                    ))
                }
            }
            None => coding_limit
                .iter()
                .position(|x| x.contains(&position))
                .map(|index| {
                    dbg!(&coding_limit);
                    dbg!(index);
                    dbg!(position);

                    Some(Coding::new(
                        (coding_limit[..index]
                            .iter()
                            .map(|x| x.end - x.start)
                            .sum::<i64>()
                            + position
                            - coding_limit[index].start) as u64,
                        0,
                    ))
                })?,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// Object that store a Transcript position
struct Transcript {
    position: i64,
}

impl Transcript {
    /// Create a new Transcript position
    pub fn new(position: i64) -> Self {
        Self { position }
    }

    /// Get position in transcript
    pub fn position(&self) -> &i64 {
        &self.position
    }

    /// Create a Transcript position from Coding value
    ///
    /// Return None if coding position isn't in transript
    pub fn from_coding(coding: Coding) -> Self {
        if coding.not_coding().is_negative() {
            Transcript::new(*coding.not_coding())
        } else {
            Transcript::new(*coding.position() as i64)
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// Object that store a Proteine position
struct Protein {
    position: u64,
}

impl Protein {
    /// Create a new Protein position
    pub fn new(position: u64) -> Self {
        Self { position }
    }

    /// Get position in protein
    pub fn position(&self) -> &u64 {
        &self.position
    }

    /// Create a Protein position from Transcript value
    ///
    /// Return None if transcript position is negative
    fn try_from_transcript(transcript: Transcript) -> Option<Self> {
        if transcript.position.is_negative() {
            None
        } else {
            Some(Protein::new(
                (transcript.position as u64)
                    .saturating_div(3)
                    .saturating_add(1),
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */
    use super::*;
    use crate::error;

    #[test]
    fn genomic() {
        let obj = Genomic::new(10);

        assert_eq!(obj.position(), &10);
    }

    #[test]
    fn coding_forward() {
        let obj = Genomic::new(10);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990],
            110,
            999,
            crate::annotation::Strand::Forward,
        );

        assert!(coding.is_none());

        let obj = Genomic::new(105);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990],
            110,
            980,
            crate::annotation::Strand::Forward,
        );

        assert!(coding.is_some());
        let tmp = coding.unwrap();
        assert_eq!(tmp.position(), &1);
        assert_eq!(tmp.not_coding(), &-5);

        let obj = Genomic::new(152);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990],
            110,
            980,
            crate::annotation::Strand::Forward,
        );

        assert!(coding.is_some());
        let tmp = coding.unwrap();
        assert_eq!(tmp.position(), &42);
        assert_eq!(tmp.not_coding(), &0);

        let obj = Genomic::new(995);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990, 1000, 1990],
            110,
            1980,
            crate::annotation::Strand::Forward,
        );

        assert!(coding.is_some());
        let tmp = coding.unwrap();
        assert_eq!(tmp.position(), &880);
        assert_eq!(tmp.not_coding(), &5);

        let obj = Genomic::new(1985);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990, 1000, 1990],
            110,
            1980,
            crate::annotation::Strand::Forward,
        );

        assert!(coding.is_some());
        let tmp = coding.unwrap();
        assert_eq!(tmp.position(), &1860);
        assert_eq!(tmp.not_coding(), &5);
    }

    #[test]
    fn coding_reverse() {
        let obj = Genomic::new(10);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990],
            110,
            980,
            crate::annotation::Strand::Reverse,
        );

        assert!(coding.is_none());
        let obj = Genomic::new(985);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990],
            110,
            980,
            crate::annotation::Strand::Reverse,
        );

        assert!(coding.is_some());
        let tmp = coding.unwrap();
        assert_eq!(tmp.position(), &1);
        assert_eq!(tmp.not_coding(), &-5);

        let obj = Genomic::new(152);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990],
            110,
            980,
            crate::annotation::Strand::Reverse,
        );

        assert!(coding.is_some());
        let tmp = coding.unwrap();
        assert_eq!(tmp.position(), &828);
        assert_eq!(tmp.not_coding(), &0);

        let obj = Genomic::new(995);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990, 1000, 1990],
            110,
            1980,
            crate::annotation::Strand::Reverse,
        );

        assert!(coding.is_some());
        let tmp = coding.unwrap();
        assert_eq!(tmp.position(), &980);
        assert_eq!(tmp.not_coding(), &5);

        let obj = Genomic::new(1985);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990, 1000, 1990],
            110,
            1980,
            crate::annotation::Strand::Reverse,
        );

        assert!(coding.is_some());
        let tmp = coding.unwrap();
        assert_eq!(tmp.position(), &1);
        assert_eq!(tmp.not_coding(), &-5);
    }

    #[test]
    fn transcript() -> error::Result<()> {
        let obj = Genomic::new(105);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990],
            110,
            980,
            crate::annotation::Strand::Forward,
        )
        .ok_or(anyhow::anyhow!(
            "positions::Genomic 2 positions::Coding failled"
        ))?;

        let transcript = Transcript::from_coding(coding);
        assert_eq!(transcript.position(), &-5);

        let obj = Genomic::new(152);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990],
            110,
            980,
            crate::annotation::Strand::Forward,
        )
        .ok_or(anyhow::anyhow!(
            "positions::Genomic 2 positions::Coding failled"
        ))?;

        let transcript = Transcript::from_coding(coding);
        assert_eq!(transcript.position(), &42);

        let obj = Genomic::new(995);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990, 1000, 1990],
            110,
            1980,
            crate::annotation::Strand::Forward,
        )
        .ok_or(anyhow::anyhow!(
            "positions::Genomic 2 positions::Coding failled"
        ))?;

        let transcript = Transcript::from_coding(coding);
        assert_eq!(transcript.position(), &880);

        let obj = Genomic::new(1985);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990, 1000, 1990],
            110,
            1980,
            crate::annotation::Strand::Forward,
        )
        .ok_or(anyhow::anyhow!(
            "positions::Genomic 2 positions::Coding failled"
        ))?;

        let transcript = Transcript::from_coding(coding);
        assert_eq!(transcript.position(), &1860);

        Ok(())
    }

    #[test]
    fn protein() -> error::Result<()> {
        let obj = Genomic::new(152);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990],
            110,
            980,
            crate::annotation::Strand::Forward,
        )
        .ok_or(anyhow::anyhow!(
            "positions::Genomic 2 positions::Coding failled"
        ))?;
        let transcript = Transcript::from_coding(coding);
        let protein = Protein::try_from_transcript(transcript);
        assert!(protein.is_some());
        let tmp = protein.unwrap();
        assert_eq!(tmp.position(), &15);

        let obj = Genomic::new(995);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990, 1000, 1990],
            110,
            1980,
            crate::annotation::Strand::Forward,
        )
        .ok_or(anyhow::anyhow!(
            "positions::Genomic 2 positions::Coding failled"
        ))?;
        let transcript = Transcript::from_coding(coding);
        let protein = Protein::try_from_transcript(transcript);
        assert!(protein.is_some());
        let tmp = protein.unwrap();
        assert_eq!(tmp.position(), &294);

        let obj = Genomic::new(1985);
        let coding = Coding::try_from_genomic(
            obj,
            vec![100, 990, 1000, 1990],
            110,
            1980,
            crate::annotation::Strand::Forward,
        )
        .ok_or(anyhow::anyhow!(
            "positions::Genomic 2 positions::Coding failled"
        ))?;
        let transcript = Transcript::from_coding(coding);
        let protein = Protein::try_from_transcript(transcript);
        assert!(protein.is_some());
        let tmp = protein.unwrap();
        assert_eq!(tmp.position(), &621);

        let protein = Protein::try_from_transcript(Transcript::new(-1));
        assert!(protein.is_none());

        Ok(())
    }
}
