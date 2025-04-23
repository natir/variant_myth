//! A variant annotator.

#![warn(missing_docs)]

/* std use */

/* crate use */
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */

/* mod declaration */
pub mod annotation;
pub mod annotations_db;
pub mod effect;
pub mod error;
pub mod memoizor;
pub mod myth;
pub mod output;
pub mod positions;
pub mod sequences_db;
pub mod translate;
pub mod variant;
pub mod variant2myth;

#[cfg(feature = "cli")]
pub mod cli;

#[cfg(test)]
pub mod test_data;

#[cfg(feature = "json")]
fn serialize_bstr<T, S>(v: T, serializer: S) -> Result<S::Ok, S::Error>
where
    T: AsRef<[u8]>,
    S: serde::Serializer,
{
    serializer.serialize_str(unsafe { std::str::from_utf8_unchecked(v.as_ref()) })
}

/// For each variants found matching annotations
#[cfg(not(feature = "parallel"))]
pub fn vcf2myth<R>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    vcf_reader: variant::VcfReader<R>,
    annotators_choices: variant2myth::AnnotatorsChoices,
    mut writer: Box<dyn output::MythWriter>,
) -> error::Result<()>
where
    R: std::io::BufRead,
{
    let variant2myth =
        variant2myth::Variant2Myth::new(annotations, translate, sequences, annotators_choices);

    for result in vcf_reader {
        let variant = result?;

        let myth = variant2myth.myth(variant);

        writer.write_myth(myth)?;
    }

    writer.close()?;

    Ok(())
}

/// For each variants found matching annotations
#[cfg(feature = "parallel")]
pub fn vcf2myth<R>(
    annotations: &annotations_db::AnnotationsDataBase,
    sequences: &sequences_db::SequencesDataBase,
    translate: &translate::Translate,
    vcf_reader: variant::VcfReader<R>,
    annotators_choices: variant2myth::AnnotatorsChoices,
    mut writer: Box<dyn output::MythWriter + std::marker::Send>,
) -> error::Result<()>
where
    R: std::io::BufRead + std::marker::Send,
{
    let (tx, rx) = std::sync::mpsc::channel::<myth::Myth>();

    let write_thread = std::thread::spawn(move || -> error::Result<()> {
        for myth in rx {
            writer.write_myth(myth)?;
        }
        writer.close()?;
        Ok(())
    });

    let variant2myth =
        variant2myth::Variant2Myth::new(annotations, translate, sequences, annotators_choices);

    let results = vcf_reader
        .par_bridge()
        .filter(Result::is_ok)
        .map(error::Result::unwrap)
        .map(|variant| tx.send(variant2myth.myth(variant)))
        .filter(|r| r.is_err())
        .collect::<Vec<core::result::Result<(), std::sync::mpsc::SendError<myth::Myth>>>>();

    for result in results {
        result?
    }

    drop(tx);

    write_thread.join().unwrap() // Err only if panic! if write thread panic all should panic !
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crate use */

    /* project use */
}
