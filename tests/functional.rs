/* std use */
use std::io::Write;

/* crate use */
use biotest::Format as _;

/* project use */

fn prepare_input() -> anyhow::Result<(
    std::path::PathBuf,
    std::path::PathBuf,
    std::path::PathBuf,
    std::path::PathBuf,
    std::path::PathBuf,
)> {
    let tmp_path = tempfile::tempdir()?.into_path();

    let mut rng = biotest::rand();

    // Generate sequence
    let mut sequence = vec![];
    biotest::Fasta::builder()
        .sequence_len(10_000)
        .comment_len(0)
        .build()?
        .records(&mut sequence, &mut rng, 5)?;
    let seq_path = tmp_path.join("sequence.fasta");
    std::fs::File::create(&seq_path)?.write_all(&sequence)?;

    let seqnames = (0..5)
        .map(|x| 1 + (10 + 1 + 10_000 + 3) * x)
        .map(|x| x..=x + 9)
        .map(|x| sequence[x].to_vec())
        .collect::<Vec<Vec<u8>>>();

    // Generate vcf
    let variant_path = tmp_path.join("variant.vcf");
    biotest::Vcf::builder()
        .header(
            biotest::format::vcf::header::Header::builder()
                .contigs(biotest::values::Chromosomes::UserDefine(vec![
                    b"GSWNPZYBHL", // content of seqnames
                    b"RTKBXIEJTM",
                    b"JRBUEQABBR",
                    b"VOYRLHIJLC",
                    b"UPFWPOWSBP",
                ]))
                .contig_length(10_000)
                .build()?,
        )
        .record(
            biotest::format::vcf::record::Record::builder()
                .contigs(biotest::values::Chromosomes::UserDefine(vec![
                    b"GSWNPZYBHL", // content of seqnames
                    b"RTKBXIEJTM",
                    b"JRBUEQABBR",
                    b"VOYRLHIJLC",
                    b"UPFWPOWSBP",
                ]))
                .position(biotest::values::Integer::UserDefine(0..10_000))
                .build()?,
        )
        .build()?
        .create(&variant_path, &mut rng, 100)?;

    // Generate translate table
    let translate_path = tmp_path.join("translate.txt");
    std::fs::File::create(&translate_path)?.write_all(variant_myth::translate::STANDARD)?;

    // Generate annotation
    let annotation_path = tmp_path.join("annotation.gff");
    let mut writer = std::fs::File::create(&annotation_path)?;
    for seq in seqnames {
        writeln!(
            writer,
            "{}\tvariant_myth\ttranscript\t1234\t4324\t.\t+\t.\ttranscript_id=transcript1",
            String::from_utf8(seq.to_vec()).unwrap()
        )?;
        writeln!(
            writer,
            "{}\tvariant_myth\tCDS\t1234\t4324\t.\t+\t.\ttranscript_id=transcript1",
            String::from_utf8(seq.to_vec()).unwrap()
        )?;
        writeln!(
            writer,
            "{}\tvariant_myth\t5UTR\t1234\t1346\t.\t+\t.\ttranscript_id=transcript1,exon_number=1",
            String::from_utf8(seq.to_vec()).unwrap()
        )?;
        writeln!(
            writer,
            "{}\tvariant_myth\texon\t1346\t1549\t.\t+\t.\ttranscript_id=transcript1,exon_number=1",
            String::from_utf8(seq.to_vec()).unwrap()
        )?;
        writeln!(
            writer,
            "{}\tvariant_myth\texon\t1623\t2624\t.\t+\t.\ttranscript_id=transcript1,exon_number=2",
            String::from_utf8(seq.to_vec()).unwrap()
        )?;
        writeln!(
            writer,
            "{}\tvariant_myth\texon\t2703\t3921\t.\t+\t.\ttranscript_id=transcript1,exon_number=3",
            String::from_utf8(seq.to_vec()).unwrap()
        )?;
        writeln!(
            writer,
            "{}\tvariant_myth\t3UTR\t3921\t4324\t.\t+\t.\ttranscript_id=transcript1,exon_number=1",
            String::from_utf8(seq.to_vec()).unwrap()
        )?;
    }

    // Output file
    let output_path = tmp_path.join("myth.json");

    Ok((
        variant_path,
        seq_path,
        annotation_path,
        translate_path,
        output_path,
    ))
}

#[test]
fn run() -> anyhow::Result<()> {
    let (variant_path, seq_path, annotation_path, translate_path, output_path) = prepare_input()?;

    let mut cmd = assert_cmd::Command::cargo_bin("variant_myth")?;
    cmd.args([
        "-i",
        &format!("{}", variant_path.display()),
        "-r",
        &format!("{}", seq_path.display()),
        "-a",
        &format!("{}", annotation_path.display()),
        "-t",
        &format!("{}", translate_path.display()),
        "-o",
        &format!("{}", output_path.display()),
    ]);

    let assert = cmd.assert();

    assert.success();

    std::fs::remove_dir_all(variant_path.parent().unwrap())?;

    Ok(())
}

#[cfg(feature = "parallel")]
#[test]
fn run_threads() -> anyhow::Result<()> {
    let (variant_path, seq_path, annotation_path, translate_path, output_path) = prepare_input()?;

    let mut cmd = assert_cmd::Command::cargo_bin("variant_myth")?;
    cmd.args([
        "-i",
        &format!("{}", variant_path.display()),
        "-r",
        &format!("{}", seq_path.display()),
        "-a",
        &format!("{}", annotation_path.display()),
        "-t",
        &format!("{}", translate_path.display()),
        "-o",
        &format!("{}", output_path.display()),
        "--threads",
        "4",
    ]);

    let assert = cmd.assert();

    assert.success();

    std::fs::remove_dir_all(variant_path.parent().unwrap())?;

    Ok(())
}

#[ignore]
#[cfg(not(feature = "parallel"))]
#[test]
fn logging_updown_setup() -> anyhow::Result<()> {
    let (variant_path, seq_path, annotation_path, translate_path, output_path) = prepare_input()?;

    let mut cmd = assert_cmd::Command::cargo_bin("variant_myth")?;
    cmd.args([
        "-vv",
        "-i",
        &format!("{}", variant_path.display()),
        "-r",
        &format!("{}", seq_path.display()),
        "-a",
        &format!("{}", annotation_path.display()),
        "-t",
        &format!("{}", translate_path.display()),
        "-o",
        &format!("{}", output_path.display()),
        "-d",
        "500",
    ]);

    let assert = cmd.assert();

    assert.success().stdout(&b""[..]).stderr(
        &b"INFO - Start read genome reference
INFO - End read genome reference
INFO - Start read annotations
INFO - End read annotations
INFO - Start read translation table
INFO - End read translation table
INFO - Start annotate variant
INFO - End annotate variant
"[..],
    );

    std::fs::remove_dir_all(variant_path.parent().unwrap())?;

    Ok(())
}
