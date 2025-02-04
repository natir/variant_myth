/* std use */

/* crate use */
#[cfg(not(feature = "out_json"))]
use arrow::array::Array as _;

/* project use */

#[cfg(feature = "out_json")]
fn compare_by_record<P, R>(truth_path: P, result_path: R) -> anyhow::Result<()>
where
    P: std::convert::AsRef<std::path::Path>,
    R: std::convert::AsRef<std::path::Path>,
{
    let truth = serde_json::Deserializer::from_reader(
        std::fs::File::open(&truth_path).map(std::io::BufReader::new)?,
    )
    .into_iter::<serde_json::Value>()
    .map(std::result::Result::unwrap)
    .collect::<std::collections::HashSet<serde_json::Value>>();

    let result = serde_json::Deserializer::from_reader(
        std::fs::File::open(&result_path).map(std::io::BufReader::new)?,
    )
    .into_iter::<serde_json::Value>()
    .map(std::result::Result::unwrap)
    .collect::<std::collections::HashSet<serde_json::Value>>();

    assert_eq!(truth, result);

    Ok(())
}

#[cfg(not(feature = "out_json"))]
fn compare_by_record<P, R>(truth_path: P, result_path: R) -> anyhow::Result<()>
where
    P: std::convert::AsRef<std::path::Path>,
    R: std::convert::AsRef<std::path::Path>,
{
    let mut truth_reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(
        std::fs::File::open(&truth_path)?,
    )?
    .build()?;

    let truth = truth_reader.next().unwrap()?;

    let mut result_reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(
        std::fs::File::open(&result_path)?,
    )?
    .build()?;

    let result = result_reader.next().unwrap()?;

    assert_eq!(truth.schema(), result.schema());
    for column in truth.schema().fields() {
        match column.data_type() {
            arrow::datatypes::DataType::Utf8 => {
                // Convert an arrow StringArray in Rust Vec<String> is hard sorry
                let proxy = truth
                    .column_by_name(column.name())
                    .unwrap()
                    .as_any()
                    .downcast_ref::<arrow::array::StringArray>()
                    .unwrap();
                let mut t = Vec::with_capacity(proxy.len());
                for i in 0..proxy.len() {
                    t.push(proxy.value(i))
                }

                let proxy = result
                    .column_by_name(column.name())
                    .unwrap()
                    .as_any()
                    .downcast_ref::<arrow::array::StringArray>()
                    .unwrap();
                let mut r = Vec::with_capacity(proxy.len());
                for i in 0..proxy.len() {
                    r.push(proxy.value(i))
                }

                t.sort();
                r.sort();
                assert_eq!(t, r);
            }
            arrow::datatypes::DataType::UInt64 => {
                // Convert an arrow UInt64 in Rust Vec<UInt64> is hard sorry
                let proxy = truth
                    .column_by_name(column.name())
                    .unwrap()
                    .as_any()
                    .downcast_ref::<arrow::array::PrimitiveArray<arrow::datatypes::UInt64Type>>()
                    .unwrap();
                let mut t = Vec::with_capacity(proxy.len());
                for i in 0..proxy.len() {
                    t.push(proxy.value(i))
                }

                let proxy = result
                    .column_by_name(column.name())
                    .unwrap()
                    .as_any()
                    .downcast_ref::<arrow::array::PrimitiveArray<arrow::datatypes::UInt64Type>>()
                    .unwrap();
                let mut r = Vec::with_capacity(proxy.len());
                for i in 0..proxy.len() {
                    r.push(proxy.value(i))
                }

                t.sort();
                r.sort();
                assert_eq!(t, r);
            }
            arrow::datatypes::DataType::UInt8 => {
                // Convert an arrow UInt8 in Rust Vec<UInt8> is hard sorry
                let proxy = truth
                    .column_by_name(column.name())
                    .unwrap()
                    .as_any()
                    .downcast_ref::<arrow::array::PrimitiveArray<arrow::datatypes::UInt8Type>>()
                    .unwrap();
                let mut t = Vec::with_capacity(proxy.len());
                for i in 0..proxy.len() {
                    t.push(proxy.value(i))
                }

                let proxy = result
                    .column_by_name(column.name())
                    .unwrap()
                    .as_any()
                    .downcast_ref::<arrow::array::PrimitiveArray<arrow::datatypes::UInt8Type>>()
                    .unwrap();
                let mut r = Vec::with_capacity(proxy.len());
                for i in 0..proxy.len() {
                    r.push(proxy.value(i))
                }

                t.sort();
                r.sort();
                assert_eq!(t, r);
            }
            a => unreachable!("This column type isn't variant_myth schema {:?}", a),
        }
    }

    Ok(())
}

#[test]
fn annotator_gene() -> anyhow::Result<()> {
    let tmp_path = tempfile::tempdir()?.into_path();
    let mut cmd = assert_cmd::Command::cargo_bin("variant_myth")?;

    let mut args = Vec::new();

    args.extend([
        "-i",
        "tests/data/variants.vcf",
        "-r",
        "tests/data/references.fasta",
        "-a",
        "tests/data/annotations.gff3",
        "-c",
        "gene",
    ]);

    if cfg!(feature = "parallel") {
        args.extend(["--threads", "2"]);
    }

    let mut output_path = tmp_path.join("myth");
    if cfg!(feature = "out_json") {
        output_path.set_extension("json");
        args.extend(["json", "-p", output_path.to_str().unwrap(), "-f", "nd-json"]);
    } else {
        output_path.set_extension("parquet");
        args.extend(["parquet", "-p", output_path.to_str().unwrap()]);
    };

    cmd.args(args);

    let assert = cmd.assert();

    assert.success();

    dbg!(&output_path);
    if cfg!(feature = "out_json") {
        compare_by_record("tests/data/truth/annotator_gene.json", &output_path)
    } else {
        compare_by_record("tests/data/truth/annotator_gene.parquet", &output_path)
    }
}

#[cfg(not(feature = "parallel"))]
#[test]
fn logging_updown_setup() -> anyhow::Result<()> {
    let tmp_path = tempfile::tempdir()?.into_path();
    let output_path = tmp_path.join("myth.parquet");

    let mut cmd = assert_cmd::Command::cargo_bin("variant_myth")?;
    cmd.args([
        "-vv",
        "-i",
        "tests/data/variants.vcf",
        "-r",
        "tests/data/references.fasta",
        "-a",
        "tests/data/annotations.gff3",
        "parquet",
        "-p",
        &format!("{}", output_path.display()),
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

    Ok(())
}
