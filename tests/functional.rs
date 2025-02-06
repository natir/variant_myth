/* std use */

/* crate use */
#[cfg(not(feature = "out_json"))]
use arrow::array::Array as _;

/* project use */

#[cfg(feature = "out_json")]
fn dive_in_tree(object: serde_json::Value, objects: &mut ahash::AHashSet<String>) -> () {
    match object {
        serde_json::Value::Null => (),
        serde_json::Value::Bool(b) => {
            objects.insert(format!("{}", b));
        }
        serde_json::Value::Number(n) => {
            objects.insert(format!("{}", n));
        }
        serde_json::Value::String(s) => {
            objects.insert(s);
        }
        serde_json::Value::Array(v) => {
            v.iter().for_each(|o| dive_in_tree(o.clone(), objects));
        }
        serde_json::Value::Object(m) => {
            m.iter().for_each(|(k, v)| {
                objects.insert(k.to_string());
                dive_in_tree(v.clone(), objects);
            });
        }
    }
}

#[cfg(feature = "out_json")]
fn compare_by_record<P, R>(truth_path: P, result_path: R) -> anyhow::Result<()>
where
    P: std::convert::AsRef<std::path::Path>,
    R: std::convert::AsRef<std::path::Path>,
{
    let mut truth = ahash::AHashSet::new();
    serde_json::Deserializer::from_reader(
        std::fs::File::open(&truth_path).map(std::io::BufReader::new)?,
    )
    .into_iter()
    .map(std::result::Result::unwrap)
    .for_each(|o| dive_in_tree(o, &mut truth));

    let mut result = ahash::AHashSet::new();
    serde_json::Deserializer::from_reader(
        std::fs::File::open(&result_path).map(std::io::BufReader::new)?,
    )
    .into_iter()
    .map(std::result::Result::unwrap)
    .for_each(|o| dive_in_tree(o, &mut result));

    if truth != result {
        return Err(anyhow::anyhow!("truth: {:?}\nresult: {:?}", truth, result));
    }

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
                if t != r {
                    return Err(anyhow::anyhow!(
                        "Column {}\n\ttruth: {:?}\n\tresult: {:?}",
                        column.name(),
                        t,
                        r
                    ));
                }
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
                if t != r {
                    return Err(anyhow::anyhow!(
                        "Column {}\n\ttruth: {:?}\n\tresult: {:?}",
                        column.name(),
                        t,
                        r
                    ));
                }
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

                if t != r {
                    return Err(anyhow::anyhow!(
                        "Column {}\n\ttruth: {:?}\n\tresult: {:?}",
                        column.name(),
                        t,
                        r
                    ));
                }
            }
            a => unreachable!("This column type isn't variant_myth schema {:?}", a),
        }
    }

    Ok(())
}

#[test]
fn annotator_gene() -> anyhow::Result<()> {
    let tmp_path = tempfile::tempdir()?.into_path();
    let mut truth_path = std::path::PathBuf::from("tests/data/truth/annotator_gene");
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

    cmd.args(&args);

    let assert = cmd.assert();

    if let Err(e) = assert.try_success() {
        eprintln!(
            "failled\n\targument: {}\n\toutput path {}",
            args.join(" "),
            output_path.to_str().unwrap()
        );
        return Err(e.into());
    }

    if cfg!(feature = "out_json") {
        truth_path.set_extension("json");
    } else {
        truth_path.set_extension("parquet");
    }

    if let Err(e) = compare_by_record(truth_path, &output_path) {
        eprintln!(
            "failled\n\targument: {}\n\toutput path {}",
            args.join(" "),
            output_path.to_str().unwrap()
        );
        Err(e)
    } else {
        Ok(())
    }
}

#[test]
fn annotator_feature() -> anyhow::Result<()> {
    let tmp_path = tempfile::tempdir()?.into_path();
    let mut truth_path = std::path::PathBuf::from("tests/data/truth/annotator_feature");
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
        "feature",
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

    cmd.args(&args);

    let assert = cmd.assert();

    if let Err(e) = assert.try_success() {
        eprintln!(
            "failled\n\targument: {}\n\toutput path {}",
            args.join(" "),
            output_path.to_str().unwrap()
        );
        return Err(e.into());
    }

    if cfg!(feature = "out_json") {
        truth_path.set_extension("json");
    } else {
        truth_path.set_extension("parquet");
    }

    if let Err(e) = compare_by_record(truth_path, &output_path) {
        eprintln!(
            "failled\n\targument: {}\n\toutput path {}",
            args.join(" "),
            output_path.to_str().unwrap()
        );
        Err(e)
    } else {
        Ok(())
    }
}

#[cfg(not(feature = "parallel"))]
#[test]
fn logging_updown_setup() -> anyhow::Result<()> {
    let tmp_path = tempfile::tempdir()?.into_path();
    let output_path = tmp_path.join("myth.parquet");

    let mut cmd = assert_cmd::Command::cargo_bin("variant_myth")?;
    let args = vec![
        "-vv",
        "-i",
        "tests/data/variants.vcf",
        "-r",
        "tests/data/references.fasta",
        "-a",
        "tests/data/annotations.gff3",
        "parquet",
        "-p",
        output_path.to_str().unwrap(),
    ];

    cmd.args(&args);

    let assert = cmd.assert();

    match assert.try_success() {
        Err(e) => {
            eprintln!(
                "failled\n\targument: {}\n\toutput path {}",
                args.join(" "),
                output_path.to_str().unwrap()
            );
            return Err(e.into());
        }
        Ok(o) => {
            o.stdout(&b""[..]).stderr(
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
        }
    }

    Ok(())
}
