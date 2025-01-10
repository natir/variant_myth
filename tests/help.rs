//! Functional test check regression in help message

/* std use */

/* crate use */

/* project use */

#[cfg(not(feature = "parallel"))]
const HELP: &[u8] = b"A variant annotater

Usage: variant_myth [OPTIONS] --input <VARIANT_PATH> --reference <REFERENCE_PATH> --annotations <ANNOTATIONS_PATH> --output <OUTPUT_PATH>

Options:
  -i, --input <VARIANT_PATH>
          Variant path
  -r, --reference <REFERENCE_PATH>
          Reference genome path
  -a, --annotations <ANNOTATIONS_PATH>
          Annotation path
  -t, --translate <TRANSLATE_PATH>
          Translate table path, if not set use human
  -o, --output <OUTPUT_PATH>
          Output path
  -d, --updown-distance <UPDOWN_DISTANCE>
          [Up|Down]stream transcript distance, default: 5,000
  -c, --annotators-choices <ANNOTATORS_CHOICES>
          Select which type of annotation you want run [possible values: gene, feature, effect, hgvs]
  -q, --quiet
          Silence all output
  -v, --verbosity...
          Verbose mode (-v, -vv, -vvv, etc)
  -T, --timestamp <TS>
          Timestamp (sec, ms, ns, none)
  -h, --help
          Print help (see more with \'--help\')
  -V, --version
          Print version
";

#[cfg(feature = "parallel")]
const HELP: &[u8] = b"A variant annotater

Usage: variant_myth [OPTIONS] --input <VARIANT_PATH> --reference <REFERENCE_PATH> --annotations <ANNOTATIONS_PATH> --output <OUTPUT_PATH>

Options:
  -i, --input <VARIANT_PATH>
          Variant path
  -r, --reference <REFERENCE_PATH>
          Reference genome path
  -a, --annotations <ANNOTATIONS_PATH>
          Annotation path
  -t, --translate <TRANSLATE_PATH>
          Translate table path, if not set use human
  -o, --output <OUTPUT_PATH>
          Output path
  -d, --updown-distance <UPDOWN_DISTANCE>
          [Up|Down]stream transcript distance, default: 5,000
  -c, --annotators-choices <ANNOTATORS_CHOICES>
          Select which type of annotation you want run [possible values: gene, feature, effect, hgvs]
      --threads <THREADS>
          Number of theard use 0 use all avaible core, default value 0
  -q, --quiet
          Silence all output
  -v, --verbosity...
          Verbose mode (-v, -vv, -vvv, etc)
  -T, --timestamp <TS>
          Timestamp (sec, ms, ns, none)
  -h, --help
          Print help (see more with \'--help\')
  -V, --version
          Print version
";

#[test]
fn help_message() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("variant_myth")?;
    cmd.args(["-h"]);

    let assert = cmd.assert();

    assert.success().stdout(HELP);

    Ok(())
}
