//! Functional test check regression in help message

/* std use */

/* crate use */

/* project use */

const USAGE: &[u8] = b"A variant annotater

Usage: variant_myth [OPTIONS] --input <VARIANT_PATHS> --reference <REFERENCE_PATH> --annotations <ANNOTATIONS_PATH> <COMMAND>

Commands:
";

#[cfg(feature = "parquet")]
const SUBCOMMAND_PARQUET: &[u8] = b"  parquet  Output are write in parquet format
";
#[cfg(feature = "json")]
const SUBCOMMAND_JSON: &[u8] = b"  json     Output are write in json format
";
const SUBCOMMAND_HELP: &[u8] =
    b"  help     Print this message or the help of the given subcommand(s)

";

const LOCAL_OPTIONS: &[u8] = b"Options:
  -i, --input <VARIANT_PATHS>
          Variants path
  -r, --reference <REFERENCE_PATH>
          Reference genome path
  -a, --annotations <ANNOTATIONS_PATH>
          Annotation path
  -t, --translate <TRANSLATE_PATH>
          Translate table path, if not set use human
  -d, --updown-distance <UPDOWN_DISTANCE>
          [Up|Down]stream transcript distance, default: 5,000
  -c, --annotators-choices <ANNOTATORS_CHOICES>
          Select which type of annotation you want run [possible values: gene, feature, effect, hgvs]
";

#[cfg(feature = "parallel")]
const THREAD_OPTIONS: &[u8] = b"      --threads <THREADS>
          Number of theard use 0 use all avaible core, default value 0
";

const GENERIC_OPTIONS: &[u8] = b"  -q, --quiet
          Silence all output
  -v, --verbosity...
          Verbose mode (-v, -vv, -vvv, etc)
  -T, --timestamp <TS>
          Timestamp (sec, ms, ns, none)
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version
";

#[test]
fn help_message() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("variant_myth")?;
    cmd.args(["-h"]);

    let assert = cmd.assert();

    let mut help = Vec::new();
    help.extend(USAGE);

    #[cfg(feature = "parquet")]
    help.extend(SUBCOMMAND_PARQUET);
    #[cfg(feature = "json")]
    help.extend(SUBCOMMAND_JSON);

    help.extend(SUBCOMMAND_HELP);
    help.extend(LOCAL_OPTIONS);

    #[cfg(feature = "parallel")]
    help.extend(THREAD_OPTIONS);

    help.extend(GENERIC_OPTIONS);

    assert.success().stdout(help);

    Ok(())
}
