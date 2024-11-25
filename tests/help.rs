//! Functional test check regression in help message

/* std use */

/* crate use */

/* project use */

#[cfg(not(feature = "parallel"))]
const HELP: &[u8] = b"A variant annotater

Usage: variant_myth [OPTIONS] --output <OUTPUT_PATH> <COMMAND>

Commands:
  var2gene  Annotate variant with only gene name
  var2full  Annotate variants with all annotation
  help      Print this message or the help of the given subcommand(s)

Options:
  -o, --output <OUTPUT_PATH>               Output path
  -d, --updown-distance <UPDOWN_DISTANCE>  [Up|Down]stream transcript distance, default: 5,000
  -q, --quiet                              Silence all output
  -v, --verbosity...                       Verbose mode (-v, -vv, -vvv, etc)
  -T, --timestamp <TS>                     Timestamp (sec, ms, ns, none)
  -h, --help                               Print help
  -V, --version                            Print version
";

#[cfg(feature = "parallel")]
const HELP: &[u8] = b"A variant annotater

Usage: variant_myth [OPTIONS] --output <OUTPUT_PATH> <COMMAND>

Commands:
  var2gene  Annotate variant with only gene name
  var2full  Annotate variants with all annotation
  help      Print this message or the help of the given subcommand(s)

Options:
  -o, --output <OUTPUT_PATH>
          Output path
  -d, --updown-distance <UPDOWN_DISTANCE>
          [Up|Down]stream transcript distance, default: 5,000
      --threads <THREADS>
          Number of theard use 0 use all avaible core, default value 0
  -q, --quiet
          Silence all output
  -v, --verbosity...
          Verbose mode (-v, -vv, -vvv, etc)
  -T, --timestamp <TS>
          Timestamp (sec, ms, ns, none)
  -h, --help
          Print help
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
