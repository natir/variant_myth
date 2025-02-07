# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

Keep in mind as you contribute, that code, docs and other material submitted to this projects are considered licensed under  license.

## Setup developement environment

We recommand to install rust with [rustup](https://rustup.rs/).
If you want perform a documentation contribution install [mdbook](https://rust-lang.github.io/mdBook/guide/installation.html).

## Contribution

Before start any modification please create a specific branch:
```bash
git switch -c fix_11         # branch create to fix issue 11
git switch -c feat_index_rc  # branch to add a new index reverse complement method
```

### Code contribution

Before submit pull request make sure you run:

```bash
cargo fmt
cargo clippy
cargo test
```

You can check your new code are covered by run:
```bash
cargo tarpaulin
```
And open `target/coverage/tarpaulin-report.html`

### Dataset acquisition and generation

Get human data:
```bash
./get_data.py download
```
This call create a `data` directory with file `annotations.gff3.gz`, `references.fasta.gz` and `variants.vcf.gz`, this file are usefull "real use case" run.

You could generate a subsample of `data` directory content for a manual check with command:
```bash
./get_data.py subsample -a data/annotations.gff3.gz -r data/references.fasta.gz -A data/sub_annotations.gff3 -R data/sub_references.fasta -V data/sub_variants.vcf --seed 1234
```

File in `tests/data/`, used for functional test, are generate by command:
```bash
mkdir -p tests/data/
./get_data.py subsample -a data/annotations.gff3.gz -r data/references.fasta.gz -A tests/data/annotations.gff3 -R tests/data/references.fasta -V tests/data/variants.vcf --seed 42
```


### Documentation pull request

After change you can run:
```
cargo doc
```
And open `target/doc/variant_myth/index.html` to check effect of your change.

### Website pull request

After change you can run:
```
mdbook serve
```
To check effect of your change.
