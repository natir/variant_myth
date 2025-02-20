//! Benchmark gene annotation

/* std use */

/* crate use */

/* project use */

fn gene_benchmark(c: &mut criterion::Criterion) {
    let tmp_path = tempfile::tempdir().unwrap().into_path();
    let mut cmd = assert_cmd::Command::cargo_bin("variant_myth").unwrap();

    let mut args = Vec::new();

    args.extend([
        "-i",
        "data/variants.vcf.gz",
        "-r",
        "data/references.fasta.gz",
        "-a",
        "data/annotations.gff3.gz",
        "-c",
        "gene",
    ]);

    if cfg!(feature = "parallel") {
        args.extend(["--threads", "4"]);
    }

    let mut output_path = tmp_path.join("myth");
    if cfg!(feature = "json") {
        output_path.set_extension("json");
        args.extend(["json", "-p", output_path.to_str().unwrap(), "-f", "nd-json"]);
    } else {
        output_path.set_extension("parquet");
        args.extend(["parquet", "-p", output_path.to_str().unwrap()]);
    };

    cmd.args(&args);

    c.bench_function("gene", |b| b.iter(|| std::hint::black_box(cmd.assert())));
}

criterion::criterion_group!(
    name = benches;
    config = criterion::Criterion::default().sample_size(10);
    targets = gene_benchmark
);
criterion::criterion_main!(benches);
