use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn create_tree() -> variant_myth::interval_tree::IntervalTree<u64, Vec<u8>> {
    let mut obj = variant_myth::interval_tree::IntervalTree::new();
    let file = std::io::BufReader::new(
        niffler::get_reader(Box::new(
            std::fs::File::open("hg38.knownGene.gff.gz").unwrap(),
        ))
        .unwrap()
        .0,
    );

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(file);

    for result in reader.byte_records() {
        let record = result.unwrap();
        let start = String::from_utf8(record.get(3).unwrap().to_vec())
            .unwrap()
            .parse::<u64>()
            .unwrap();
        let end = String::from_utf8(record.get(4).unwrap().to_vec())
            .unwrap()
            .parse::<u64>()
            .unwrap();
        obj.insert(start..end, record.get(8).unwrap().to_vec())
    }

    obj.index();

    obj
}

pub fn benchmark(c: &mut Criterion) {
    let mut g = c.benchmark_group("build");

    g.bench_function("local", |b| b.iter(|| black_box(create_tree())));
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
