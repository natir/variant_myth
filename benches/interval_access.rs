use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::Rng;
use rand::SeedableRng;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

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
        if record.get(1) != Some(b"chr1") {
            continue;
        }

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
    let mut g = c.benchmark_group("request");

    g.sample_size(10);

    let mut requests = Vec::new();
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);

    for _ in 0..7_000_000 {
        let a = rng.gen_range(1..248_937_043);
        let b = rng.gen_range(1..5);
        requests.push(a..a + b);
    }

    let array_tree = create_tree();

    g.bench_function("array", |b| {
        b.iter(|| {
            for req in requests.iter() {
                black_box(array_tree.find(req.clone()));
            }
        })
    });

    #[cfg(feature = "parallel")]
    g.bench_function("array_rayon", |b| {
        b.iter(|| {
            black_box(
                requests
                    .par_iter()
                    .map(|req| {
                        array_tree.find(req.clone());
                    })
                    .collect::<Vec<_>>(),
            );
        })
    });
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
