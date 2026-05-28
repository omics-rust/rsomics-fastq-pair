use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;
use tempfile;

fn bench_fastq_pair(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-fastq-pair");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let r1 = manifest.join("tests/golden/r1.fq");
    let r2 = manifest.join("tests/golden/r2.fq");
    let dir = tempfile::tempdir().unwrap();
    let out1 = dir.path().join("out1.fq");
    let out2 = dir.path().join("out2.fq");
    c.bench_function("rsomics-fastq-pair golden", |b| {
        b.iter(|| {
            let status = Command::new(black_box(bin))
                .args([
                    r1.to_str().unwrap(),
                    r2.to_str().unwrap(),
                    "--out1",
                    out1.to_str().unwrap(),
                    "--out2",
                    out2.to_str().unwrap(),
                ])
                .status()
                .unwrap();
            assert!(status.success());
        });
    });
}

criterion_group!(benches, bench_fastq_pair);
criterion_main!(benches);
