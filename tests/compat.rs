use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-fastq-pair"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn seqkit_available() -> bool {
    Command::new("seqkit")
        .arg("version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

// Re-pairing must match `seqkit pair`: keep reads present in both R1 and R2
// (in R1's order). Uses matching-ID fixtures (the format seqkit pair expects).
#[test]
fn pair_matches_seqkit() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    let dir = std::env::temp_dir().join("rsomics-fastq-pair-compat");
    let _ = std::fs::create_dir_all(&dir);
    let o1 = dir.join("o1.fq");
    let o2 = dir.join("o2.fq");
    assert!(
        ours()
            .arg(golden("mr1.fq"))
            .arg(golden("mr2.fq"))
            .arg("--out1")
            .arg(&o1)
            .arg("--out2")
            .arg(&o2)
            .status()
            .unwrap()
            .success()
    );

    let sk = dir.join("sk");
    let _ = std::fs::remove_dir_all(&sk);
    assert!(
        Command::new("seqkit")
            .args(["pair", "-1"])
            .arg(golden("mr1.fq"))
            .arg("-2")
            .arg(golden("mr2.fq"))
            .arg("-O")
            .arg(&sk)
            .status()
            .unwrap()
            .success()
    );

    let read = |p: PathBuf| std::fs::read(p).unwrap();
    assert_eq!(read(o1), read(sk.join("mr1.fq")), "R1 mismatch");
    assert_eq!(read(o2), read(sk.join("mr2.fq")), "R2 mismatch");
}
