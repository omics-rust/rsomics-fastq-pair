use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub struct PairStats {
    pub paired: u64,
    pub singletons: u64,
}

pub fn pair_fastq(
    r1_path: &Path,
    r2_path: &Path,
    out1: &mut dyn Write,
    out2: &mut dyn Write,
    mut singles: Option<&mut dyn Write>,
) -> Result<PairStats> {
    // Single pass over R2: name -> its 4-line record.
    let file2 = File::open(r2_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", r2_path.display())))?;
    let mut r2: HashMap<String, [String; 4]> = HashMap::new();
    let mut lines2 = BufReader::new(file2).lines();
    while let Some(h) = lines2.next() {
        let h = h.map_err(RsomicsError::Io)?;
        let s = next_line(&mut lines2)?;
        let p = next_line(&mut lines2)?;
        let q = next_line(&mut lines2)?;
        let name = read_name(&h);
        r2.insert(name, [h, s, p, q]);
    }

    let file1 = File::open(r1_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", r1_path.display())))?;
    let mut o1 = BufWriter::with_capacity(256 * 1024, out1);
    let mut o2 = BufWriter::with_capacity(256 * 1024, out2);
    let mut stats = PairStats {
        paired: 0,
        singletons: 0,
    };

    let mut lines1 = BufReader::new(file1).lines();
    while let Some(header) = lines1.next() {
        let header = header.map_err(RsomicsError::Io)?;
        let seq = next_line(&mut lines1)?;
        let plus = next_line(&mut lines1)?;
        let qual = next_line(&mut lines1)?;

        let name = read_name(&header);
        if let Some(mate) = r2.get(&name) {
            writeln!(o1, "{header}\n{seq}\n{plus}\n{qual}").map_err(RsomicsError::Io)?;
            writeln!(o2, "{}\n{}\n{}\n{}", mate[0], mate[1], mate[2], mate[3])
                .map_err(RsomicsError::Io)?;
            stats.paired += 1;
        } else if let Some(s) = singles.as_deref_mut() {
            writeln!(s, "{header}\n{seq}\n{plus}\n{qual}").map_err(RsomicsError::Io)?;
            stats.singletons += 1;
        }
    }

    o1.flush().map_err(RsomicsError::Io)?;
    o2.flush().map_err(RsomicsError::Io)?;
    Ok(stats)
}

fn read_name(header: &str) -> String {
    header
        .split_once(|c: char| c.is_whitespace() || c == '/')
        .map_or(header, |(name, _)| name)
        .trim_start_matches('@')
        .to_string()
}

fn next_line<B: BufRead>(lines: &mut std::io::Lines<B>) -> Result<String> {
    lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("truncated FASTQ".into()))?
        .map_err(RsomicsError::Io)
}
