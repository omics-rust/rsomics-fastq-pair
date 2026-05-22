use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub struct PairStats {
    pub paired: u64,
    pub singletons: u64,
}

/// Read ID = bytes after a leading `@`, up to the first whitespace or `/`.
fn read_name(line: &[u8]) -> &[u8] {
    let l = line.strip_prefix(b"@").unwrap_or(line);
    let end = l
        .iter()
        .position(|&c| c == b' ' || c == b'\t' || c == b'/' || c == b'\r')
        .unwrap_or(l.len());
    &l[..end]
}

pub fn pair_fastq(
    r1_path: &Path,
    r2_path: &Path,
    out1: &mut dyn Write,
    out2: &mut dyn Write,
    mut singles: Option<&mut dyn Write>,
) -> Result<PairStats> {
    // R2 is read once into one buffer; the map borrows name + full-record byte
    // slices out of it (zero per-record allocation). R1 is streamed.
    let r2_buf = fs::read(r2_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", r2_path.display())))?;
    let mut map: HashMap<&[u8], &[u8]> = HashMap::new();
    {
        let n = r2_buf.len();
        let mut i = 0;
        let mut line_no = 0u8;
        let mut rec_start = 0;
        let mut name: &[u8] = &[];
        while i < n {
            let nl = r2_buf[i..]
                .iter()
                .position(|&c| c == b'\n')
                .map_or(n, |p| i + p);
            if line_no == 0 {
                name = read_name(&r2_buf[i..nl]);
                rec_start = i;
            }
            line_no += 1;
            i = nl + 1;
            if line_no == 4 {
                map.insert(name, &r2_buf[rec_start..i.min(n)]);
                line_no = 0;
            }
        }
    }

    let reader = BufReader::new(
        File::open(r1_path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", r1_path.display())))?,
    );
    let mut reader = reader;
    let mut o1 = BufWriter::with_capacity(256 * 1024, out1);
    let mut o2 = BufWriter::with_capacity(256 * 1024, out2);
    let mut stats = PairStats {
        paired: 0,
        singletons: 0,
    };

    let mut buf: Vec<u8> = Vec::with_capacity(512);
    loop {
        buf.clear();
        let mut lines_read = 0;
        for _ in 0..4 {
            if reader
                .read_until(b'\n', &mut buf)
                .map_err(RsomicsError::Io)?
                == 0
            {
                break;
            }
            lines_read += 1;
        }
        if lines_read == 0 {
            break;
        }
        if lines_read < 4 {
            return Err(RsomicsError::InvalidInput("truncated FASTQ".into()));
        }
        let line0_end = buf.iter().position(|&c| c == b'\n').unwrap_or(buf.len());
        let name = read_name(&buf[..line0_end]);
        if let Some(mate) = map.get(name) {
            o1.write_all(&buf).map_err(RsomicsError::Io)?;
            o2.write_all(mate).map_err(RsomicsError::Io)?;
            stats.paired += 1;
        } else if let Some(s) = singles.as_deref_mut() {
            s.write_all(&buf).map_err(RsomicsError::Io)?;
            stats.singletons += 1;
        }
    }

    o1.flush().map_err(RsomicsError::Io)?;
    o2.flush().map_err(RsomicsError::Io)?;
    Ok(stats)
}
