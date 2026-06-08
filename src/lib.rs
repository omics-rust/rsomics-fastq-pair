use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
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

/// FNV-1a over the read name. Names are short, so this is cheaper than siphash
/// and avoids storing the name bytes for the divergent (shuffled) index path.
fn name_hash(name: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in name {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// One FASTQ record's raw bytes (the four lines, newline-terminated as on disk)
/// plus its name, read from a buffered stream.
struct Record {
    bytes: Vec<u8>,
    name_end: usize,
    name_start: usize,
}

impl Record {
    fn name(&self) -> &[u8] {
        &self.bytes[self.name_start..self.name_end]
    }
}

/// Pull the next 4-line record off `reader`, appending its bytes. Returns `None`
/// at clean EOF, errors on a truncated record.
fn next_record<R: BufRead>(reader: &mut R) -> Result<Option<Record>> {
    let mut bytes = Vec::with_capacity(512);
    let mut lines = 0;
    let mut name_start = 0;
    let mut name_end = 0;
    for li in 0..4 {
        let before = bytes.len();
        if reader
            .read_until(b'\n', &mut bytes)
            .map_err(RsomicsError::Io)?
            == 0
        {
            break;
        }
        if li == 0 {
            let line0 = &bytes[before..];
            let line0_end = line0
                .iter()
                .position(|&c| c == b'\n')
                .unwrap_or(line0.len());
            let name = read_name(&line0[..line0_end]);
            name_start = before + (name.as_ptr() as usize - line0.as_ptr() as usize);
            name_end = name_start + name.len();
        }
        lines += 1;
    }
    if lines == 0 {
        return Ok(None);
    }
    if lines < 4 {
        return Err(RsomicsError::InvalidInput("truncated FASTQ".into()));
    }
    Ok(Some(Record {
        bytes,
        name_end,
        name_start,
    }))
}

pub fn pair_fastq(
    r1_path: &Path,
    r2_path: &Path,
    out1: &mut dyn Write,
    out2: &mut dyn Write,
    mut singles: Option<&mut dyn Write>,
) -> Result<PairStats> {
    let mut r1 = BufReader::with_capacity(
        256 * 1024,
        File::open(r1_path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", r1_path.display())))?,
    );
    let mut r2 = BufReader::with_capacity(
        256 * 1024,
        File::open(r2_path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", r2_path.display())))?,
    );
    let mut o1 = BufWriter::with_capacity(256 * 1024, out1);
    let mut o2 = BufWriter::with_capacity(256 * 1024, out2);
    let mut stats = PairStats {
        paired: 0,
        singletons: 0,
    };

    // Lockstep pass: while R1 and R2 sit in the same order, each mate is the
    // record the other reader is already pointing at, so we emit straight
    // through — no name index, no buffering, RSS stays at the I/O buffers.
    // The shuffled case is rare in practice; the first mismatch hands off to
    // the offset-index path below.
    let mut pending_r1: Option<Record> = None;
    while let Some(rec1) = next_record(&mut r1)? {
        let Some(rec2) = next_record(&mut r2)? else {
            // R2 ran out first — every remaining R1 read is mateless, so the
            // tail drains straight to singletons below.
            pending_r1 = Some(rec1);
            break;
        };
        if rec1.name() == rec2.name() {
            o1.write_all(&rec1.bytes).map_err(RsomicsError::Io)?;
            o2.write_all(&rec2.bytes).map_err(RsomicsError::Io)?;
            stats.paired += 1;
        } else {
            return drain_shuffled(
                rec1,
                &mut r1,
                r2_path,
                &mut o1,
                &mut o2,
                &mut singles,
                stats,
            );
        }
    }

    // R1 may still hold an unmatched tail when R2 ended early.
    if let Some(rec1) = pending_r1 {
        emit_or_single(&rec1, &mut singles, &mut stats)?;
        while let Some(rec1) = next_record(&mut r1)? {
            emit_or_single(&rec1, &mut singles, &mut stats)?;
        }
    }

    o1.flush().map_err(RsomicsError::Io)?;
    o2.flush().map_err(RsomicsError::Io)?;
    Ok(stats)
}

/// Write `rec` to the singletons file if one is configured, else drop it.
fn emit_or_single(
    rec: &Record,
    singles: &mut Option<&mut dyn Write>,
    stats: &mut PairStats,
) -> Result<()> {
    if let Some(s) = singles.as_deref_mut() {
        s.write_all(&rec.bytes).map_err(RsomicsError::Io)?;
        stats.singletons += 1;
    }
    Ok(())
}

/// R1 and R2 stopped lining up at `rec1`. Index R2 by name-hash → (offset, len),
/// then stream the rest of R1, reading each mate's R2 record back from disk on a
/// hit. Storing the hash + offset + length — not the record — keeps memory in
/// seqkit's class even when every read is shuffled.
#[allow(clippy::too_many_arguments)]
fn drain_shuffled<R1: BufRead, W1: Write, W2: Write>(
    rec1: Record,
    r1: &mut R1,
    r2_path: &Path,
    o1: &mut W1,
    o2: &mut W2,
    singles: &mut Option<&mut dyn Write>,
    mut stats: PairStats,
) -> Result<PairStats> {
    // R2's BufReader already over-read past the divergence point, so its file
    // cursor isn't a record boundary. Re-open R2 and index it from the top; the
    // handful of records already emitted in lockstep are simply re-found by hash.
    let mut r2_file = File::open(r2_path).map_err(RsomicsError::Io)?;
    let mut index: HashMap<u64, (u64, u32)> = HashMap::new();
    {
        let mut scan = BufReader::with_capacity(256 * 1024, &mut r2_file);
        let mut offset: u64 = 0;
        while let Some(rec) = next_record(&mut scan)? {
            let len = rec.bytes.len() as u32;
            index.insert(name_hash(rec.name()), (offset, len));
            offset += rec.bytes.len() as u64;
        }
    }

    let mut mate_buf: Vec<u8> = Vec::with_capacity(512);
    let mut emit_pair = |rec: &Record,
                         o1: &mut W1,
                         o2: &mut W2,
                         singles: &mut Option<&mut dyn Write>,
                         stats: &mut PairStats|
     -> Result<()> {
        match index.get(&name_hash(rec.name())) {
            Some(&(off, len)) => {
                mate_buf.resize(len as usize, 0);
                r2_file
                    .seek(SeekFrom::Start(off))
                    .map_err(RsomicsError::Io)?;
                r2_file
                    .read_exact(&mut mate_buf)
                    .map_err(RsomicsError::Io)?;
                let mate_line0 = &mate_buf[..mate_buf
                    .iter()
                    .position(|&c| c == b'\n')
                    .unwrap_or(mate_buf.len())];
                if read_name(mate_line0) != rec.name() {
                    return Err(RsomicsError::InvalidInput(
                        "read-name hash collision while pairing FASTQ".into(),
                    ));
                }
                o1.write_all(&rec.bytes).map_err(RsomicsError::Io)?;
                o2.write_all(&mate_buf).map_err(RsomicsError::Io)?;
                stats.paired += 1;
            }
            None => {
                if let Some(s) = singles.as_deref_mut() {
                    s.write_all(&rec.bytes).map_err(RsomicsError::Io)?;
                    stats.singletons += 1;
                }
            }
        }
        Ok(())
    };

    emit_pair(&rec1, o1, o2, singles, &mut stats)?;
    while let Some(rec) = next_record(r1)? {
        emit_pair(&rec, o1, o2, singles, &mut stats)?;
    }

    o1.flush().map_err(RsomicsError::Io)?;
    o2.flush().map_err(RsomicsError::Io)?;
    Ok(stats)
}
