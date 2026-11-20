//! Databases of the tools that read hits: every sequence of a FASTA, FASTQ
//! or UDB file (known by its magic number), with its header, in order, so
//! that the target numbers of UC hits find them.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use crate::fastx::{Reader, Record};
use crate::udb::{self, Udb};

/// reads the (header, sequence) pairs of the database at `path`.
pub fn read(path: &str) -> Result<Vec<(String, Vec<u8>)>, String> {
    let file = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    let size = file
        .metadata()
        .ok()
        .filter(|m| m.is_file())
        .map(|m| m.len());
    let mut fin = BufReader::new(file);
    let start = fin
        .fill_buf()
        .map_err(|e| format!("cannot read {path}: {e}"))?;
    if start.starts_with(udb::MAGIC) {
        let ignore = |_: &[u32], _: u64, _: &[u32]| {};
        let read = match size {
            Some(size) => Udb::read(&mut fin, size, ignore),
            None => {
                let mut bytes = Vec::new();
                fin.read_to_end(&mut bytes)
                    .map_err(|e| format!("cannot read {path}: {e}"))?;
                Udb::read(&mut &bytes[..], bytes.len() as u64, ignore)
            }
        };
        let udb = read.map_err(|e| e.to_string())?;
        let mut at = 0;
        return Ok(udb
            .header_starts
            .windows(2)
            .zip(&udb.lengths)
            .map(|(pair, &length)| {
                let header = &udb.headers[pair[0] as usize..pair[1] as usize - 1];
                let seq = udb.sequences[at..at + length as usize].to_vec();
                at += length as usize;
                (String::from_utf8_lossy(header).into_owned(), seq)
            })
            .collect());
    }
    let mut reader = Reader::new(fin);
    let mut record = Record::new();
    let mut db = Vec::new();
    loop {
        reader.read(&mut record).map_err(|e| e.to_string())?;
        if record.is_empty() {
            return Ok(db);
        }
        db.push((record.header().to_string(), record.seq().to_vec()));
    }
}
