# rsomics-fastq-pair

Re-pair two shuffled paired-end FASTQ files by read name, emitting R1/R2 in
sync — the reads present in both mates, in R1's order.

## Install

```
cargo install rsomics-fastq-pair
```

## Usage

```
rsomics-fastq-pair R1.fq R2.fq --out1 paired_R1.fq --out2 paired_R2.fq

rsomics-fastq-pair shuffled_R1.fastq shuffled_R2.fastq \
  --out1 sync_R1.fastq --out2 sync_R2.fastq
```

- positional `R1` `R2` — the two input mates.
- `--out1` / `--out2` — paired output paths (required).

## Origin

Independent Rust reimplementation of paired-end re-pairing, based on black-box
comparison against `seqkit pair`: the output keeps the reads present in both
mates, in R1's order. This is the same operation provided by `fastq_pair` and
BBTools `repair.sh`.

License: MIT OR Apache-2.0.
Upstream credit: [seqkit](https://github.com/shenwei356/seqkit) (MIT), the
verified oracle; the operation is also provided by fastq_pair and BBTools
repair.sh.
