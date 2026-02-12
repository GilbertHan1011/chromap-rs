# Split Alignment Usage Guide

This guide explains how to use the new split alignment feature in chromap-rs for stitched Hi-C reads.

## Overview

The split alignment feature allows mapping single stitched reads that may align to multiple locations on the genome. This is particularly useful for Hi-C data where reads are stitched together and may contain chimeric sequences from different genomic locations.

## API

### Rust API (chromap-rust)

```rust
use chromap_rust::ChromapAligner;

// Create aligner
let aligner = ChromapAligner::new(
    "reference.index",
    "reference.fa",
    4  // num_threads
)?;

// Prepare stitched sequences
let seqs = vec![
    "ACGTACGTACGTACGTACGTACGTACGTACGT",  // Stitched read 1
    "TGCATGCATGCATGCATGCATGCATGCATGCA",  // Stitched read 2
];
let quals = vec![
    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
];

// Map with split alignment
let alignments = aligner.map_split_batch(&seqs, &quals)?;

// Process results
for aln in alignments {
    println!("Read {}: chr={}, pos={}, strand={:?}, mapq={}, primary={}",
        aln.read_id,
        aln.rid,
        aln.ref_pos,
        aln.strand,
        aln.mapq,
        aln.is_primary
    );
}
```

### Split Alignment Result Structure

```rust
pub struct SplitAlignmentResult {
    pub read_id: u32,       // Which read this belongs to
    pub rid: i32,           // Reference ID (-1 = unmapped)
    pub ref_pos: u32,       // Reference position
    pub strand: Strand,     // Forward or Reverse
    pub query_start: u32,   // Start position in the read
    pub query_end: u32,     // End position in the read
    pub mapq: u8,           // Mapping quality
    pub is_primary: bool,   // Is this the primary alignment?
}
```

## Using in hic-tailor

The `ChromapAligner` in hic-tailor now has an `align_stitched_batch` method:

```rust
use hic_tailor::align::chromap::ChromapAligner;

let aligner = ChromapAligner::new(
    "reference.index",
    "reference.fa",
    4
)?;

// Align stitched fragments
let split_alignments = aligner.align_stitched_batch(&stitched_fragments)?;

// Process split alignments grouped by read
for (read_idx, alns) in split_alignments.iter().enumerate() {
    if alns.is_empty() {
        println!("Read {} has no alignments", read_idx);
        continue;
    }

    // Alignments are sorted by query_start
    let five_prime = alns.first().unwrap();  // 5' end (leftmost in read)
    let three_prime = alns.last().unwrap();  // 3' end (rightmost in read)

    println!("Read {}: 5' at chr{}:{}, 3' at chr{}:{}",
        read_idx,
        five_prime.rid,
        five_prime.ref_pos,
        three_prime.rid,
        three_prime.ref_pos
    );
}
```

## Key Features

1. **Thread-Safe**: The implementation creates per-batch processors on the stack, making it safe to use from multiple threads with Rayon.

2. **Multiple Alignments per Read**: A single stitched read can produce multiple alignment segments, allowing proper handling of chimeric reads.

3. **Query Position Tracking**: Each alignment includes `query_start` and `query_end` to identify which part of the read it corresponds to.

4. **Primary/Supplementary Flags**: Alignments are marked as primary or supplementary, helping identify the best alignment.

## Performance Optimization with Rayon

For parallel processing of large batches:

```rust
use rayon::prelude::*;

let chunk_size = 5000;
let all_alignments: Vec<Vec<SplitAlignmentResult>> = stitched_fragments
    .par_chunks(chunk_size)
    .map(|chunk| {
        aligner.align_stitched_batch(chunk).unwrap_or_default()
    })
    .flatten()
    .collect();
```

## Differences from Paired-End Mapping

| Feature | Paired-End (`map_batch`) | Split Alignment (`map_split_batch`) |
|---------|-------------------------|-------------------------------------|
| Input | R1 and R2 sequences | Single stitched sequence |
| Output | One pair per read | Multiple alignments per read |
| Use Case | Standard Hi-C | Stitched Hi-C reads |
| Chimeric Handling | Limited | Full support |

## Example: Selecting 5' and 3' Alignments

```rust
// Group alignments by read_id
let mut grouped: HashMap<u32, Vec<SplitAlignmentResult>> = HashMap::new();
for aln in alignments {
    grouped.entry(aln.read_id).or_default().push(aln);
}

// For each read, select 5' and 3' alignments
for (read_id, mut alns) in grouped {
    // Sort by query position
    alns.sort_by_key(|a| a.query_start);

    let five_prime = alns.first();  // Leftmost alignment (5' end)
    let three_prime = alns.last();  // Rightmost alignment (3' end)

    // These two alignments form your Hi-C contact
    if let (Some(fp), Some(tp)) = (five_prime, three_prime) {
        println!("Contact: {}:{} <-> {}:{}",
            fp.rid, fp.ref_pos,
            tp.rid, tp.ref_pos
        );
    }
}
```

## Notes

- The split alignment feature is enabled by default in the C++ wrapper
- Buffer size is automatically set to `n_reads * 4` to accommodate multiple alignments per read
- Alignments are returned in the order they are generated, not necessarily sorted by read_id
- For best performance, process reads in batches of 1000-5000
