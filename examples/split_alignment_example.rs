// Example: Using split alignment for stitched Hi-C reads
use chromap_rust::ChromapAligner;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create aligner
    let aligner = ChromapAligner::new(
        "reference.index",
        "reference.fa",
        4  // threads
    )?;

    // Example stitched reads (chimeric sequences)
    let seqs = vec![
        "ACGTACGTACGTACGTACGTACGTACGTACGT",  // Read 0
        "TGCATGCATGCATGCATGCATGCATGCATGCA",  // Read 1
    ];
    let quals = vec![
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    ];

    // Map with split alignment - returns ALL alignments
    let all_alignments = aligner.map_split_batch(&seqs, &quals)?;

    println!("Total alignments: {}", all_alignments.len());

    // Group by read_id
    use std::collections::HashMap;
    let mut grouped: HashMap<u32, Vec<_>> = HashMap::new();
    for aln in all_alignments {
        grouped.entry(aln.read_id).or_default().push(aln);
    }

    // Process each read's alignments
    for (read_id, mut alns) in grouped {
        println!("\nRead {}: {} alignments", read_id, alns.len());

        // Sort by query position
        alns.sort_by_key(|a| a.query_start);

        for (i, aln) in alns.iter().enumerate() {
            println!("  Alignment {}: chr={}, pos={}, strand={:?}, query={}..{}, mapq={}, primary={}",
                i,
                aln.rid,
                aln.ref_pos,
                aln.strand,
                aln.query_start,
                aln.query_end,
                aln.mapq,
                aln.is_primary
            );
        }

        // For Hi-C: select 5' and 3' ends
        if alns.len() >= 2 {
            let five_prime = &alns[0];
            let three_prime = &alns[alns.len() - 1];
            println!("  Hi-C Contact: chr{}:{} <-> chr{}:{}",
                five_prime.rid, five_prime.ref_pos,
                three_prime.rid, three_prime.ref_pos
            );
        }
    }

    Ok(())
}
