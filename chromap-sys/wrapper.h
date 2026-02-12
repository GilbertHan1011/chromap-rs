#pragma once

#include <cstdint>
#include <cstddef>

#ifdef __cplusplus
extern "C" {
#endif

// Opaque pointer for the persistent index/mapper object
typedef struct ChromapMapper ChromapMapper;

// Output struct for Hi-C pairs mapping results
// Matches Chromap's PairsMapping format for Hi-C data
struct RustPairsMapping {
    uint32_t read_id;
    int32_t rid1;      // Reference ID for first end
    int32_t rid2;      // Reference ID for second end
    uint32_t pos1;     // Position for first end
    uint32_t pos2;     // Position for second end
    int32_t strand1;   // 1: +, 0: -
    int32_t strand2;   // 1: +, 0: -
    uint8_t mapq;
    uint8_t is_unique;
    uint8_t num_dups;
};

// Create a new ChromapMapper instance
// index_path: path to the chromap index file
// ref_path: path to the reference FASTA file
// num_threads: number of threads to use for mapping
// min_num_seeds: minimum number of seeds for mapping (-s, use -1 for default: 2)
// min_read_length: minimum read length (--min-read-length, use -1 for default: 30)
ChromapMapper* mapper_new(const char* index_path, const char* ref_path, int num_threads,
                          int min_num_seeds, int min_read_length);

// Free the ChromapMapper instance
void mapper_free(ChromapMapper* mapper);

// Map a batch of paired-end reads in Hi-C mode
// mapper: the ChromapMapper instance
// seqs1: array of pointers to R1 DNA strings (null-terminated)
// seqs2: array of pointers to R2 DNA strings (null-terminated)
// quals1: array of pointers to R1 quality strings (null-terminated)
// quals2: array of pointers to R2 quality strings (null-terminated)
// n_reads: number of read pairs in the batch
// out_buffer: pointer to output buffer (allocated by caller)
// buffer_size: size of the output buffer
// Returns: number of pairs mappings written to out_buffer (one mapping per pair)
int mapper_map_batch(ChromapMapper* mapper,
                     const char** seqs1, const char** seqs2,
                     const char** quals1, const char** quals2,
                     int n_reads,
                     RustPairsMapping* out_buffer,
                     int buffer_size);

// Get the number of reference sequences
int mapper_get_num_references(ChromapMapper* mapper);

// Get the name of a reference sequence by index
const char* mapper_get_reference_name(ChromapMapper* mapper, int index);

// Get the length of a reference sequence by index
int mapper_get_reference_length(ChromapMapper* mapper, int index);

// Output struct for split alignment results (single-end reads with multiple alignments)
struct RustSplitAlignment {
    uint32_t read_id;       // Which read this belongs to
    int32_t rid;            // Reference ID (-1 = unmapped)
    uint32_t ref_pos;       // Reference start position
    int32_t strand;         // 0: +, 1: -
    uint32_t query_start;   // Where in the read this alignment starts
    uint32_t query_end;     // Where in the read this alignment ends
    uint8_t mapq;
    uint8_t is_primary;     // Is this the "best" alignment?
};

// Map a batch of single-end reads with split alignment support
// mapper: the ChromapMapper instance
// seqs: array of pointers to DNA strings (null-terminated)
// quals: array of pointers to quality strings (null-terminated)
// n_reads: number of reads in the batch
// out_buffer: pointer to output buffer (allocated by caller)
// buffer_size: size of the output buffer (should be n_reads * max_splits_per_read)
// Returns: number of split alignments written to out_buffer (may be > n_reads)
int mapper_map_split_batch(ChromapMapper* mapper,
                           const char** seqs,
                           const char** quals,
                           int n_reads,
                           RustSplitAlignment* out_buffer,
                           int buffer_size);

#ifdef __cplusplus
}
#endif
