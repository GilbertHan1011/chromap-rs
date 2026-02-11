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
ChromapMapper* mapper_new(const char* index_path, const char* ref_path, int num_threads);

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

#ifdef __cplusplus
}
#endif
