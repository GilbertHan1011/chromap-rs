#pragma once

#include <cstdint>
#include <cstddef>

#ifdef __cplusplus
extern "C" {
#endif

// Opaque pointer for the persistent index/mapper object
typedef struct ChromapMapper ChromapMapper;

// Output struct for mapping results
struct RustMapping {
    uint32_t read_id;
    uint32_t ref_id;
    uint32_t ref_pos;
    uint8_t mapq;
    int8_t strand;  // 0: +, 1: -
    uint16_t num_dups;
    uint16_t is_unique;
    uint16_t is_primary;
};

// Create a new ChromapMapper instance
// index_path: path to the chromap index file
// ref_path: path to the reference FASTA file
// num_threads: number of threads to use for mapping
ChromapMapper* mapper_new(const char* index_path, const char* ref_path, int num_threads);

// Free the ChromapMapper instance
void mapper_free(ChromapMapper* mapper);

// Map a batch of paired-end reads
// mapper: the ChromapMapper instance
// seqs1: array of pointers to R1 DNA strings (null-terminated)
// seqs2: array of pointers to R2 DNA strings (null-terminated)
// quals1: array of pointers to R1 quality strings (null-terminated)
// quals2: array of pointers to R2 quality strings (null-terminated)
// n_reads: number of read pairs in the batch
// out_buffer: pointer to output buffer (allocated by caller)
// buffer_size: size of the output buffer
// Returns: number of mappings written to out_buffer
int mapper_map_batch(ChromapMapper* mapper,
                     const char** seqs1, const char** seqs2,
                     const char** quals1, const char** quals2,
                     int n_reads,
                     RustMapping* out_buffer,
                     int buffer_size);

#ifdef __cplusplus
}
#endif
