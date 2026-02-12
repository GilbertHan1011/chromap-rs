#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]

use std::os::raw::{c_char, c_int};

/// Opaque pointer for the ChromapMapper C++ object
#[repr(C)]
pub struct ChromapMapper {
    _private: [u8; 0],
}

/// Legacy single-end mapping result structure (no longer used in Hi-C mode)
/// Kept for backwards compatibility and tests.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct RustMapping {
    pub read_id: u32,
    pub ref_id: u32,
    pub ref_pos: u32,
    pub mapq: u8,
    pub strand: i8,
    pub num_dups: u16,
    pub is_unique: u16,
    pub is_primary: u16,
}

/// Hi-C pairs mapping result structure matching the C++ RustPairsMapping struct.
/// This is the primary mapping struct used by the current integration.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct RustPairsMapping {
    pub read_id: u32,
    pub rid1: i32,
    pub rid2: i32,
    pub pos1: u32,
    pub pos2: u32,
    pub strand1: i32,
    pub strand2: i32,
    pub mapq: u8,
    pub is_unique: u8,
    pub num_dups: u8,
}

/// Split alignment result structure for single-end reads with multiple alignments.
/// Used for stitched Hi-C reads that may map to multiple locations.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct RustSplitAlignment {
    pub read_id: u32,
    pub rid: i32,
    pub ref_pos: u32,
    pub strand: i32,
    pub query_start: u32,
    pub query_end: u32,
    pub mapq: u8,
    pub is_primary: u8,
}

extern "C" {
    /// Create a new ChromapMapper instance
    ///
    /// # Arguments
    /// * `index_path` - Path to the chromap index file
    /// * `ref_path` - Path to the reference FASTA file
    /// * `num_threads` - Number of threads to use for mapping
    ///
    /// # Returns
    /// Pointer to ChromapMapper or null on failure
    ///
    /// # Optional parameters (use -1 for default)
    /// * `min_num_seeds` - minimum number of seeds for mapping (-s), default 2
    /// * `min_read_length` - minimum read length (--min-read-length), default 30
    pub fn mapper_new(
        index_path: *const c_char,
        ref_path: *const c_char,
        num_threads: c_int,
        min_num_seeds: c_int,
        min_read_length: c_int,
    ) -> *mut ChromapMapper;

    /// Free a ChromapMapper instance
    ///
    /// # Arguments
    /// * `mapper` - Pointer to ChromapMapper to free
    pub fn mapper_free(mapper: *mut ChromapMapper);

    /// Map a batch of paired-end reads
    ///
    /// # Arguments
    /// * `mapper` - Pointer to ChromapMapper instance
    /// * `seqs1` - Array of pointers to R1 DNA strings
    /// * `seqs2` - Array of pointers to R2 DNA strings
    /// * `quals1` - Array of pointers to R1 quality strings
    /// * `quals2` - Array of pointers to R2 quality strings
    /// * `n_reads` - Number of read pairs
    /// * `out_buffer` - Output buffer for mappings
    /// * `buffer_size` - Size of output buffer
    ///
    /// # Returns
    /// Number of mappings written or -1 on error
    pub fn mapper_map_batch(
        mapper: *mut ChromapMapper,
        seqs1: *const *const c_char,
        seqs2: *const *const c_char,
        quals1: *const *const c_char,
        quals2: *const *const c_char,
        n_reads: c_int,
        out_buffer: *mut RustPairsMapping,
        buffer_size: c_int,
    ) -> c_int;

    /// Get the number of reference sequences
    ///
    /// # Arguments
    /// * `mapper` - Pointer to ChromapMapper instance
    ///
    /// # Returns
    /// Number of reference sequences or -1 on error
    pub fn mapper_get_num_references(mapper: *mut ChromapMapper) -> c_int;

    /// Get the name of a reference sequence by index
    ///
    /// # Arguments
    /// * `mapper` - Pointer to ChromapMapper instance
    /// * `index` - Index of the reference sequence
    ///
    /// # Returns
    /// Pointer to null-terminated string or null on error
    pub fn mapper_get_reference_name(mapper: *mut ChromapMapper, index: c_int) -> *const c_char;

    /// Get the length of a reference sequence by index
    ///
    /// # Arguments
    /// * `mapper` - Pointer to ChromapMapper instance
    /// * `index` - Index of the reference sequence
    ///
    /// # Returns
    /// Length of the reference sequence or -1 on error
    pub fn mapper_get_reference_length(mapper: *mut ChromapMapper, index: c_int) -> c_int;

    /// Map a batch of single-end reads with split alignment support
    ///
    /// # Arguments
    /// * `mapper` - Pointer to ChromapMapper instance
    /// * `seqs` - Array of pointers to DNA strings
    /// * `quals` - Array of pointers to quality strings
    /// * `n_reads` - Number of reads
    /// * `out_buffer` - Output buffer for split alignments
    /// * `buffer_size` - Size of output buffer (should be n_reads * max_splits_per_read)
    ///
    /// # Returns
    /// Number of split alignments written or -1 on error
    pub fn mapper_map_split_batch(
        mapper: *mut ChromapMapper,
        seqs: *const *const c_char,
        quals: *const *const c_char,
        n_reads: c_int,
        out_buffer: *mut RustSplitAlignment,
        buffer_size: c_int,
    ) -> c_int;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_struct_sizes() {
        // Verify struct sizes match expectations
        assert_eq!(std::mem::size_of::<RustMapping>(), 20);
        // Sanity-check Hi-C pairs mapping struct size as well
        assert!(
            std::mem::size_of::<RustPairsMapping>() >= 32,
            "RustPairsMapping size unexpectedly small"
        );
    }
}
