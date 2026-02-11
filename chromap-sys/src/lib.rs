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

/// Mapping result structure matching the C++ RustMapping struct
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
    pub fn mapper_new(
        index_path: *const c_char,
        ref_path: *const c_char,
        num_threads: c_int,
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
        out_buffer: *mut RustMapping,
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
    }
}
