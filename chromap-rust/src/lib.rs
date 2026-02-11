use chromap_sys::{ChromapMapper, RustMapping};
use std::ffi::{CString, NulError};
use std::path::Path;

/// Error types for chromap operations
#[derive(Debug)]
pub enum ChromapError {
    /// Failed to create CString (contains null byte)
    NullByteError(NulError),
    /// Failed to initialize mapper
    InitializationError,
    /// Mapping operation failed
    MappingError,
    /// Invalid input
    InvalidInput(String),
}

impl From<NulError> for ChromapError {
    fn from(err: NulError) -> Self {
        ChromapError::NullByteError(err)
    }
}

impl std::fmt::Display for ChromapError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ChromapError::NullByteError(e) => write!(f, "Null byte error: {}", e),
            ChromapError::InitializationError => write!(f, "Failed to initialize mapper"),
            ChromapError::MappingError => write!(f, "Mapping operation failed"),
            ChromapError::InvalidInput(msg) => write!(f, "Invalid input: {}", msg),
        }
    }
}

impl std::error::Error for ChromapError {}

/// A mapping result for a single read
#[derive(Debug, Clone)]
pub struct MappingResult {
    pub read_id: u32,
    pub ref_id: u32,
    pub ref_pos: u32,
    pub mapq: u8,
    pub strand: Strand,
    pub num_dups: u16,
    pub is_unique: bool,
    pub is_primary: bool,
}

/// Strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl From<RustMapping> for MappingResult {
    fn from(mapping: RustMapping) -> Self {
        MappingResult {
            read_id: mapping.read_id,
            ref_id: mapping.ref_id,
            ref_pos: mapping.ref_pos,
            mapq: mapping.mapq,
            strand: if mapping.strand == 0 {
                Strand::Forward
            } else {
                Strand::Reverse
            },
            num_dups: mapping.num_dups,
            is_unique: mapping.is_unique != 0,
            is_primary: mapping.is_primary != 0,
        }
    }
}

/// High-level chromap aligner interface
pub struct ChromapAligner {
    inner: *mut ChromapMapper,
}

impl ChromapAligner {
    /// Create a new ChromapAligner
    ///
    /// # Arguments
    /// * `index_path` - Path to the chromap index file
    /// * `ref_path` - Path to the reference FASTA file
    /// * `num_threads` - Number of threads to use for mapping
    ///
    /// # Example
    /// ```no_run
    /// use chromap_rust::ChromapAligner;
    ///
    /// let aligner = ChromapAligner::new(
    ///     "reference.index",
    ///     "reference.fa",
    ///     4
    /// ).expect("Failed to create aligner");
    /// ```
    pub fn new<P1: AsRef<Path>, P2: AsRef<Path>>(
        index_path: P1,
        ref_path: P2,
        num_threads: usize,
    ) -> Result<Self, ChromapError> {
        let index_path_str = index_path
            .as_ref()
            .to_str()
            .ok_or_else(|| ChromapError::InvalidInput("Invalid index path".to_string()))?;
        let ref_path_str = ref_path
            .as_ref()
            .to_str()
            .ok_or_else(|| ChromapError::InvalidInput("Invalid reference path".to_string()))?;

        let c_index_path = CString::new(index_path_str)?;
        let c_ref_path = CString::new(ref_path_str)?;

        unsafe {
            let inner = chromap_sys::mapper_new(
                c_index_path.as_ptr(),
                c_ref_path.as_ptr(),
                num_threads as i32,
            );

            if inner.is_null() {
                return Err(ChromapError::InitializationError);
            }

            Ok(ChromapAligner { inner })
        }
    }

    /// Map a batch of paired-end reads
    ///
    /// # Arguments
    /// * `r1_seqs` - R1 sequences (DNA strings)
    /// * `r2_seqs` - R2 sequences (DNA strings)
    /// * `r1_quals` - R1 quality strings
    /// * `r2_quals` - R2 quality strings
    ///
    /// # Returns
    /// Vector of mapping results
    ///
    /// # Example
    /// ```no_run
    /// # use chromap_rust::ChromapAligner;
    /// # let aligner = ChromapAligner::new("ref.index", "ref.fa", 4).unwrap();
    /// let r1_seqs = vec!["ACGTACGT", "TGCATGCA"];
    /// let r2_seqs = vec!["GCTAGCTA", "ATCGATCG"];
    /// let r1_quals = vec!["IIIIIIII", "IIIIIIII"];
    /// let r2_quals = vec!["IIIIIIII", "IIIIIIII"];
    ///
    /// let results = aligner.map_batch(&r1_seqs, &r2_seqs, &r1_quals, &r2_quals)
    ///     .expect("Mapping failed");
    /// ```
    pub fn map_batch(
        &self,
        r1_seqs: &[&str],
        r2_seqs: &[&str],
        r1_quals: &[&str],
        r2_quals: &[&str],
    ) -> Result<Vec<MappingResult>, ChromapError> {
        if r1_seqs.len() != r2_seqs.len()
            || r1_seqs.len() != r1_quals.len()
            || r1_seqs.len() != r2_quals.len()
        {
            return Err(ChromapError::InvalidInput(
                "All input arrays must have the same length".to_string(),
            ));
        }

        let n = r1_seqs.len();
        if n == 0 {
            return Ok(Vec::new());
        }

        // Convert Rust strings to C strings
        let c_r1_seqs: Vec<CString> = r1_seqs
            .iter()
            .map(|s| CString::new(*s))
            .collect::<Result<Vec<_>, _>>()?;
        let c_r2_seqs: Vec<CString> = r2_seqs
            .iter()
            .map(|s| CString::new(*s))
            .collect::<Result<Vec<_>, _>>()?;
        let c_r1_quals: Vec<CString> = r1_quals
            .iter()
            .map(|s| CString::new(*s))
            .collect::<Result<Vec<_>, _>>()?;
        let c_r2_quals: Vec<CString> = r2_quals
            .iter()
            .map(|s| CString::new(*s))
            .collect::<Result<Vec<_>, _>>()?;

        // Create pointer arrays
        let ptrs_r1: Vec<*const i8> = c_r1_seqs.iter().map(|c| c.as_ptr()).collect();
        let ptrs_r2: Vec<*const i8> = c_r2_seqs.iter().map(|c| c.as_ptr()).collect();
        let ptrs_q1: Vec<*const i8> = c_r1_quals.iter().map(|c| c.as_ptr()).collect();
        let ptrs_q2: Vec<*const i8> = c_r2_quals.iter().map(|c| c.as_ptr()).collect();

        // Allocate output buffer
        let mut results = vec![
            RustMapping {
                read_id: 0,
                ref_id: 0,
                ref_pos: 0,
                mapq: 0,
                strand: 0,
                num_dups: 0,
                is_unique: 0,
                is_primary: 0,
            };
            n
        ];

        unsafe {
            let ret = chromap_sys::mapper_map_batch(
                self.inner,
                ptrs_r1.as_ptr(),
                ptrs_r2.as_ptr(),
                ptrs_q1.as_ptr(),
                ptrs_q2.as_ptr(),
                n as i32,
                results.as_mut_ptr(),
                n as i32,
            );

            if ret < 0 {
                return Err(ChromapError::MappingError);
            }

            // Convert to high-level results
            Ok(results.into_iter().map(|r| r.into()).collect())
        }
    }
}

impl Drop for ChromapAligner {
    fn drop(&mut self) {
        unsafe {
            chromap_sys::mapper_free(self.inner);
        }
    }
}

// Safety: ChromapMapper is thread-safe (assuming proper internal synchronization)
unsafe impl Send for ChromapAligner {}
unsafe impl Sync for ChromapAligner {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_conversion() {
        let mapping = RustMapping {
            read_id: 0,
            ref_id: 1,
            ref_pos: 100,
            mapq: 60,
            strand: 0,
            num_dups: 0,
            is_unique: 1,
            is_primary: 1,
        };

        let result: MappingResult = mapping.into();
        assert_eq!(result.strand, Strand::Forward);
    }
}

