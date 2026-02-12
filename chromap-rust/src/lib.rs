use chromap_sys::{ChromapMapper, RustPairsMapping};
use std::ffi::{CStr, CString, NulError};
use std::path::Path;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::fastq;

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

/// A Hi-C pairs mapping result containing both ends of the pair
#[derive(Debug, Clone)]
pub struct MappingResult {
    pub read_id: u32,
    pub rid1: i32,      // Reference ID for end 1 (-1 = unmapped)
    pub rid2: i32,      // Reference ID for end 2 (-1 = unmapped)
    pub pos1: u32,      // Position for end 1
    pub pos2: u32,      // Position for end 2
    pub strand1: Strand, // Strand for end 1
    pub strand2: Strand, // Strand for end 2
    pub mapq: u8,
    pub is_unique: bool,
    pub num_dups: u8,
}

/// Strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl From<RustPairsMapping> for MappingResult {
    fn from(mapping: RustPairsMapping) -> Self {
        MappingResult {
            read_id: mapping.read_id,
            rid1: mapping.rid1,
            rid2: mapping.rid2,
            pos1: mapping.pos1,
            pos2: mapping.pos2,
            strand1: if mapping.strand1 == 1 {
                Strand::Forward
            } else {
                Strand::Reverse
            },
            strand2: if mapping.strand2 == 1 {
                Strand::Forward
            } else {
                Strand::Reverse
            },
            mapq: mapping.mapq,
            is_unique: mapping.is_unique != 0,
            num_dups: mapping.num_dups,
        }
    }
}

/// High-level chromap aligner interface
pub struct ChromapAligner {
    inner: *mut ChromapMapper,
    header: sam::Header,
    contig_names: Vec<String>,
    contig_lengths: Vec<usize>,
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
        Self::new_with_options(index_path, ref_path, num_threads, None, None)
    }

    /// Create a new ChromapAligner with optional mapping parameters.
    pub fn new_with_options<P1: AsRef<Path>, P2: AsRef<Path>>(
        index_path: P1,
        ref_path: P2,
        num_threads: usize,
        min_num_seeds: Option<u32>,
        min_read_length: Option<u32>,
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

        let min_num_seeds_c = min_num_seeds.map(|v| v as i32).unwrap_or(-1);
        let min_read_length_c = min_read_length.map(|v| v as i32).unwrap_or(-1);

        unsafe {
            let inner = chromap_sys::mapper_new(
                c_index_path.as_ptr(),
                c_ref_path.as_ptr(),
                num_threads as i32,
                min_num_seeds_c,
                min_read_length_c,
            );

            if inner.is_null() {
                return Err(ChromapError::InitializationError);
            }

            // Get reference sequence information
            let num_refs = chromap_sys::mapper_get_num_references(inner);
            if num_refs < 0 {
                chromap_sys::mapper_free(inner);
                return Err(ChromapError::InitializationError);
            }

            let mut contig_names = Vec::new();
            let mut contig_lengths = Vec::new();

            for i in 0..num_refs {
                let name_ptr = chromap_sys::mapper_get_reference_name(inner, i);
                if name_ptr.is_null() {
                    chromap_sys::mapper_free(inner);
                    return Err(ChromapError::InitializationError);
                }
                let name = CStr::from_ptr(name_ptr)
                    .to_str()
                    .map_err(|_| ChromapError::InvalidInput("Invalid reference name".to_string()))?
                    .to_string();

                let length = chromap_sys::mapper_get_reference_length(inner, i);
                if length < 0 {
                    chromap_sys::mapper_free(inner);
                    return Err(ChromapError::InitializationError);
                }

                contig_names.push(name);
                contig_lengths.push(length as usize);
            }

            // Create SAM header
            let ref_seqs = contig_names
                .iter()
                .zip(contig_lengths.iter())
                .map(|(name, len)| {
                    (
                        bstr::BString::from(name.as_str()),
                        Map::<ReferenceSequence>::new(std::num::NonZeroUsize::try_from(*len).unwrap()),
                    )
                })
                .collect();
            let header = sam::Header::builder()
                .set_reference_sequences(ref_seqs)
                .build();

            Ok(ChromapAligner {
                inner,
                header,
                contig_names,
                contig_lengths,
            })
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

        // Allocate output buffer for Hi-C pairs mappings
        let mut results = vec![
            RustPairsMapping {
                read_id: 0,
                rid1: -1,
                rid2: -1,
                pos1: 0,
                pos2: 0,
                strand1: 0,
                strand2: 0,
                mapq: 0,
                is_unique: 0,
                num_dups: 0,
            };
            n
        ];

        unsafe {
            eprintln!("[ChromapAligner] Calling mapper_map_batch with {} pairs", n);
            if !r1_seqs.is_empty() {
                let sample_seq: String = r1_seqs[0].chars().take(50).collect();
                let sample_qual: String = r1_quals[0].chars().take(50).collect();
                eprintln!("[ChromapAligner] Sample R1 seq (first 50bp): {}", sample_seq);
                eprintln!("[ChromapAligner] Sample R1 qual (first 50bp): {}", sample_qual);
            } else {
                eprintln!("[ChromapAligner] No sequences to map (empty batch)");
            }
            
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

            eprintln!("[ChromapAligner] mapper_map_batch returned: {}", ret);

            if ret < 0 {
                eprintln!("[ChromapAligner] ERROR: mapper_map_batch failed with code {}", ret);
                return Err(ChromapError::MappingError);
            }

            // Convert to high-level results
            Ok(results
                .into_iter()
                .map(MappingResult::from)
                .collect::<Vec<MappingResult>>())
        }
    }

    /// Get the SAM header
    pub fn get_sam_header(&self) -> sam::Header {
        self.header.clone()
    }

    /* Commented out - not used by hic-tailor integration
    /// Align reads and return SAM records
    ///
    /// # Arguments
    /// * `records` - Vector of FASTQ record pairs (R1, R2)
    ///
    /// # Returns
    /// Vector of SAM RecordBuf objects
    pub fn align_reads(
        &self,
        records: &[(fastq::Record, fastq::Record)],
    ) -> Result<Vec<RecordBuf>, ChromapError> {
        // Implementation commented out - use map_batch instead
        unimplemented!("Use map_batch for hic-tailor integration")
    }

    /// Convert a MappingResult to a SAM RecordBuf
    fn mapping_to_sam_record(
        &self,
        mapping: &MappingResult,
        name: &[u8],
        sequence: &[u8],
        quality: &[u8],
        is_r1: bool,
    ) -> Result<RecordBuf, ChromapError> {
        // Implementation commented out - not needed for hic-tailor
        unimplemented!("Not used in hic-tailor integration")
    }
    */
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
        let mapping = RustPairsMapping {
            read_id: 0,
            rid1: 1,
            rid2: 2,
            pos1: 100,
            pos2: 200,
            strand1: 1,
            strand2: 0,
            mapq: 60,
            is_unique: 1,
            num_dups: 0,
        };

        let result: MappingResult = mapping.into();
        assert_eq!(result.strand1, Strand::Forward);
        assert_eq!(result.strand2, Strand::Reverse);
    }
}

