#include "wrapper.h"
#include "ext/chromap/src/index.h"
#include "ext/chromap/src/sequence_batch.h"
#include "ext/chromap/src/minimizer_generator.h"
#include "ext/chromap/src/candidate_processor.h"
#include "ext/chromap/src/mapping_processor.h"
#include "ext/chromap/src/mapping_generator.h"
#include "ext/chromap/src/mapping_parameters.h"
#include "ext/chromap/src/pairs_mapping.h"
#include "ext/chromap/src/draft_mapping_generator.h"
#include "ext/chromap/src/paired_end_mapping_metadata.h"
#include "ext/chromap/src/mapping_metadata.h"
#include <cstring>
#include <vector>
#include <memory>
#include <random>
#include <algorithm>

using namespace chromap;

// Helper function to populate a kseq_t structure from C strings
static void populate_kseq(SequenceBatch::kseq_t* seq, uint32_t id, const char* name,
                         const char* sequence, const char* qual) {
    // Set ID
    seq->id = id;

    // Set name
    size_t name_len = strlen(name);
    seq->name.l = name_len;
    seq->name.m = name_len + 1;
    seq->name.s = (char*)malloc(seq->name.m);
    memcpy(seq->name.s, name, name_len);
    seq->name.s[name_len] = '\0';

    // Set sequence
    size_t seq_len = strlen(sequence);
    seq->seq.l = seq_len;
    seq->seq.m = seq_len + 1;
    seq->seq.s = (char*)malloc(seq->seq.m);
    memcpy(seq->seq.s, sequence, seq_len);
    seq->seq.s[seq_len] = '\0';

    // Set quality
    size_t qual_len = strlen(qual);
    seq->qual.l = qual_len;
    seq->qual.m = qual_len + 1;
    seq->qual.s = (char*)malloc(seq->qual.m);
    memcpy(seq->qual.s, qual, qual_len);
    seq->qual.s[qual_len] = '\0';

    // Set comment to empty
    seq->comment.l = 0;
    seq->comment.m = 1;
    seq->comment.s = (char*)malloc(1);
    seq->comment.s[0] = '\0';

    seq->last_char = 0;
    seq->f = nullptr;
}

// Internal structure to hold the mapper state
struct ChromapMapper {
    Index* index;
    SequenceBatch* reference;
    MappingParameters params;
    int num_threads;
    int kmer_size;
    int window_size;
    uint32_t num_reference_sequences;

    // Pre-allocated processors
    MinimizerGenerator* minimizer_generator;
    CandidateProcessor* candidate_processor;
    DraftMappingGenerator* draft_mapping_generator;
    MappingGenerator<PairsMapping>* mapping_generator;

    ChromapMapper() : index(nullptr), reference(nullptr),
                      minimizer_generator(nullptr), candidate_processor(nullptr),
                      draft_mapping_generator(nullptr), mapping_generator(nullptr) {}

    ~ChromapMapper() {
        delete index;
        delete reference;
        delete minimizer_generator;
        delete candidate_processor;
        delete draft_mapping_generator;
        delete mapping_generator;
    }
};

extern "C" {

ChromapMapper* mapper_new(const char* index_path, const char* ref_path, int num_threads) {
    ChromapMapper* mapper = new ChromapMapper();

    try {
        // Set up mapping parameters with defaults
        mapper->params.index_file_path = std::string(index_path);
        mapper->params.reference_file_path = std::string(ref_path);
        mapper->params.num_threads = num_threads;
        mapper->num_threads = num_threads;

        // Apply Hi-C preset parameters (equivalent to chromap --preset hic)
        mapper->params.error_threshold = 4;              // Stricter than default (was 8)
        mapper->params.mapq_threshold = 1;               // Very permissive (default is 30)
        mapper->params.split_alignment = true;            // Allow chimeric alignments at ligation junctions
        mapper->params.low_memory_mode = true;           // Hi-C preset uses low memory mode
        mapper->params.mapping_output_format = MAPPINGFORMAT_PAIRS;  // Hi-C pairs format
        
        // Keep standard minimizer/seed settings
        mapper->params.min_num_seeds_required_for_mapping = 2;
        mapper->params.max_seed_frequencies.resize(2);
        mapper->params.max_seed_frequencies[0] = 500;
        mapper->params.max_seed_frequencies[1] = 1000;
        mapper->params.max_num_best_mappings = 1;
        mapper->params.max_insert_size = 1000;           // Not strictly enforced in Hi-C mode
        mapper->params.is_bulk_data = true;
        
        // Initialize string fields that might be accessed during mapping
        mapper->params.read_format = "";
        mapper->params.mapping_output_file_path = "";
        mapper->params.barcode_whitelist_file_path = "";
        mapper->params.custom_rid_order_file_path = "";
        mapper->params.pairs_flipping_custom_rid_order_file_path = "";
        mapper->params.barcode_translate_table_file_path = "";
        mapper->params.summary_metadata_file_path = "";
        mapper->params.matrix_output_prefix = "";

        // Load reference
        mapper->reference = new SequenceBatch();
        mapper->reference->InitializeLoading(mapper->params.reference_file_path);
        mapper->reference->LoadAllSequences();
        mapper->num_reference_sequences = mapper->reference->GetNumSequences();

        // Load index
        mapper->index = new Index(mapper->params.index_file_path);
        mapper->index->Load();
        mapper->kmer_size = mapper->index->GetKmerSize();
        mapper->window_size = mapper->index->GetWindowSize();

        // Initialize processors
        std::cerr << "[mapper_new] Initializing MinimizerGenerator..." << std::endl;
        mapper->minimizer_generator = new MinimizerGenerator(mapper->kmer_size, mapper->window_size);
        std::cerr << "[mapper_new] Initializing CandidateProcessor..." << std::endl;
        mapper->candidate_processor = new CandidateProcessor(
            mapper->params.min_num_seeds_required_for_mapping,
            mapper->params.max_seed_frequencies);
        std::cerr << "[mapper_new] Initializing DraftMappingGenerator..." << std::endl;
        mapper->draft_mapping_generator = new DraftMappingGenerator(mapper->params);

        std::cerr << "[mapper_new] Initializing MappingGenerator (Hi-C PairsMapping mode)..." << std::endl;
        // For Hi-C pairs format, we need to provide chromosome rank ordering
        // Empty rank vector means use natural ordering (chr1, chr2, ...)
        std::vector<int> chr_ranks(mapper->num_reference_sequences);
        for (uint32_t i = 0; i < mapper->num_reference_sequences; ++i) {
            chr_ranks[i] = i;  // Natural ordering by index
        }
        mapper->mapping_generator = new MappingGenerator<PairsMapping>(mapper->params, chr_ranks);
        std::cerr << "[mapper_new] All processors initialized successfully" << std::endl;

        std::cerr << "[mapper_new] Mapper initialization completed successfully" << std::endl;
        return mapper;
    } catch (const std::exception& e) {
        std::cerr << "[mapper_new] ERROR: Exception during initialization: " << e.what() << std::endl;
        delete mapper;
        return nullptr;
    } catch (...) {
        std::cerr << "[mapper_new] ERROR: Unknown exception during initialization" << std::endl;
        delete mapper;
        return nullptr;
    }
}

void mapper_free(ChromapMapper* mapper) {
    if (mapper) {
        delete mapper;
    }
}

int mapper_map_batch(ChromapMapper* mapper,
                     const char** seqs1, const char** seqs2,
                     const char** quals1, const char** quals2,
                     int n_reads,
                     RustPairsMapping* out_buffer,
                     int buffer_size) {
    if (!mapper || !seqs1 || !seqs2 || n_reads <= 0 || !out_buffer || buffer_size < n_reads) {
        std::cerr << "[C++ Wrapper] Invalid input parameters" << std::endl;
        return -1;
    }

    std::cerr << "[C++ Wrapper] Starting mapper_map_batch with " << n_reads << " reads" << std::endl;

    try {
        // Create SequenceBatch objects for the reads
        std::cerr << "[C++ Wrapper] Creating SequenceBatch objects..." << std::endl;
        SequenceEffectiveRange full_range;  // Full sequence range
        SequenceBatch read_batch1(n_reads, full_range);
        SequenceBatch read_batch2(n_reads, full_range);
        std::cerr << "[C++ Wrapper] SequenceBatch objects created successfully" << std::endl;

        // Populate the sequence batches
        std::cerr << "[C++ Wrapper] Populating sequence batches..." << std::endl;
        auto& batch1_seqs = read_batch1.GetSequenceBatch();
        auto& batch2_seqs = read_batch2.GetSequenceBatch();

        for (int i = 0; i < n_reads; ++i) {
            char read_name[32];
            snprintf(read_name, sizeof(read_name), "read_%d", i);

            populate_kseq(batch1_seqs[i], i * 2, read_name, seqs1[i], quals1[i]);
            populate_kseq(batch2_seqs[i], i * 2 + 1, read_name, seqs2[i], quals2[i]);
            
            if (i == 0) {
                std::cerr << "[C++ Wrapper] First pair - R1 len: " << strlen(seqs1[i]) 
                          << ", R2 len: " << strlen(seqs2[i]) << std::endl;
            }
        }
        std::cerr << "[C++ Wrapper] Populated " << n_reads << " pairs successfully" << std::endl;

        // Prepare negative sequences (reverse complement)
        std::cerr << "[C++ Wrapper] Preparing negative sequences..." << std::endl;
        for (int i = 0; i < n_reads; ++i) {
            read_batch1.PrepareNegativeSequenceAt(i);
            read_batch2.PrepareNegativeSequenceAt(i);
        }
        std::cerr << "[C++ Wrapper] Negative sequences prepared successfully" << std::endl;

        // Process each read pair
        std::cerr << "[C++ Wrapper] Starting main mapping loop for " << n_reads << " pairs..." << std::endl;
        std::mt19937 generator(11);
        std::vector<int> best_mapping_indices(mapper->params.max_num_best_mappings);

        for (int pair_idx = 0; pair_idx < n_reads; ++pair_idx) {
            if (pair_idx == 0 || pair_idx % 5000 == 0) {
                std::cerr << "[C++ Wrapper] Processing pair " << pair_idx << "/" << n_reads << std::endl;
            }
            
            // Skip if reads are too short
            if (read_batch1.GetSequenceLengthAt(pair_idx) < 20 ||
                read_batch2.GetSequenceLengthAt(pair_idx) < 20) {
                out_buffer[pair_idx].read_id = pair_idx;
                out_buffer[pair_idx].rid1 = -1;  // Unmapped
                out_buffer[pair_idx].rid2 = -1;  // Unmapped
                out_buffer[pair_idx].pos1 = 0;
                out_buffer[pair_idx].pos2 = 0;
                out_buffer[pair_idx].strand1 = 0;
                out_buffer[pair_idx].strand2 = 0;
                out_buffer[pair_idx].mapq = 0;
                out_buffer[pair_idx].is_unique = 0;
                out_buffer[pair_idx].num_dups = 0;
                continue;
            }

          try {
            // Create metadata for this read pair
            if (pair_idx == 0) std::cerr << "[C++ Wrapper] Creating PairedEndMappingMetadata..." << std::endl;
            PairedEndMappingMetadata paired_metadata;
            if (pair_idx == 0) std::cerr << "[C++ Wrapper] Calling PreparedForMappingNextReadPair..." << std::endl;
            paired_metadata.PreparedForMappingNextReadPair(
                mapper->params.max_seed_frequencies[0]);

            // Generate minimizers for both reads
            if (pair_idx == 0) std::cerr << "[C++ Wrapper] Generating minimizers for R1..." << std::endl;
            mapper->minimizer_generator->GenerateMinimizers(
                read_batch1, pair_idx,
                paired_metadata.mapping_metadata1_.minimizers_);
            if (pair_idx == 0) std::cerr << "[C++ Wrapper] Generating minimizers for R2..." << std::endl;
            mapper->minimizer_generator->GenerateMinimizers(
                read_batch2, pair_idx,
                paired_metadata.mapping_metadata2_.minimizers_);

            // Check if both ends have minimizers
            if (!paired_metadata.BothEndsHaveMinimizers()) {
                out_buffer[pair_idx].read_id = pair_idx;
                out_buffer[pair_idx].rid1 = -1;  // Unmapped
                out_buffer[pair_idx].rid2 = -1;  // Unmapped
                out_buffer[pair_idx].pos1 = 0;
                out_buffer[pair_idx].pos2 = 0;
                out_buffer[pair_idx].strand1 = 0;
                out_buffer[pair_idx].strand2 = 0;
                out_buffer[pair_idx].mapq = 0;
                out_buffer[pair_idx].is_unique = 0;
                out_buffer[pair_idx].num_dups = 0;
                continue;
            }

            // Generate candidates for both reads
            if (pair_idx == 0) std::cerr << "[C++ Wrapper] Generating candidates for R1..." << std::endl;
            mapper->candidate_processor->GenerateCandidates(
                mapper->params.error_threshold,
                *(mapper->index),
                paired_metadata.mapping_metadata1_);
            if (pair_idx == 0) std::cerr << "[C++ Wrapper] Generating candidates for R2..." << std::endl;
            mapper->candidate_processor->GenerateCandidates(
                mapper->params.error_threshold,
                *(mapper->index),
                paired_metadata.mapping_metadata2_);

            size_t num_candidates1 = paired_metadata.mapping_metadata1_.GetNumCandidates();
            size_t num_candidates2 = paired_metadata.mapping_metadata2_.GetNumCandidates();

            // If we have candidates for both reads, process them
            if (num_candidates1 > 0 && num_candidates2 > 0) {
                // Supplement candidates with mate information
                if (!mapper->params.split_alignment) {
                    if (pair_idx == 0) std::cerr << "[C++ Wrapper] Supplementing candidates..." << std::endl;
                    mapper->candidate_processor->SupplementCandidates(
                        mapper->params.error_threshold,
                        2 * mapper->params.max_insert_size,
                        *(mapper->index),
                        paired_metadata);

                    // Move candidates to buffer and reduce based on insert size
                    if (pair_idx == 0) std::cerr << "[C++ Wrapper] Moving candidates to buffer..." << std::endl;
                    paired_metadata.MoveCandidiatesToBuffer();
                    if (pair_idx == 0) std::cerr << "[C++ Wrapper] Reducing candidates..." << std::endl;
                    mapper->candidate_processor->ReduceCandidatesForPairedEndRead(
                        mapper->params.max_insert_size,
                        paired_metadata);

                    if (pair_idx == 0) std::cerr << "[C++ Wrapper] Getting candidate counts..." << std::endl;
                    num_candidates1 = paired_metadata.mapping_metadata1_.GetNumCandidates();
                    if (pair_idx == 0) std::cerr << "[C++ Wrapper] R1 candidates: " << num_candidates1 << std::endl;
                    num_candidates2 = paired_metadata.mapping_metadata2_.GetNumCandidates();
                    if (pair_idx == 0) std::cerr << "[C++ Wrapper] R2 candidates: " << num_candidates2 << std::endl;
                }

                // Generate draft mappings if we still have candidates
                if (pair_idx == 0) std::cerr << "[C++ Wrapper] Checking if we have candidates for draft mapping..." << std::endl;
                if (num_candidates1 > 0 && num_candidates2 > 0) {
                    if (pair_idx == 0) std::cerr << "[C++ Wrapper] Generating draft mappings for R1..." << std::endl;
                    mapper->draft_mapping_generator->GenerateDraftMappings(
                        read_batch1, pair_idx, *(mapper->reference),
                        paired_metadata.mapping_metadata1_);
                    if (pair_idx == 0) std::cerr << "[C++ Wrapper] Generating draft mappings for R2..." << std::endl;
                    mapper->draft_mapping_generator->GenerateDraftMappings(
                        read_batch2, pair_idx, *(mapper->reference),
                        paired_metadata.mapping_metadata2_);

                    size_t num_draft1 = paired_metadata.mapping_metadata1_.GetNumDraftMappings();
                    size_t num_draft2 = paired_metadata.mapping_metadata2_.GetNumDraftMappings();

                    // Generate best mappings if we have draft mappings
                    if (num_draft1 > 0 && num_draft2 > 0) {
                        // Sort mappings by position for paired-end processing
                        if (!mapper->params.split_alignment) {
                            if (pair_idx == 0) std::cerr << "[C++ Wrapper] Sorting mappings by positions..." << std::endl;
                            paired_metadata.SortMappingsByPositions();
                        }

                        // Create temporary storage for Hi-C pairs mappings
                        if (pair_idx == 0) std::cerr << "[C++ Wrapper] Creating storage for " << mapper->num_reference_sequences << " references..." << std::endl;
                        std::vector<std::vector<PairsMapping>> mappings_on_refs(
                            mapper->num_reference_sequences);

                        // Dummy barcode batch (not used for bulk data)
                        if (pair_idx == 0) std::cerr << "[C++ Wrapper] Creating dummy barcode batch..." << std::endl;
                        SequenceBatch dummy_barcode(0, full_range);

                        // Generate best mappings for Hi-C pairs
                        if (pair_idx == 0) std::cerr << "[C++ Wrapper] Generating best Hi-C pairs mappings..." << std::endl;
                        mapper->mapping_generator->GenerateBestMappingsForPairedEndRead(
                            pair_idx, read_batch1, read_batch2, dummy_barcode,
                            *(mapper->reference), best_mapping_indices, generator,
                            -1,  // force_mapq
                            paired_metadata, mappings_on_refs);

                        // Extract the first pairs mapping if available
                        // Note: PairsMapping contains BOTH ends' information in one record
                        bool found_mapping = false;
                        for (const auto& ref_mappings : mappings_on_refs) {
                            if (!ref_mappings.empty()) {
                                const PairsMapping& mapping = ref_mappings[0];
                                out_buffer[pair_idx].read_id = pair_idx;
                                out_buffer[pair_idx].rid1 = mapping.rid1_;
                                out_buffer[pair_idx].rid2 = mapping.rid2_;
                                out_buffer[pair_idx].pos1 = mapping.pos1_;
                                out_buffer[pair_idx].pos2 = mapping.pos2_;
                                out_buffer[pair_idx].strand1 = mapping.strand1_;
                                out_buffer[pair_idx].strand2 = mapping.strand2_;
                                out_buffer[pair_idx].mapq = mapping.mapq_;
                                out_buffer[pair_idx].is_unique = mapping.is_unique_;
                                out_buffer[pair_idx].num_dups = mapping.num_dups_;
                                found_mapping = true;
                                break;
                            }
                        }

                        if (!found_mapping) {
                            // No mapping found
                            out_buffer[pair_idx].read_id = pair_idx;
                            out_buffer[pair_idx].rid1 = -1;
                            out_buffer[pair_idx].rid2 = -1;
                            out_buffer[pair_idx].pos1 = 0;
                            out_buffer[pair_idx].pos2 = 0;
                            out_buffer[pair_idx].strand1 = 0;
                            out_buffer[pair_idx].strand2 = 0;
                            out_buffer[pair_idx].mapq = 0;
                            out_buffer[pair_idx].is_unique = 0;
                            out_buffer[pair_idx].num_dups = 0;
                        }
                    } else {
                        // No draft mappings
                        out_buffer[pair_idx].read_id = pair_idx;
                        out_buffer[pair_idx].rid1 = -1;
                        out_buffer[pair_idx].rid2 = -1;
                        out_buffer[pair_idx].pos1 = 0;
                        out_buffer[pair_idx].pos2 = 0;
                        out_buffer[pair_idx].strand1 = 0;
                        out_buffer[pair_idx].strand2 = 0;
                        out_buffer[pair_idx].mapq = 0;
                        out_buffer[pair_idx].is_unique = 0;
                        out_buffer[pair_idx].num_dups = 0;
                    }
                } else {
                    // No candidates after filtering
                    out_buffer[pair_idx].read_id = pair_idx;
                    out_buffer[pair_idx].rid1 = -1;
                    out_buffer[pair_idx].rid2 = -1;
                    out_buffer[pair_idx].pos1 = 0;
                    out_buffer[pair_idx].pos2 = 0;
                    out_buffer[pair_idx].strand1 = 0;
                    out_buffer[pair_idx].strand2 = 0;
                    out_buffer[pair_idx].mapq = 0;
                    out_buffer[pair_idx].is_unique = 0;
                    out_buffer[pair_idx].num_dups = 0;
                }
            } else {
                // No candidates for one or both reads
                out_buffer[pair_idx].read_id = pair_idx;
                out_buffer[pair_idx].rid1 = -1;
                out_buffer[pair_idx].rid2 = -1;
                out_buffer[pair_idx].pos1 = 0;
                out_buffer[pair_idx].pos2 = 0;
                out_buffer[pair_idx].strand1 = 0;
                out_buffer[pair_idx].strand2 = 0;
                out_buffer[pair_idx].mapq = 0;
                out_buffer[pair_idx].is_unique = 0;
                out_buffer[pair_idx].num_dups = 0;
            }
          } catch (const std::exception& e) {
              std::cerr << "[C++ Wrapper] ERROR at pair_idx=" << pair_idx << ": " << e.what() << std::endl;
              // Log read lengths for the failing pair
              std::cerr << "[C++ Wrapper] Failing pair R1 len: " << read_batch1.GetSequenceLengthAt(pair_idx)
                        << ", R2 len: " << read_batch2.GetSequenceLengthAt(pair_idx) << std::endl;
              return -1;
          } catch (...) {
              std::cerr << "[C++ Wrapper] ERROR at pair_idx=" << pair_idx << ": unknown exception" << std::endl;
              return -1;
          }
        }

        std::cerr << "[C++ Wrapper] Mapping completed successfully, returning " << n_reads << std::endl;
        return n_reads;
    } catch (const std::exception& e) {
        std::cerr << "[C++ Wrapper] ERROR: Exception caught: " << e.what() << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "[C++ Wrapper] ERROR: Unknown exception caught" << std::endl;
        return -1;
    }
}

// Get the number of reference sequences
int mapper_get_num_references(ChromapMapper* mapper) {
    if (!mapper || !mapper->reference) {
        return -1;
    }
    return mapper->num_reference_sequences;
}

// Get the name of a reference sequence by index
const char* mapper_get_reference_name(ChromapMapper* mapper, int index) {
    if (!mapper || !mapper->reference || index < 0 || index >= (int)mapper->num_reference_sequences) {
        return nullptr;
    }
    return mapper->reference->GetSequenceNameAt(index);
}

// Get the length of a reference sequence by index
int mapper_get_reference_length(ChromapMapper* mapper, int index) {
    if (!mapper || !mapper->reference || index < 0 || index >= (int)mapper->num_reference_sequences) {
        return -1;
    }
    return mapper->reference->GetSequenceLengthAt(index);
}

} // extern "C"
