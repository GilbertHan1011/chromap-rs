#include "wrapper.h"
#include "../chromap/src/index.h"
#include "../chromap/src/sequence_batch.h"
#include "../chromap/src/minimizer_generator.h"
#include "../chromap/src/candidate_processor.h"
#include "../chromap/src/mapping_processor.h"
#include "../chromap/src/mapping_generator.h"
#include "../chromap/src/mapping_parameters.h"
#include "../chromap/src/sam_mapping.h"
#include "../chromap/src/draft_mapping_generator.h"
#include "../chromap/src/paired_end_mapping_metadata.h"
#include "../chromap/src/mapping_metadata.h"
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
    MappingGenerator<SAMMapping>* mapping_generator;

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

        // Set default mapping parameters
        mapper->params.error_threshold = 8;
        mapper->params.min_num_seeds_required_for_mapping = 2;
        mapper->params.max_seed_frequencies.resize(2);
        mapper->params.max_seed_frequencies[0] = 500;
        mapper->params.max_seed_frequencies[1] = 1000;
        mapper->params.max_num_best_mappings = 1;
        mapper->params.max_insert_size = 1000;
        mapper->params.split_alignment = false;

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
        mapper->minimizer_generator = new MinimizerGenerator(mapper->kmer_size, mapper->window_size);
        mapper->candidate_processor = new CandidateProcessor(
            mapper->params.min_num_seeds_required_for_mapping,
            mapper->params.max_seed_frequencies);
        mapper->draft_mapping_generator = new DraftMappingGenerator(mapper->params);

        std::vector<int> empty_ranks;  // Empty custom rid ranks
        mapper->mapping_generator = new MappingGenerator<SAMMapping>(mapper->params, empty_ranks);

        return mapper;
    } catch (...) {
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
                     RustMapping* out_buffer,
                     int buffer_size) {
    if (!mapper || !seqs1 || !seqs2 || n_reads <= 0 || !out_buffer || buffer_size < n_reads) {
        return -1;
    }

    try {
        // Create SequenceBatch objects for the reads
        SequenceEffectiveRange full_range;  // Full sequence range
        SequenceBatch read_batch1(n_reads, full_range);
        SequenceBatch read_batch2(n_reads, full_range);

        // Populate the sequence batches
        auto& batch1_seqs = read_batch1.GetSequenceBatch();
        auto& batch2_seqs = read_batch2.GetSequenceBatch();

        for (int i = 0; i < n_reads; ++i) {
            char read_name[32];
            snprintf(read_name, sizeof(read_name), "read_%d", i);

            populate_kseq(batch1_seqs[i], i * 2, read_name, seqs1[i], quals1[i]);
            populate_kseq(batch2_seqs[i], i * 2 + 1, read_name, seqs2[i], quals2[i]);
        }

        // Prepare negative sequences (reverse complement)
        for (int i = 0; i < n_reads; ++i) {
            read_batch1.PrepareNegativeSequenceAt(i);
            read_batch2.PrepareNegativeSequenceAt(i);
        }

        // Process each read pair
        std::vector<SAMMapping> all_mappings;
        std::mt19937 generator(11);
        std::vector<int> best_mapping_indices(mapper->params.max_num_best_mappings);

        for (int pair_idx = 0; pair_idx < n_reads; ++pair_idx) {
            // Skip if reads are too short
            if (read_batch1.GetSequenceLengthAt(pair_idx) < 20 ||
                read_batch2.GetSequenceLengthAt(pair_idx) < 20) {
                out_buffer[pair_idx].read_id = pair_idx;
                out_buffer[pair_idx].ref_id = 0xFFFFFFFF;  // Unmapped
                out_buffer[pair_idx].ref_pos = 0;
                out_buffer[pair_idx].mapq = 0;
                out_buffer[pair_idx].strand = 0;
                out_buffer[pair_idx].num_dups = 0;
                out_buffer[pair_idx].is_unique = 0;
                out_buffer[pair_idx].is_primary = 1;
                continue;
            }

            // Create metadata for this read pair
            PairedEndMappingMetadata paired_metadata;
            paired_metadata.PreparedForMappingNextReadPair(
                mapper->params.max_seed_frequencies[0]);

            // Generate minimizers for both reads
            mapper->minimizer_generator->GenerateMinimizers(
                read_batch1, pair_idx,
                paired_metadata.mapping_metadata1_.minimizers_);
            mapper->minimizer_generator->GenerateMinimizers(
                read_batch2, pair_idx,
                paired_metadata.mapping_metadata2_.minimizers_);

            // Check if both ends have minimizers
            if (!paired_metadata.BothEndsHaveMinimizers()) {
                out_buffer[pair_idx].read_id = pair_idx;
                out_buffer[pair_idx].ref_id = 0xFFFFFFFF;  // Unmapped
                out_buffer[pair_idx].ref_pos = 0;
                out_buffer[pair_idx].mapq = 0;
                out_buffer[pair_idx].strand = 0;
                out_buffer[pair_idx].num_dups = 0;
                out_buffer[pair_idx].is_unique = 0;
                out_buffer[pair_idx].is_primary = 1;
                continue;
            }

            // Generate candidates for both reads
            mapper->candidate_processor->GenerateCandidates(
                mapper->params.error_threshold,
                *(mapper->index),
                paired_metadata.mapping_metadata1_);
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
                    mapper->candidate_processor->SupplementCandidates(
                        mapper->params.error_threshold,
                        2 * mapper->params.max_insert_size,
                        *(mapper->index),
                        paired_metadata);

                    // Move candidates to buffer and reduce based on insert size
                    paired_metadata.MoveCandidiatesToBuffer();
                    mapper->candidate_processor->ReduceCandidatesForPairedEndRead(
                        mapper->params.max_insert_size,
                        paired_metadata);

                    num_candidates1 = paired_metadata.mapping_metadata1_.GetNumCandidates();
                    num_candidates2 = paired_metadata.mapping_metadata2_.GetNumCandidates();
                }

                // Generate draft mappings if we still have candidates
                if (num_candidates1 > 0 && num_candidates2 > 0) {
                    mapper->draft_mapping_generator->GenerateDraftMappings(
                        read_batch1, pair_idx, *(mapper->reference),
                        paired_metadata.mapping_metadata1_);
                    mapper->draft_mapping_generator->GenerateDraftMappings(
                        read_batch2, pair_idx, *(mapper->reference),
                        paired_metadata.mapping_metadata2_);

                    size_t num_draft1 = paired_metadata.mapping_metadata1_.GetNumDraftMappings();
                    size_t num_draft2 = paired_metadata.mapping_metadata2_.GetNumDraftMappings();

                    // Generate best mappings if we have draft mappings
                    if (num_draft1 > 0 && num_draft2 > 0) {
                        // Sort mappings by position for paired-end processing
                        if (!mapper->params.split_alignment) {
                            paired_metadata.SortMappingsByPositions();
                        }

                        // Create temporary storage for mappings
                        std::vector<std::vector<SAMMapping>> mappings_on_refs(
                            mapper->num_reference_sequences);

                        // Dummy barcode batch (not used for bulk data)
                        SequenceBatch dummy_barcode(0, full_range);

                        // Generate best mappings
                        mapper->mapping_generator->GenerateBestMappingsForPairedEndRead(
                            pair_idx, read_batch1, read_batch2, dummy_barcode,
                            *(mapper->reference), best_mapping_indices, generator,
                            -1,  // force_mapq
                            paired_metadata, mappings_on_refs);

                        // Extract the first mapping if available
                        bool found_mapping = false;
                        for (const auto& ref_mappings : mappings_on_refs) {
                            if (!ref_mappings.empty()) {
                                const SAMMapping& mapping = ref_mappings[0];
                                out_buffer[pair_idx].read_id = pair_idx;
                                out_buffer[pair_idx].ref_id = mapping.rid_;
                                out_buffer[pair_idx].ref_pos = mapping.pos_;
                                out_buffer[pair_idx].mapq = mapping.mapq_;
                                out_buffer[pair_idx].strand = (mapping.is_rev_ ? 1 : 0);
                                out_buffer[pair_idx].num_dups = 0;
                                out_buffer[pair_idx].is_unique = (paired_metadata.GetNumBestMappings() == 1 ? 1 : 0);
                                out_buffer[pair_idx].is_primary = 1;
                                found_mapping = true;
                                break;
                            }
                        }

                        if (!found_mapping) {
                            // No mapping found
                            out_buffer[pair_idx].read_id = pair_idx;
                            out_buffer[pair_idx].ref_id = 0xFFFFFFFF;
                            out_buffer[pair_idx].ref_pos = 0;
                            out_buffer[pair_idx].mapq = 0;
                            out_buffer[pair_idx].strand = 0;
                            out_buffer[pair_idx].num_dups = 0;
                            out_buffer[pair_idx].is_unique = 0;
                            out_buffer[pair_idx].is_primary = 1;
                        }
                    } else {
                        // No draft mappings
                        out_buffer[pair_idx].read_id = pair_idx;
                        out_buffer[pair_idx].ref_id = 0xFFFFFFFF;
                        out_buffer[pair_idx].ref_pos = 0;
                        out_buffer[pair_idx].mapq = 0;
                        out_buffer[pair_idx].strand = 0;
                        out_buffer[pair_idx].num_dups = 0;
                        out_buffer[pair_idx].is_unique = 0;
                        out_buffer[pair_idx].is_primary = 1;
                    }
                } else {
                    // No candidates after filtering
                    out_buffer[pair_idx].read_id = pair_idx;
                    out_buffer[pair_idx].ref_id = 0xFFFFFFFF;
                    out_buffer[pair_idx].ref_pos = 0;
                    out_buffer[pair_idx].mapq = 0;
                    out_buffer[pair_idx].strand = 0;
                    out_buffer[pair_idx].num_dups = 0;
                    out_buffer[pair_idx].is_unique = 0;
                    out_buffer[pair_idx].is_primary = 1;
                }
            } else {
                // No candidates for one or both reads
                out_buffer[pair_idx].read_id = pair_idx;
                out_buffer[pair_idx].ref_id = 0xFFFFFFFF;
                out_buffer[pair_idx].ref_pos = 0;
                out_buffer[pair_idx].mapq = 0;
                out_buffer[pair_idx].strand = 0;
                out_buffer[pair_idx].num_dups = 0;
                out_buffer[pair_idx].is_unique = 0;
                out_buffer[pair_idx].is_primary = 1;
            }
        }

        return n_reads;
    } catch (...) {
        return -1;
    }
}

} // extern "C"
