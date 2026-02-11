# chromap-rs

Rust FFI bindings for [chromap](https://github.com/haowenz/chromap), a fast alignment and preprocessing tool for chromatin profiles.

This project provides in-memory FFI bindings to chromap. Instead of running chromap as a CLI subprocess with file I/O, it links the C++ code as a static library and provides a Rust API for direct memory-to-memory mapping.

## Project Structure

```
chromap-rs/
├── chromap/              # Upstream chromap source code (git submodule)
├── chromap-sys/          # Low-level FFI bindings
│   ├── wrapper.h         # C API header
│   ├── wrapper.cpp       # C++ wrapper implementation
│   ├── build.rs          # Build script
│   └── src/lib.rs        # Rust FFI declarations
├── chromap-rust/         # High-level Rust API
│   └── src/lib.rs        # Safe Rust interface
└── Cargo.toml            # Workspace configuration
```

## Architecture

### chromap-sys (Low-Level Bindings)

The `chromap-sys` crate provides low-level C FFI bindings:

- **wrapper.h/cpp**: Exposes a C-compatible API that wraps chromap's C++ classes
- **ChromapMapper**: Opaque struct that holds the index, reference, and mapping state
- **RustMapping**: C-compatible struct for mapping results
- **build.rs**: Compiles chromap C++ source files with proper flags (OpenMP, SSE4.1/AVX2)

### chromap-rust (High-Level API)

The `chromap-rust` crate provides a safe, idiomatic Rust interface:

- **ChromapAligner**: Main struct for performing alignments
- **MappingResult**: Rust-native mapping result type
- **ChromapError**: Error handling with proper Rust error types

## Usage

```rust
use chromap_rust::ChromapAligner;

// Create aligner (loads index once, keeps in memory)
let aligner = ChromapAligner::new(
    "reference.index",
    "reference.fa",
    4  // num_threads
)?;

// Map a batch of paired-end reads
let r1_seqs = vec!["ACGTACGT", "TGCATGCA"];
let r2_seqs = vec!["GCTAGCTA", "ATCGATCG"];
let r1_quals = vec!["IIIIIIII", "IIIIIIII"];
let r2_quals = vec!["IIIIIIII", "IIIIIIII"];

let results = aligner.map_batch(&r1_seqs, &r2_seqs, &r1_quals, &r2_quals)?;

for result in results {
    println!("Read {} mapped to ref {} at position {}",
             result.read_id, result.ref_id, result.ref_pos);
}
```

## Building

### Prerequisites

- Rust 1.70+
- C++ compiler with C++11 support
- OpenMP library (libgomp on Linux, libomp on macOS)
- zlib development headers

### Build Commands

```bash
# Build the entire workspace
cargo build --release

# Build only chromap-sys
cargo build -p chromap-sys --release

# Build only chromap-rust
cargo build -p chromap-rust --release

# Run tests
cargo test
```

## Implementation Status

### ✅ Completed

- [x] Project structure setup
- [x] C++ wrapper header and implementation
- [x] Build system configuration
- [x] Low-level FFI bindings (chromap-sys)
- [x] High-level Rust API (chromap-rust)
- [x] Successful compilation
- [x] **Full mapping pipeline implementation**:
  - [x] Creating `SequenceBatch` objects from input sequences
  - [x] Generating minimizers for each read
  - [x] Finding candidates using the index
  - [x] Supplementing candidates with mate information
  - [x] Generating draft mappings
  - [x] Selecting best mappings with MAPQ scores
  - [x] Converting chromap's internal mapping format to `RustMapping`

### ⚠️ TODO

- [ ] **Testing**: Add integration tests with real chromap index and reference data
- [ ] **Validation**: Verify mapping results match chromap CLI output
- [ ] **Documentation**: Add more examples and API documentation
- [ ] **Performance optimization**: Profile and optimize the FFI boundary
- [ ] **Error handling**: Improve error reporting from C++ to Rust (currently returns nullptr on failure)

## Key Differences from bwa-mem2-rust

1. **Chromap Architecture**: Chromap is more tightly coupled with file I/O and has a more complex internal structure than bwa-mem2
2. **Mapping Logic**: Extracting the core mapping loop from chromap requires more refactoring due to its integrated pipeline design
3. **State Management**: Chromap uses more stateful components (caches, barcode processing) that need careful handling

## Performance Considerations

1. **Index Loading**: The index is loaded once during `ChromapAligner::new()` and kept in memory
2. **Thread Safety**: OpenMP handles parallelism internally in C++. Avoid calling `map_batch` from multiple Rust threads simultaneously
3. **Batch Size**: For optimal performance, map reads in batches of 5,000-10,000 pairs
4. **Memory**: The aligner keeps the full reference and index in memory

## Contributing

This is a functional implementation with a complete mapping pipeline. Contributions are welcome, especially for:

- Adding comprehensive integration tests with real data
- Validating mapping accuracy against chromap CLI
- Performance benchmarking and optimization
- Documentation improvements
- Better error handling and reporting

## License

This project follows the same license as chromap (MIT).

## References

- [chromap](https://github.com/haowenz/chromap) - Original chromap implementation
- [bwa-mem2-rust](https://github.com/PacificBiosciences/bwa-mem2-rust) - Reference implementation strategy
