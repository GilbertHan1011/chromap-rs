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
