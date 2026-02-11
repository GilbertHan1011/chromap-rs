use std::env;
use std::path::Path;

fn main() {
    let mut build = cc::Build::new();

    // Determine the architecture flags based on environment variables
    let arch = env::var("CARGO_CFG_TARGET_FEATURE").unwrap_or_default();
    println!("cargo:warning=arch: {}", arch);

    let arch_flags = if arch.contains("avx2") {
        vec!["-mavx2"]
    } else if arch.contains("sse4.1") {
        vec!["-msse4.1"]
    } else {
        vec!["-msse4.1"]  // Default to SSE4.1
    };

    for f in arch_flags.into_iter() {
        build.flag(f);
    }

    // Set compiler flags matching chromap's Makefile
    build
        .warnings(false)  // Disable warnings for cleaner output
        .cpp(true)
        .flag("-std=c++11")
        .flag("-O3")
        .flag("-fopenmp")  // Enable OpenMP for parallel processing
        .flag("-fPIC")
        .include("ext/chromap/src");

    // List of chromap source files
    let src_dir = "ext/chromap/src/";
    let files = [
        "sequence_batch.cc",
        "index.cc",
        "minimizer_generator.cc",
        "candidate_processor.cc",
        "alignment.cc",
        "feature_barcode_matrix.cc",
        "ksw.cc",
        "draft_mapping_generator.cc",
        "mapping_generator.cc",
        "mapping_writer.cc",
        "chromap.cc",
    ];

    // Add wrapper and chromap source files to build
    build.file("wrapper.cpp");
    for file in &files {
        build.file(Path::new(src_dir).join(file));
    }

    // Compile to static library
    build.compile("libchromap.a");

    println!("cargo:rustc-link-search=native={}", env::var("OUT_DIR").unwrap());

    // Link required libraries
    println!("cargo:rustc-link-lib=z");      // zlib for compression
    println!("cargo:rustc-link-lib=gomp");   // GNU OpenMP library (use "omp" on macOS)
    println!("cargo:rustc-link-lib=m");      // Math library

    // Rerun if wrapper files change
    println!("cargo:rerun-if-changed=wrapper.cpp");
    println!("cargo:rerun-if-changed=wrapper.h");
}
