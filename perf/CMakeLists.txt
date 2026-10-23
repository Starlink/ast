# Performance benchmarks

set(PERF_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/.."
    "${CMAKE_CURRENT_SOURCE_DIR}/../src"
    "${CMAKE_BINARY_DIR}"
    "${CMAKE_BINARY_DIR}/src"
)

set(PERF_LIBS
    ast
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_err>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_2.0>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_3.2>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_5.6>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf3d>
    m
)

add_executable(bench_simd bench_simd.c)
target_include_directories(bench_simd PRIVATE ${PERF_INCLUDES})
target_link_libraries(bench_simd PRIVATE ${PERF_LIBS})
target_compile_definitions(bench_simd PRIVATE HAVE_CONFIG_H)
# Match the warning/sanitizer/coverage flags applied to the rest of the build
# so the benchmark links correctly against a sanitizer-instrumented libast.
ast_apply_dev_options(bench_simd)

