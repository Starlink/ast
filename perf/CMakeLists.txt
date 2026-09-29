# Performance benchmarks

set(PERF_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/.."
    "${CMAKE_CURRENT_SOURCE_DIR}/../src"
    "${CMAKE_BINARY_DIR}"
    "${CMAKE_BINARY_DIR}/src"
)

set(PERF_LIBS
    ast
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_err-static>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_2.0-static>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_3.2-static>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_5.6-static>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf3d-static>
    m
)

add_executable(bench_lock bench_lock.c)
target_include_directories(bench_lock PRIVATE ${PERF_INCLUDES})
target_link_libraries(bench_lock PRIVATE ${PERF_LIBS})
target_compile_definitions(bench_lock PRIVATE HAVE_CONFIG_H)
# Match the warning/sanitizer/coverage flags applied to the rest of the build
# so the benchmark links correctly against a sanitizer-instrumented libast.
ast_apply_dev_options(bench_lock)
