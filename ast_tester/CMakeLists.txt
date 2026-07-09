# Test programs from the ast_tester directory.
#
# These are the existing C test programs, plus Fortran tests that have
# been converted to C.  All tests use CTest and do not require any
# Starlink infrastructure.
include(CheckCSourceRuns)

# Ensure test executables find the locally-built libast before any
# system-installed version (e.g. from a conda environment whose lib
# directory may be injected into RPATH by the compiler toolchain via
# $LDFLAGS, and the --disable-new-dtags flag forces RPATH which is
# searched before LD_LIBRARY_PATH).  Prepending -rpath for the build
# directory puts it first in the RPATH list.
set(CMAKE_EXE_LINKER_FLAGS
    "-Wl,-rpath,${CMAKE_BINARY_DIR} ${CMAKE_EXE_LINKER_FLAGS}")

# Common include directories for all tests
set(TEST_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/.."          # top-level (pal2ast.h etc.)
    "${CMAKE_CURRENT_SOURCE_DIR}/../src"      # internal headers
    "${CMAKE_CURRENT_SOURCE_DIR}/../wcslib"   # wcslib headers
    "${CMAKE_BINARY_DIR}"                     # generated ast.h, ast_err.h, config.h
    "${CMAKE_BINARY_DIR}/src"                 # generated version.h, object.h
    "${CMAKE_SOURCE_DIR}/cmake"               # sae_par.h stub
    "${CMAKE_SOURCE_DIR}/cmake/star"          # mers.h stub
)

# Common link libraries for all tests.
# The satellite libraries (ast_err, ast_grf_*) are static archives whose
# symbols are referenced by libast at runtime (via -undefined dynamic_lookup
# on macOS).  We must force them into the executable with whole-archive
# linkage, otherwise the linker drops the unreferenced archive members and
# libast crashes on null function pointers at runtime.
set(TEST_LIBS
    ast
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_err>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_2.0>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_3.2>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_5.6>
    $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf3d>
    m
)

# Helper function to add a test executable
function(ast_add_test name)
    cmake_parse_arguments(T "INTERNAL_HEADERS;NO_DEFAULT_GRF_LIBS" "" "SOURCES;EXTRA_LIBS;TEST_ARGS" ${ARGN})
    if(NOT T_SOURCES)
        set(T_SOURCES ${name}.c)
    endif()
    add_executable(${name} ${T_SOURCES})
    target_include_directories(${name} PRIVATE ${TEST_INCLUDES})
    
    if(T_NO_DEFAULT_GRF_LIBS)
        set(TARGET_LIBS ast $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_err> ${T_EXTRA_LIBS} m)
    else()
        set(TARGET_LIBS ${TEST_LIBS} ${T_EXTRA_LIBS})
    endif()
    
    target_link_libraries(${name} PRIVATE ${TARGET_LIBS})
    target_compile_definitions(${name} PRIVATE HAVE_CONFIG_H)
    ast_apply_dev_options(${name})
    add_test(NAME ${name} COMMAND ${name} ${T_TEST_ARGS}
             WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

    set(test_env "")
    # AST's memory system intentionally caches freed blocks for reuse,
    # so LeakSanitizer reports false positives. Disable leak detection.
    if(AST_ENABLE_SANITIZERS)
        list(APPEND test_env "ASAN_OPTIONS=detect_leaks=0")
    endif()

    # Tests locate their input fixture files relative to $srcdir, so point it
    # at the source directory where those fixtures live.
    list(APPEND test_env "srcdir=${CMAKE_CURRENT_SOURCE_DIR}")

    set_tests_properties(${name} PROPERTIES
        ENVIRONMENT "${test_env}")
endfunction()

# Read a pipe-delimited text file into a list of encoded records.
# Skips blank lines and lines starting with #.
# Each element of <out_var> is a pipe-separated string of stripped fields.
function(ast_read_test_data file out_var nfields)
    file(STRINGS "${file}" _raw_lines)
    set(_result)
    foreach(_line IN LISTS _raw_lines)
        string(STRIP "${_line}" _line)
        if(_line STREQUAL "" OR _line MATCHES "^#")
            continue()
        endif()
        string(REPLACE "|" ";" _fields "${_line}")
        list(LENGTH _fields _nf)
        if(NOT _nf EQUAL ${nfields})
            message(FATAL_ERROR
                "${file}: expected ${nfields} fields, got ${_nf} in: ${_line}")
        endif()
        set(_stripped)
        foreach(_f IN LISTS _fields)
            string(STRIP "${_f}" _f)
            list(APPEND _stripped "${_f}")
        endforeach()
        string(REPLACE ";" "|" _encoded "${_stripped}")
        list(APPEND _result "${_encoded}")
    endforeach()
    set(${out_var} "${_result}" PARENT_SCOPE)
endfunction()

# --- Tests using public API only ---
ast_add_test(testerror)
ast_add_test(testobject)
ast_add_test(testconvert)

# --- Transform output oracle ---
ast_add_test(test_oracle_util SOURCES test_oracle_util.c transform_oracle_util.c)
set_target_properties(test_oracle_util PROPERTIES C_STANDARD 11)
# Append so we keep the ASAN_OPTIONS=detect_leaks=0 that ast_add_test sets
# under sanitizer builds, rather than replacing it.
set_property(TEST test_oracle_util APPEND PROPERTY
    ENVIRONMENT "ORACLE_TEST_ROOT=${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(gen_transform_oracle gen_transform_oracle.c transform_oracle_util.c)
target_include_directories(gen_transform_oracle PRIVATE ${TEST_INCLUDES})
target_link_libraries(gen_transform_oracle PRIVATE ${TEST_LIBS})
target_compile_definitions(gen_transform_oracle PRIVATE HAVE_CONFIG_H)
set_target_properties(gen_transform_oracle PROPERTIES C_STANDARD 11)
ast_apply_dev_options(gen_transform_oracle)

add_executable(check_transform_oracle check_transform_oracle.c transform_oracle_util.c)
target_include_directories(check_transform_oracle PRIVATE ${TEST_INCLUDES})
target_link_libraries(check_transform_oracle PRIVATE ${TEST_LIBS})
target_compile_definitions(check_transform_oracle PRIVATE HAVE_CONFIG_H)
set_target_properties(check_transform_oracle PROPERTIES C_STANDARD 11)
ast_apply_dev_options(check_transform_oracle)

add_test(NAME transform_oracle_selftest
         COMMAND check_transform_oracle --selftest)
set(_oracle_overrides
    "${CMAKE_CURRENT_SOURCE_DIR}/transform_oracle_overrides.txt")
add_test(NAME transform_oracle_simplify
         COMMAND check_transform_oracle
                 "${CMAKE_CURRENT_SOURCE_DIR}"
                 "${CMAKE_CURRENT_SOURCE_DIR}/simplify_fixtures.oracle"
                 "${_oracle_overrides}")
add_test(NAME transform_oracle_headers
         COMMAND check_transform_oracle
                 "${CMAKE_CURRENT_SOURCE_DIR}"
                 "${CMAKE_CURRENT_SOURCE_DIR}/headers.oracle"
                 "${_oracle_overrides}")
add_test(NAME transform_oracle_framesets
         COMMAND check_transform_oracle
                 "${CMAKE_CURRENT_SOURCE_DIR}"
                 "${CMAKE_CURRENT_SOURCE_DIR}/framesets.oracle"
                 "${_oracle_overrides}")
if(AST_ENABLE_SANITIZERS)
    set_tests_properties(transform_oracle_simplify transform_oracle_headers
                         transform_oracle_framesets
        PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
endif()

# --- Tests using internal headers (#define astCLASS) ---
ast_add_test(testresimp INTERNAL_HEADERS)

# --- Tests needing SAE_PAR stub ---
ast_add_test(testaxis)
ast_add_test(testframe)
ast_add_test(testunitnorm)
ast_add_test(testobjectkeymap)

# --- Tests needing data files ---
ast_add_test(testsplinemap_c)

# --- Conditional tests ---
# --- Fortran tests converted to C (Batch 1) ---
ast_add_test(testzoommap_c    SOURCES testzoommap.c)
ast_add_test(testnormmap_c    SOURCES testnormmap.c)
ast_add_test(testmapping_c    SOURCES testmapping.c)
ast_add_test(testskyframe_c   SOURCES testskyframe.c)
ast_add_test(testcmpframe_c   SOURCES testcmpframe.c)
ast_add_test(testlutmap_c     SOURCES testlutmap.c)
ast_add_test(testratemap_c    SOURCES testratemap.c)
ast_add_test(testchannel_c    SOURCES testchannel.c)

# --- Fortran tests converted to C (Batch 2) ---
ast_add_test(testrate_c       SOURCES testrate.c)
ast_add_test(testspecframe_c  SOURCES testspecframe.c)
ast_add_test(testflux_c       SOURCES testflux.c)
ast_add_test(testspecflux_c   SOURCES testspecflux.c)
ast_add_test(testcmpmap_c     SOURCES testcmpmap.c)
ast_add_test(testpolymap_c    SOURCES testpolymap.c)
ast_add_test(testchebymap_c   SOURCES testchebymap.c)
ast_add_test(testunitnormmap_c SOURCES testunitnormmap.c)
ast_add_test(testtrangrid_c   SOURCES testtrangrid.c)
ast_add_test(testmoc_c        SOURCES testmoc.c)
ast_add_test(testfitstable_c  SOURCES testfitstable.c)
ast_add_test(testframeset_c   SOURCES testframeset.c)
ast_add_test(testswitchmap_c  SOURCES testswitchmap.c)
ast_add_test(testtime_c       SOURCES testtime.c)
ast_add_test(testkeymap_c     SOURCES testkeymap.c)
ast_add_test(testmocchan_c    SOURCES testmocchan.c)
ast_add_test(testxmlchan_c SOURCES testxmlchan.c)
ast_add_test(testtable_c SOURCES testtable.c)
ast_add_test(teststcschan_c SOURCES teststcschan.c)
ast_add_test(testregions_c SOURCES testregions.c)
ast_add_test(testregionmasking_c SOURCES testregionmasking.c)
ast_add_test(testrebinseq_c SOURCES testrebinseq.c)
ast_add_test(testrebin_c SOURCES testrebin.c)
ast_add_test(testresample)
ast_add_test(teststc_c SOURCES teststc.c)
ast_add_test(testfitschan_c SOURCES testfitschan.c)
ast_add_test(testfitschanerrors_c SOURCES testfitschanerrors.c)

# --- New C-only regression tests (no Fortran original) ---
ast_add_test(testslamap)

# --- wcsconverter: regression-diff harness ported from wcsconverter.f ---
add_executable(wcsconverter wcsconverter.c)
target_include_directories(wcsconverter PRIVATE ${TEST_INCLUDES})
target_link_libraries(wcsconverter PRIVATE ${TEST_LIBS})
target_compile_definitions(wcsconverter PRIVATE HAVE_CONFIG_H)
ast_apply_dev_options(wcsconverter)

# ast_astequal: semantic equivalence check for AST object files.
# Used alongside string diffs for platform-robust verification.
add_executable(ast_astequal ast_astequal.c)
target_include_directories(ast_astequal PRIVATE ${TEST_INCLUDES})
target_link_libraries(ast_astequal PRIVATE ${TEST_LIBS})
target_compile_definitions(ast_astequal PRIVATE HAVE_CONFIG_H)
ast_apply_dev_options(ast_astequal)

# Helper to register wcsconverter regression tests.
# Each call copies the input + reference file into the build tree and
# registers two ctest entries:
#   wcsconv_<NAME>          -- byte-level diff of the produced output
#                              against the committed reference.  Catches
#                              unintentional changes to AST's output
#                              format.
#   wcsconv_<NAME>_astequal -- astEqual-based semantic check between the
#                              output and the reference.  Tolerates
#                              platform-specific serialisation and ulp
#                              drift.
# Adding a new fixture is a one-line add_wcsconv_test(...) call below,
# plus dropping the input/reference files into ast_tester/.
#
# Options:
#   SKIP_STRING_COMPARE -- register the astEqual test only; do not
#                          register the byte-level string-diff test.
#                          Use for fixtures whose serialised output is
#                          known to vary at the ulp level between
#                          platforms' math libraries, where a strict
#                          byte-diff is not a reliable regression gate
#                          on any single platform.
function(add_wcsconv_test)
    cmake_parse_arguments(W
        "SKIP_STRING_COMPARE"
        "NAME;INPUT;REFERENCE;ENCODING;ATTRS"
        ""
        ${ARGN})
    if(NOT W_NAME OR NOT W_INPUT OR NOT W_REFERENCE OR NOT W_ENCODING)
        message(FATAL_ERROR
            "add_wcsconv_test: NAME, INPUT, REFERENCE, ENCODING all required")
    endif()
    if(NOT DEFINED W_ATTRS)
        set(W_ATTRS "")
    endif()
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${W_INPUT}"
                   "${CMAKE_CURRENT_BINARY_DIR}/${W_INPUT}" COPYONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${W_REFERENCE}"
                   "${CMAKE_CURRENT_BINARY_DIR}/${W_REFERENCE}" COPYONLY)

    if(NOT W_SKIP_STRING_COMPARE)
        add_test(NAME wcsconv_${W_NAME}
            COMMAND ${CMAKE_COMMAND}
                -DWCSCONVERTER=$<TARGET_FILE:wcsconverter>
                -DIN_FILE=${W_INPUT}
                -DREF_FILE=${W_REFERENCE}
                -DOUT_FILE=${W_NAME}.out
                -DENCODING=${W_ENCODING}
                -DATTRS=${W_ATTRS}
                -P "${CMAKE_SOURCE_DIR}/cmake/run_wcsconverter_test.cmake"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
        if(AST_ENABLE_SANITIZERS)
            set_tests_properties(wcsconv_${W_NAME} PROPERTIES
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
        endif()
    endif()

    add_test(NAME wcsconv_${W_NAME}_astequal
        COMMAND ${CMAKE_COMMAND}
            -DWCSCONVERTER=$<TARGET_FILE:wcsconverter>
            -DAST_ASTEQUAL=$<TARGET_FILE:ast_astequal>
            -DIN_FILE=${W_INPUT}
            -DREF_FILE=${W_REFERENCE}
            -DOUT_FILE=${W_NAME}.astequal.out
            -DENCODING=${W_ENCODING}
            -DATTRS=${W_ATTRS}
            -P "${CMAKE_SOURCE_DIR}/cmake/run_wcsconverter_astequal_test.cmake"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    if(AST_ENABLE_SANITIZERS)
        set_tests_properties(wcsconv_${W_NAME}_astequal PROPERTIES
            ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
    endif()
endfunction()

ast_read_test_data("${CMAKE_CURRENT_SOURCE_DIR}/wcsconv_tests.txt" _wcsconv_tests 6)
foreach(_entry IN LISTS _wcsconv_tests)
    string(REPLACE "|" ";" _fields "${_entry}")
    list(GET _fields 0 _name)
    list(GET _fields 1 _input)
    list(GET _fields 2 _reference)
    list(GET _fields 3 _encoding)
    list(GET _fields 4 _attrs)
    list(GET _fields 5 _skip)
    set(_skip_flag "")
    if(_skip STREQUAL "yes")
        set(_skip_flag "SKIP_STRING_COMPARE")
    endif()
    set(_attrs_arg "")
    if(NOT _attrs STREQUAL "")
        set(_attrs_arg "ATTRS" "${_attrs}")
    endif()
    add_wcsconv_test(NAME ${_name}
                     INPUT ${_input}
                     REFERENCE ${_reference}
                     ENCODING ${_encoding}
                     ${_attrs_arg}
                     ${_skip_flag})
endforeach()

# --- simplify: regression-diff harness ported from simplify.f ---
# Reads a Mapping from an AST dump, simplifies it, and diffs the result
# against a committed reference file.  Same pattern as wcsconverter above.
add_executable(simplify simplify.c)
target_include_directories(simplify PRIVATE ${TEST_INCLUDES})
target_link_libraries(simplify PRIVATE ${TEST_LIBS})
target_compile_definitions(simplify PRIVATE HAVE_CONFIG_H)
ast_apply_dev_options(simplify)

add_executable(normalize_fixture normalize_fixture.c)
target_include_directories(normalize_fixture PRIVATE ${TEST_INCLUDES})
target_link_libraries(normalize_fixture PRIVATE ${TEST_LIBS})
target_compile_definitions(normalize_fixture PRIVATE HAVE_CONFIG_H)
ast_apply_dev_options(normalize_fixture)

add_executable(gen_simplify_fixtures gen_simplify_fixtures.c)
target_include_directories(gen_simplify_fixtures PRIVATE ${TEST_INCLUDES})
target_link_libraries(gen_simplify_fixtures PRIVATE ${TEST_LIBS})
target_compile_definitions(gen_simplify_fixtures PRIVATE HAVE_CONFIG_H)
ast_apply_dev_options(gen_simplify_fixtures)

# Helper to register simplify regression tests.
# Each call copies the input .map and reference .simp files into the build
# tree and registers two ctest entries:
#   simplify_<NAME>          -- byte-level diff of the produced output
#                               against the committed reference.
#   simplify_<NAME>_astequal -- astEqual-based semantic check between the
#                               output and the reference (uses ast_astequal
#                               with encoding "ast").
#
# Options:
#   SKIP_STRING_COMPARE -- register the astEqual test only; do not
#                          register the byte-level string-diff test.
#                          Use for fixtures whose serialised output is
#                          known to vary at the ulp level between
#                          platforms' math libraries.
function(add_simplify_test)
    cmake_parse_arguments(S
        "SKIP_STRING_COMPARE"
        "NAME;INPUT;REFERENCE"
        ""
        ${ARGN})
    if(NOT S_NAME OR NOT S_INPUT OR NOT S_REFERENCE)
        message(FATAL_ERROR
            "add_simplify_test: NAME, INPUT, REFERENCE all required")
    endif()
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${S_INPUT}"
                   "${CMAKE_CURRENT_BINARY_DIR}/${S_INPUT}" COPYONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${S_REFERENCE}"
                   "${CMAKE_CURRENT_BINARY_DIR}/${S_REFERENCE}" COPYONLY)

    if(NOT S_SKIP_STRING_COMPARE)
        add_test(NAME simplify_${S_NAME}
            COMMAND ${CMAKE_COMMAND}
                -DSIMPLIFY=$<TARGET_FILE:simplify>
                -DIN_FILE=${S_INPUT}
                -DREF_FILE=${S_REFERENCE}
                -DOUT_FILE=${S_NAME}.out
                -P "${CMAKE_SOURCE_DIR}/cmake/run_simplify_test.cmake"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
        if(AST_ENABLE_SANITIZERS)
            set_tests_properties(simplify_${S_NAME} PROPERTIES
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
        endif()
    endif()

    add_test(NAME simplify_${S_NAME}_astequal
        COMMAND ${CMAKE_COMMAND}
            -DSIMPLIFY=$<TARGET_FILE:simplify>
            -DAST_ASTEQUAL=$<TARGET_FILE:ast_astequal>
            -DIN_FILE=${S_INPUT}
            -DREF_FILE=${S_REFERENCE}
            -DOUT_FILE=${S_NAME}.astequal.out
            -P "${CMAKE_SOURCE_DIR}/cmake/run_simplify_astequal_test.cmake"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    if(AST_ENABLE_SANITIZERS)
        set_tests_properties(simplify_${S_NAME}_astequal PROPERTIES
            ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
    endif()
endfunction()

ast_read_test_data("${CMAKE_CURRENT_SOURCE_DIR}/simplify_tests.txt" _simplify_tests 4)
foreach(_entry IN LISTS _simplify_tests)
    string(REPLACE "|" ";" _fields "${_entry}")
    list(GET _fields 0 _name)
    list(GET _fields 1 _input)
    list(GET _fields 2 _reference)
    list(GET _fields 3 _skip)
    set(_skip_flag "")
    if(_skip STREQUAL "yes")
        set(_skip_flag "SKIP_STRING_COMPARE")
    endif()
    add_simplify_test(NAME ${_name}
                      INPUT ${_input}
                      REFERENCE ${_reference}
                      ${_skip_flag})
endforeach()

if(AST_ENABLE_HUGE_TEST)
    ast_add_test(testhuge_c SOURCES testhuge.c)
endif()

# --- Conditional tests ---
if(HAVE_YAML)
    ast_add_test(testyamlchan)
endif()

if(HAVE_PTHREADS)
    ast_add_test(testthreads EXTRA_LIBS Threads::Threads)
endif()

# --- Grid smoke tests (no PLplot dependency) ---
# Uses the logging GRF module (grf_log.c) to exercise astGrid on each
# .head file without requiring any external graphics library.
ast_add_test(testgrid
    NO_DEFAULT_GRF_LIBS
    EXTRA_LIBS $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_log>
               $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf3d>
    TEST_ARGS "lsst.fits-wcs" " " " "
)
target_include_directories(testgrid PRIVATE "${CMAKE_SOURCE_DIR}/src")

function(add_grid_test)
    cmake_parse_arguments(P "SKIP_TEXT_COMPARE"
        "NAME;HEAD;ATTR;FATTR;BOX;REFERENCE"
        ""
        ${ARGN})
    if(NOT P_NAME OR NOT P_HEAD)
        message(FATAL_ERROR
            "add_grid_test: NAME, HEAD required")
    endif()
    if(NOT DEFINED P_ATTR)
        set(P_ATTR " ")
    endif()
    if(NOT DEFINED P_FATTR)
        set(P_FATTR " ")
    endif()
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${P_HEAD}"
                   "${CMAKE_CURRENT_BINARY_DIR}/${P_HEAD}" COPYONLY)

    set(_box_arg "")
    if(DEFINED P_BOX)
        set(_box_arg "${P_BOX}")
    endif()

    set(_ref_arg "")
    if(P_REFERENCE AND NOT P_SKIP_TEXT_COMPARE)
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${P_REFERENCE}"
                       "${CMAKE_CURRENT_BINARY_DIR}/${P_REFERENCE}" COPYONLY)
        set(_ref_arg "-DREF_FILE=${P_REFERENCE}")
    endif()

    add_test(NAME grid_${P_NAME}
        COMMAND ${CMAKE_COMMAND}
            -DTESTGRID=$<TARGET_FILE:testgrid>
            -DHEAD_FILE=${P_HEAD}
            "-DATTR=${P_ATTR}"
            "-DFATTR=${P_FATTR}"
            -DOUT_FILE=grid_${P_NAME}.svg
            "${_ref_arg}"
            "-DBOX=${_box_arg}"
            -P "${CMAKE_SOURCE_DIR}/cmake/run_grid_test.cmake"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    if(AST_ENABLE_SANITIZERS)
        set_tests_properties(grid_${P_NAME} PROPERTIES
            ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
    endif()
endfunction()

ast_read_test_data("${CMAKE_CURRENT_SOURCE_DIR}/grid_tests.txt" _grid_tests 7)
foreach(_entry IN LISTS _grid_tests)
    string(REPLACE "|" ";" _fields "${_entry}")
    list(GET _fields 0 _name)
    list(GET _fields 1 _head)
    list(GET _fields 2 _attr)
    list(GET _fields 3 _fattr)
    list(GET _fields 4 _box)
    list(GET _fields 5 _reference)
    list(GET _fields 6 _skip)
    set(_args "NAME" "${_name}" "HEAD" "${_head}")
    if(NOT _attr STREQUAL "")
        list(APPEND _args "ATTR" "${_attr}")
    endif()
    if(NOT _fattr STREQUAL "")
        list(APPEND _args "FATTR" "${_fattr}")
    endif()
    if(NOT _box STREQUAL "")
        list(APPEND _args "BOX" "${_box}")
    endif()
    if(NOT _reference STREQUAL "")
        list(APPEND _args "REFERENCE" "${_reference}")
    endif()
    if(_skip STREQUAL "yes")
        list(APPEND _args "SKIP_TEXT_COMPARE")
    endif()
    add_grid_test(${_args})
endforeach()

if(PLPLOT_FOUND)
    # --- Determine if pdfcairo is available ---
    # If not, plplot falls back on the "pdf" device which uses libHaru, which
    # fails with some of the pdf tests; if the pdfcairo device is not availble
    # then we just mark these tests as expected failures.

    set(CMAKE_REQUIRED_INCLUDES ${PLPLOT_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${PLPLOT_LIBRARIES})
    set(CMAKE_REQUIRED_LINK_DIRECTORIES ${PLPLOT_LIBRARY_DIRS})
    check_c_source_runs("
    #include <string.h>
    #include <plplot.h>

    int main(void) {
        const char *menustr[32];
        const char *devnames[32];
        const char **pmenustr = menustr;
        const char **pdevnames = devnames;
        int idx;
        int pl_ndev = 32;
        plgDevs(&pmenustr, &pdevnames, &pl_ndev);
        for (idx = 0; idx < pl_ndev; idx++) {
            if (devnames[idx] && strcmp(devnames[idx], \"pdfcairo\") == 0)
                return 0;
        }
        return 1;
    }
    " HAVE_PLPLOT_PDFCAIRO)
    unset(CMAKE_REQUIRED_INCLUDES)
    unset(CMAKE_REQUIRED_LIBRARIES)
    unset(CMAKE_REQUIRED_LINK_DIRECTORIES)

    # Plotter tests that are known to fail without pdfcairo
    set(PLOTTER_XFAIL_NO_PDFCAIRO
        car4
        car6
        specflux
        timeplot
    )

    configure_file(lsst.fits-wcs "${CMAKE_CURRENT_BINARY_DIR}/lsst.fits-wcs" COPYONLY)
    ast_add_test(testplotter
        SOURCES testplotter.c plplotutil.c
        NO_DEFAULT_GRF_LIBS
        EXTRA_LIBS $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_plplot> $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf3d_plplot> ${PLPLOT_LIBRARIES}
        TEST_ARGS "${CMAKE_CURRENT_SOURCE_DIR}/lsst.fits-wcs" " " " " "testplotter.ps"
    )
    target_include_directories(testplotter PRIVATE ${PLPLOT_INCLUDE_DIRS})
    target_link_directories(testplotter PRIVATE ${PLPLOT_LIBRARY_DIRS})

    ast_add_test(testplot3d_c
        SOURCES testplot3d.c plplotutil.c
        NO_DEFAULT_GRF_LIBS
        SOURCES testplot3d.c
        EXTRA_LIBS $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_plplot> $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf3d_plplot> ${PLPLOT_LIBRARIES}
        TEST_ARGS "testplot3d.pdf"
    )
    target_include_directories(testplot3d_c PRIVATE ${PLPLOT_INCLUDE_DIRS})
    target_link_directories(testplot3d_c PRIVATE ${PLPLOT_LIBRARY_DIRS})

    # --- Plotter smoke tests ---
    # Run testplotter on each .head file with its .attr/.fattr/.box
    # companions.  Pure smoke test: verifies astGrid does not crash
    # on the given WCS.  No output comparison.
    function(add_plotter_test)
        cmake_parse_arguments(P ""
            "NAME;HEAD;ATTR;FATTR;BOX"
            ""
            ${ARGN})
        if(NOT P_NAME OR NOT P_HEAD)
            message(FATAL_ERROR
                "add_plotter_test: NAME, HEAD required")
        endif()
        if(NOT DEFINED P_ATTR)
            set(P_ATTR " ")
        endif()
        if(NOT DEFINED P_FATTR)
            set(P_FATTR " ")
        endif()
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${P_HEAD}"
                       "${CMAKE_CURRENT_BINARY_DIR}/${P_HEAD}" COPYONLY)

        set(_box_arg "")
        if(DEFINED P_BOX)
            set(_box_arg "${P_BOX}")
        endif()

        add_test(NAME plotter_${P_NAME}
            COMMAND ${CMAKE_COMMAND}
                -DTESTPLOTTER=$<TARGET_FILE:testplotter>
                -DHEAD_FILE=${P_HEAD}
                "-DATTR=${P_ATTR}"
                "-DFATTR=${P_FATTR}"
                -DOUT_FILE=plotter_${P_NAME}.pdf
                "-DBOX=${_box_arg}"
                -P "${CMAKE_SOURCE_DIR}/cmake/run_plotter_test.cmake"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

        if(AST_ENABLE_SANITIZERS)
            set_tests_properties(plotter_${P_NAME} PROPERTIES
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
        endif()

        # Mark test as expected failure if it's one known not to work
        # without pdfcairo
        if(NOT HAVE_PLPLOT_PDFCAIRO)
            list(FIND PLOTTER_XFAIL_NO_PDFCAIRO "${P_NAME}" _idx)
            if(NOT _idx EQUAL -1)
                set_tests_properties(plotter_${P_NAME} PROPERTIES
                    WILL_FAIL TRUE)
            endif()
        endif()
    endfunction()

    ast_read_test_data("${CMAKE_CURRENT_SOURCE_DIR}/grid_tests.txt" _plotter_tests 7)
    foreach(_entry IN LISTS _plotter_tests)
        string(REPLACE "|" ";" _fields "${_entry}")
        list(GET _fields 0 _name)
        list(GET _fields 1 _head)
        list(GET _fields 2 _attr)
        list(GET _fields 3 _fattr)
        list(GET _fields 4 _box)
        set(_args "NAME" "${_name}" "HEAD" "${_head}")
        if(NOT _attr STREQUAL "")
            list(APPEND _args "ATTR" "${_attr}")
        endif()
        if(NOT _fattr STREQUAL "")
            list(APPEND _args "FATTR" "${_fattr}")
        endif()
        if(NOT _box STREQUAL "")
            list(APPEND _args "BOX" "${_box}")
        endif()
        add_plotter_test(${_args})
    endforeach()
endif()


# ============================================================
# Fortran tests (AST_BUILD_FORTRAN=ON only)
# ============================================================
# These run the original Fortran test programs against the Fortran binding,
# complementing the C equivalents above.  Each test is suffixed _f to
# distinguish it from the C version in ctest output.
if(AST_BUILD_FORTRAN)

    # Minimal C implementations of the Starlink APIs called by the Fortran
    # tests (EMS error handling, PSX memory, CNF C/Fortran interface, etc.).
    add_library(ast_star_compat STATIC
        star_compat/err_standalone.c
        star_compat/psx_standalone.c
        star_compat/chr_standalone.c
        star_compat/cnf_standalone.c
        star_compat/num_standalone.c
    )
    target_include_directories(ast_star_compat PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/star_compat"
        "${CMAKE_BINARY_DIR}"       # config.h, ast_err.h
        "${CMAKE_BINARY_DIR}/src"   # f77.h (generated from f77.h.in)
    )
    target_compile_definitions(ast_star_compat PRIVATE HAVE_CONFIG_H)

    # Helper: add a Fortran test program and register it with ctest.
    # Test names get a _f suffix to distinguish from C equivalents.
    # The Fortran INCLUDE search path covers:
    #   - star_compat/  : SAE_PAR, SAI_ERR, CNF_PAR, ERR_PAR, etc.
    #   - CMAKE_BINARY_DIR : generated AST_PAR and AST_ERR
    function(ast_add_fortran_test name)
        cmake_parse_arguments(FT "NO_DEFAULT_GRF_LIBS" "" "EXTRA_LIBS;TEST_ARGS" ${ARGN})
        add_executable(${name}_f ${name}.f)
        if(AST_FORTRAN_TEST_FLAGS)
            target_compile_options(${name}_f PRIVATE
                $<$<COMPILE_LANGUAGE:Fortran>:${AST_FORTRAN_TEST_FLAGS}>
            )
        endif()
        target_include_directories(${name}_f PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/star_compat"
            "${CMAKE_BINARY_DIR}"   # AST_PAR, AST_ERR
        )
        if(FT_NO_DEFAULT_GRF_LIBS)
            set(_grf_libs
                $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_err>
                ${FT_EXTRA_LIBS}
            )
        else()
            set(_grf_libs
                $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_err>
                $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_2.0>
                $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_3.2>
                $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_5.6>
                $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf3d>
                ${FT_EXTRA_LIBS}
            )
        endif()
        target_link_libraries(${name}_f PRIVATE
            ast
            ${_grf_libs}
            ast_star_compat
            m
        )
        add_test(NAME ${name}_f COMMAND ${name}_f ${FT_TEST_ARGS}
                 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

        set(test_env "")
        if(AST_ENABLE_SANITIZERS)
            list(APPEND test_env "ASAN_OPTIONS=detect_leaks=0")
        endif()

        # Tests locate their input fixture files relative to $srcdir, so point
        # it at the source directory where those fixtures live.
        list(APPEND test_env "srcdir=${CMAKE_CURRENT_SOURCE_DIR}")

        set_tests_properties(${name}_f PROPERTIES
            ENVIRONMENT "${test_env}")
        endfunction()

    # Custom target that ensures AST_PAR and AST_ERR exist before any
    # Fortran test is compiled (they are Fortran INCLUDE files).
    # ast_par and ast_codegen are defined in the top-level CMakeLists.txt.
    add_custom_target(ast_fortran_includes)
    add_dependencies(ast_fortran_includes ast_par ast_codegen)

    # All 33 unconditional Fortran test programs
    set(_fortran_tests
        testzoommap
        testnormmap
        testmapping
        testskyframe
        testcmpframe
        testlutmap
        testratemap
        testchannel
        testrate
        testspecframe
        testflux
        testspecflux
        testcmpmap
        testpolymap
        testchebymap
        testunitnormmap
        testtrangrid
        testmoc
        testfitstable
        testframeset
        testswitchmap
        testtime
        testkeymap
        testmocchan
        testxmlchan
        testtable
        teststcschan
        testregions
        testrebinseq
        testrebin
        teststc
        testfitschan
        testsplinemap
    )
    foreach(_f_test IN LISTS _fortran_tests)
        ast_add_fortran_test(${_f_test})
        add_dependencies(${_f_test}_f ast_fortran_includes)
    endforeach()

    if(AST_ENABLE_HUGE_TEST)
        ast_add_fortran_test(testhuge)
        add_dependencies(testhuge_f ast_fortran_includes)
    endif()
endif() # AST_BUILD_FORTRAN
