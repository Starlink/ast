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

# All committed test data lives below a single fixture root.  Stage every
# fixture at the same relative path below the build-tree ast_tester directory
# so source-tree and CTest runs use the same "fixtures/<relative-path>"
# convention.
set(AST_FIXTURE_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
set(AST_FIXTURE_BINARY_ROOT "${CMAKE_CURRENT_BINARY_DIR}/fixtures")

# Copy the whole tree in one build-time rule rather than one configure_file()
# per fixture: a configure_file() would make every one of the ~1200 fixtures a
# configure dependency, so regenerating a single .oracle would reconfigure the
# project.  The stamp file keeps the copy off the critical path of an
# unchanged build.  Adding or removing a fixture changes the glob and so does
# reconfigure, which is wanted -- a new fixture normally needs a new test
# registered as well.
file(GLOB_RECURSE _ast_fixture_files CONFIGURE_DEPENDS
     RELATIVE "${AST_FIXTURE_SOURCE_ROOT}"
     "${AST_FIXTURE_SOURCE_ROOT}/*")
set(_ast_fixture_deps)
foreach(_fixture IN LISTS _ast_fixture_files)
    list(APPEND _ast_fixture_deps "${AST_FIXTURE_SOURCE_ROOT}/${_fixture}")
endforeach()

set(_ast_fixture_stamp "${CMAKE_CURRENT_BINARY_DIR}/fixtures.stamp")
add_custom_command(OUTPUT "${_ast_fixture_stamp}"
    COMMAND "${CMAKE_COMMAND}" -E copy_directory
            "${AST_FIXTURE_SOURCE_ROOT}" "${AST_FIXTURE_BINARY_ROOT}"
    COMMAND "${CMAKE_COMMAND}" -E touch "${_ast_fixture_stamp}"
    DEPENDS ${_ast_fixture_deps}
    COMMENT "Staging ast_tester fixtures"
    VERBATIM)
add_custom_target(ast_stage_fixtures ALL DEPENDS "${_ast_fixture_stamp}")

# Fail at configure time, naming the manifest, when a registered fixture path
# does not exist.  The per-test configure_file() calls that used to build
# these paths gave this check as a side effect; without it a typo in a
# cases.txt surfaces only as a harness failure part-way through a ctest run,
# with nothing to say the file was never there.
function(ast_require_fixture manifest relative_path)
    if(NOT EXISTS "${AST_FIXTURE_SOURCE_ROOT}/${relative_path}")
        message(FATAL_ERROR
            "${manifest}: no such fixture: fixtures/${relative_path}")
    endif()
endfunction()

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
    # So that building a single test target also stages the fixtures it reads.
    add_dependencies(${name} ast_stage_fixtures)
    add_test(NAME ${name} COMMAND ${name} ${T_TEST_ARGS}
             WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    # AST's memory system intentionally caches freed blocks for reuse,
    # so LeakSanitizer reports false positives. Disable leak detection.
    if(AST_ENABLE_SANITIZERS)
        set_tests_properties(${name} PROPERTIES
            ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
    endif()
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
    ENVIRONMENT "ORACLE_TEST_ROOT=${AST_FIXTURE_BINARY_ROOT}")

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
    "${AST_FIXTURE_BINARY_ROOT}/oracle/transform_oracle_overrides.txt")
add_test(NAME transform_oracle_simplify
         COMMAND check_transform_oracle
                 "${AST_FIXTURE_BINARY_ROOT}"
                 "${AST_FIXTURE_BINARY_ROOT}/oracle/simplify_fixtures.oracle"
                 "${_oracle_overrides}")
add_test(NAME transform_oracle_headers
         COMMAND check_transform_oracle
                 "${AST_FIXTURE_BINARY_ROOT}"
                 "${AST_FIXTURE_BINARY_ROOT}/oracle/headers.oracle"
                 "${_oracle_overrides}")
add_test(NAME transform_oracle_framesets
         COMMAND check_transform_oracle
                 "${AST_FIXTURE_BINARY_ROOT}"
                 "${AST_FIXTURE_BINARY_ROOT}/oracle/framesets.oracle"
                 "${_oracle_overrides}")
if(AST_ENABLE_SANITIZERS)
    set_tests_properties(transform_oracle_simplify transform_oracle_headers
                         transform_oracle_framesets
        PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
endif()

# --- Tests using internal headers (#define astCLASS) ---
ast_add_test(testresimp INTERNAL_HEADERS)
ast_add_test(testissimp INTERNAL_HEADERS)

# --- Tests needing SAE_PAR stub ---
ast_add_test(testaxis)
ast_add_test(testframe)
ast_add_test(testunitnorm)
ast_add_test(testobjectkeymap)

# --- Tests needing data files ---
ast_add_test(testsplinemap_c)

# --- Conditional tests ---
# --- Fortran tests converted to C (Batch 1) ---
ast_add_test(testzoommap)
ast_add_test(testimmutable)
ast_add_test(testnormmap)
ast_add_test(testmapping)
ast_add_test(testskyframe)
ast_add_test(testcmpframe)
ast_add_test(testlutmap)
ast_add_test(testratemap)
ast_add_test(testchannel)

# --- Fortran tests converted to C (Batch 2) ---
ast_add_test(testrate)
ast_add_test(testspecframe)
ast_add_test(testflux)
ast_add_test(testspecflux)
ast_add_test(testcmpmap)
ast_add_test(testpolymap)
ast_add_test(testchebymap)
ast_add_test(testunitnormmap_c SOURCES testunitnormmap.c)
ast_add_test(testtrangrid)
ast_add_test(testmoc)
ast_add_test(testfitstable)
ast_add_test(testframeset)
ast_add_test(testswitchmap)
ast_add_test(testtime)
ast_add_test(testkeymap)
ast_add_test(testmocchan)
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
# Each call registers two ctest entries using fixture-root-relative paths:
#   wcsconv_<NAME>          -- byte-level diff of the produced output
#                              against the committed reference.  Catches
#                              unintentional changes to AST's output
#                              format.
#   wcsconv_<NAME>_astequal -- astEqual-based semantic check between the
#                              output and the reference.  Tolerates
#                              platform-specific serialisation and ulp
#                              drift.
# Add fixtures below fixtures/ and list their fixture-root-relative paths in
# fixtures/wcsconv/cases.txt.
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
        "NAME;INPUT;REFERENCE;ENCODING;ATTRS;TOLERANCE"
        ""
        ${ARGN})
    if(NOT W_NAME OR NOT W_INPUT OR NOT W_REFERENCE OR NOT W_ENCODING)
        message(FATAL_ERROR
            "add_wcsconv_test: NAME, INPUT, REFERENCE, ENCODING all required")
    endif()
    if(NOT DEFINED W_ATTRS)
        set(W_ATTRS "")
    endif()
    ast_require_fixture("wcsconv/cases.txt" "${W_INPUT}")
    ast_require_fixture("wcsconv/cases.txt" "${W_REFERENCE}")
    set(_input "fixtures/${W_INPUT}")
    set(_reference "fixtures/${W_REFERENCE}")

    if(NOT W_SKIP_STRING_COMPARE)
        add_test(NAME wcsconv_${W_NAME}
            COMMAND ${CMAKE_COMMAND}
                -DWCSCONVERTER=$<TARGET_FILE:wcsconverter>
                -DCOMPARE=$<TARGET_FILE:compare_dumps>
                -DTOLERANCE=${W_TOLERANCE}
                -DIN_FILE=${_input}
                -DREF_FILE=${_reference}
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
            -DIN_FILE=${_input}
            -DREF_FILE=${_reference}
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

ast_read_test_data("${AST_FIXTURE_SOURCE_ROOT}/wcsconv/cases.txt"
                   _wcsconv_tests 6)
foreach(_entry IN LISTS _wcsconv_tests)
    string(REPLACE "|" ";" _fields "${_entry}")
    list(GET _fields 0 _name)
    list(GET _fields 1 _input)
    list(GET _fields 2 _reference)
    list(GET _fields 3 _encoding)
    list(GET _fields 4 _attrs)
    list(GET _fields 5 _skip)
    # Three values.  Blank compares the output with numeric tokens equal to
    # within 4 ULP.  "astequal" widens that to AST's own astEQUAL, for a value
    # amplified by cancellation -- OBSALT recovered from OBSGEO-X/Y/Z is the case
    # that needs it.  "yes" runs no string comparison at all, and needs a comment
    # above the row saying what varies.
    set(_skip_flag "")
    set(_tolerance "")
    if(_skip STREQUAL "yes")
        set(_skip_flag "SKIP_STRING_COMPARE")
    elseif(_skip STREQUAL "astequal")
        set(_tolerance "astequal")
    elseif(NOT _skip STREQUAL "")
        message(FATAL_ERROR
            "wcsconv/cases.txt: row ${_name} has skip_string_compare "
            "\"${_skip}\"; expected blank, astequal or yes")
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
                     TOLERANCE ${_tolerance}
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

# --- roundtrip: read a dump and write it straight back ---
#
# The output must be byte-identical to the input.  This is the only family that
# can see an attribute lost on load: both simplify families call astSimplify
# after reading, and astSimplify sets IsSimple on every node it visits, so a
# record a loader dropped is put straight back before anything compares it.
add_executable(roundtrip roundtrip.c)
target_include_directories(roundtrip PRIVATE ${TEST_INCLUDES})
target_link_libraries(roundtrip PRIVATE ${TEST_LIBS})
target_compile_definitions(roundtrip PRIVATE HAVE_CONFIG_H)
ast_apply_dev_options(roundtrip)
add_dependencies(roundtrip ast_stage_fixtures)
add_dependencies(simplify compare_dumps)

function(add_roundtrip_test)
    cmake_parse_arguments(R "" "NAME;FIXTURE" "" ${ARGN})
    if(NOT R_NAME OR NOT R_FIXTURE)
        message(FATAL_ERROR "add_roundtrip_test: NAME and FIXTURE required")
    endif()
    ast_require_fixture("roundtrip" "${R_FIXTURE}")
    add_test(NAME roundtrip_${R_NAME}
        COMMAND "${CMAKE_COMMAND}"
                -DROUNDTRIP=$<TARGET_FILE:roundtrip>
                -DIN_FILE=${AST_FIXTURE_BINARY_ROOT}/${R_FIXTURE}
                -DOUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/roundtrip_${R_NAME}.out
                -P "${CMAKE_SOURCE_DIR}/cmake/run_roundtrip_test.cmake"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endfunction()

# compare_dumps: a dump comparator that tolerates last-bit variation.  It links
# no AST library -- it is text and arithmetic only -- so a cmake -P driver can
# call it without the library on the path.
add_executable(compare_dumps compare_dumps.c)
target_link_libraries(compare_dumps PRIVATE m)
set_target_properties(compare_dumps PROPERTIES C_STANDARD 11)
ast_apply_dev_options(compare_dumps)

add_executable(compare_dumps_selftest compare_dumps.c)
target_compile_definitions(compare_dumps_selftest PRIVATE SELFTEST)
target_link_libraries(compare_dumps_selftest PRIVATE m)
set_target_properties(compare_dumps_selftest PROPERTIES C_STANDARD 11)
ast_apply_dev_options(compare_dumps_selftest)
add_test(NAME compare_dumps_selftest COMMAND compare_dumps_selftest)

add_executable(gen_keymap_mapsz gen_keymap_mapsz.c)
target_include_directories(gen_keymap_mapsz PRIVATE ${TEST_INCLUDES})
target_link_libraries(gen_keymap_mapsz PRIVATE ${TEST_LIBS})
target_compile_definitions(gen_keymap_mapsz PRIVATE HAVE_CONFIG_H)
ast_apply_dev_options(gen_keymap_mapsz)

add_executable(gen_frameset_fixtures gen_frameset_fixtures.c)
target_include_directories(gen_frameset_fixtures PRIVATE ${TEST_INCLUDES})
target_link_libraries(gen_frameset_fixtures PRIVATE ${TEST_LIBS})
target_compile_definitions(gen_frameset_fixtures PRIVATE HAVE_CONFIG_H)
set_target_properties(gen_frameset_fixtures PROPERTIES C_STANDARD 11)
ast_apply_dev_options(gen_frameset_fixtures)

add_executable(gen_native_fixtures gen_native_fixtures.c)
target_include_directories(gen_native_fixtures PRIVATE ${TEST_INCLUDES})
target_link_libraries(gen_native_fixtures PRIVATE ${TEST_LIBS})
target_compile_definitions(gen_native_fixtures PRIVATE HAVE_CONFIG_H)
set_target_properties(gen_native_fixtures PROPERTIES C_STANDARD 11)
ast_apply_dev_options(gen_native_fixtures)

add_executable(gen_simplify_fixtures gen_simplify_fixtures.c gen_simplify_flags.c)
target_include_directories(gen_simplify_fixtures PRIVATE ${TEST_INCLUDES})
target_link_libraries(gen_simplify_fixtures PRIVATE ${TEST_LIBS})
target_compile_definitions(gen_simplify_fixtures PRIVATE HAVE_CONFIG_H)
# gen_simplify_flags.c reaches the protected Mapping flags through the internal
# headers, so it must see the same AstObject / AstMapping layout the library was
# compiled with -- otherwise the flag macros write at the wrong offset.
if(HAVE_PTHREADS)
    target_compile_definitions(gen_simplify_fixtures PRIVATE THREAD_SAFE)
endif()
ast_apply_dev_options(gen_simplify_fixtures)

# Helper to register simplify regression tests.
# Each call registers tests using fixture-root-relative input/reference paths:
#   simplify_<NAME>          -- byte-level diff of the produced output
#                               against the committed reference.
#   simplify_<NAME>_astequal -- astEqual-based semantic check between the
#                               output and the reference (uses ast_astequal
#                               with encoding "ast").
#
# A fixture whose REFERENCE is its own INPUT is a negative control: it asserts
# astSimplify leaves the Mapping alone.  Such a fixture gets a third entry:
#   simplify_<NAME>_noop     -- byte-level diff of the output against the
#                               input with the IsSimp markers removed.
# This is registered even under SKIP_STRING_COMPARE, because astEqual compares
# meaning and so passes a structural rewrite that preserves the transformation
# -- exactly the failure a negative control exists to catch.
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
        "NAME;INPUT;REFERENCE;TOLERANCE"
        ""
        ${ARGN})
    if(NOT S_NAME OR NOT S_INPUT OR NOT S_REFERENCE)
        message(FATAL_ERROR
            "add_simplify_test: NAME, INPUT, REFERENCE all required")
    endif()
    ast_require_fixture("simplify/cases.txt" "${S_INPUT}")
    ast_require_fixture("simplify/cases.txt" "${S_REFERENCE}")
    set(_input "fixtures/${S_INPUT}")
    set(_reference "fixtures/${S_REFERENCE}")

    if(NOT S_SKIP_STRING_COMPARE)
        add_test(NAME simplify_${S_NAME}
            COMMAND ${CMAKE_COMMAND}
                -DSIMPLIFY=$<TARGET_FILE:simplify>
                -DCOMPARE=$<TARGET_FILE:compare_dumps>
                -DTOLERANCE=${S_TOLERANCE}
                -DIN_FILE=${_input}
                -DREF_FILE=${_reference}
                -DOUT_FILE=${S_NAME}.out
                -P "${CMAKE_SOURCE_DIR}/cmake/run_simplify_test.cmake"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
        if(AST_ENABLE_SANITIZERS)
            set_tests_properties(simplify_${S_NAME} PROPERTIES
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
        endif()
    endif()

    if(S_INPUT STREQUAL S_REFERENCE)
        add_test(NAME simplify_${S_NAME}_noop
            COMMAND ${CMAKE_COMMAND}
                -DSIMPLIFY=$<TARGET_FILE:simplify>
                -DIN_FILE=${_input}
                -DOUT_FILE=${S_NAME}.noop.out
                -P "${CMAKE_SOURCE_DIR}/cmake/run_simplify_noop_test.cmake"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
        if(AST_ENABLE_SANITIZERS)
            set_tests_properties(simplify_${S_NAME}_noop PROPERTIES
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
        endif()
    endif()

    add_test(NAME simplify_${S_NAME}_astequal
        COMMAND ${CMAKE_COMMAND}
            -DSIMPLIFY=$<TARGET_FILE:simplify>
            -DAST_ASTEQUAL=$<TARGET_FILE:ast_astequal>
            -DIN_FILE=${_input}
            -DREF_FILE=${_reference}
            -DOUT_FILE=${S_NAME}.astequal.out
            -P "${CMAKE_SOURCE_DIR}/cmake/run_simplify_astequal_test.cmake"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    if(AST_ENABLE_SANITIZERS)
        set_tests_properties(simplify_${S_NAME}_astequal PROPERTIES
            ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
    endif()
endfunction()

ast_read_test_data("${AST_FIXTURE_SOURCE_ROOT}/simplify/cases.txt"
                   _simplify_tests 4)
foreach(_entry IN LISTS _simplify_tests)
    string(REPLACE "|" ";" _fields "${_entry}")
    list(GET _fields 0 _name)
    list(GET _fields 1 _input)
    list(GET _fields 2 _reference)
    list(GET _fields 3 _skip)
    # The fourth column takes three values.  Blank compares the dump with
    # numeric tokens to 4 ULP.  "astequal" widens that to AST's own astEQUAL,
    # for a value amplified by cancellation.  "yes" runs no string comparison at
    # all, which is right for a negative fixture: its reference is the
    # unsimplified input, and astSimplify always stamps IsSimp on its result, so
    # the two can never match by string.  Those are gated by the _noop test
    # instead, which strips IsSimp before comparing.
    set(_skip_flag "")
    set(_tolerance "")
    if(_skip STREQUAL "yes")
        set(_skip_flag "SKIP_STRING_COMPARE")
    elseif(_skip STREQUAL "astequal")
        set(_tolerance "astequal")
    elseif(NOT _skip STREQUAL "")
        message(FATAL_ERROR
            "simplify/cases.txt: row ${_name} has skip_string_compare "
            "\"${_skip}\"; expected blank, astequal or yes")
    endif()
    add_simplify_test(NAME ${_name}
                      INPUT ${_input}
                      REFERENCE ${_reference}
                      TOLERANCE ${_tolerance}
                      ${_skip_flag})

# Round-trip both halves of the pair.  Reading a dump and writing it back is
# the identity whatever the fixture is, so this needs no reference of its own
# and no skip flag: a fixture that cannot round-trip is a serialisation defect.
    add_roundtrip_test(NAME ${_name}_map      FIXTURE ${_input})
    if(NOT _reference STREQUAL _input)
        add_roundtrip_test(NAME ${_name}_simp FIXTURE ${_reference})
    endif()
endforeach()

# Round-trip every other native-encoded fixture in the tree.  These are inputs
# to individual programs rather than manifest cases, so they are discovered by
# glob: a new one is picked up on the next reconfigure.
#
# wcsconv/ is deliberately excluded.  Its references are FITS-encoded, so
# reading one back needs a FitsChan and the encoding it was written with, which
# is what the wcsconv harness already does.
file(GLOB_RECURSE _ast_native_fixtures CONFIGURE_DEPENDS
     RELATIVE "${AST_FIXTURE_SOURCE_ROOT}"
     "${AST_FIXTURE_SOURCE_ROOT}/programs/*.ast"
     "${AST_FIXTURE_SOURCE_ROOT}/serialisation/*.ast"
     "${AST_FIXTURE_SOURCE_ROOT}/wcsconv/framesets/*.ast")
foreach(_fixture IN LISTS _ast_native_fixtures)
    string(REGEX REPLACE "[/.]" "_" _rt_name "${_fixture}")
    add_roundtrip_test(NAME ${_rt_name} FIXTURE ${_fixture})
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
    TEST_ARGS "fixtures/wcsconv/inputs/lsst.fits-wcs" " " " "
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
    ast_require_fixture("plot/cases.txt" "${P_HEAD}")
    set(_head "fixtures/${P_HEAD}")

    set(_box_arg "")
    if(DEFINED P_BOX)
        set(_box_arg "${P_BOX}")
    endif()

    set(_ref_arg "")
    if(P_REFERENCE AND NOT P_SKIP_TEXT_COMPARE)
        ast_require_fixture("plot/cases.txt" "${P_REFERENCE}")
        set(_ref_arg "-DREF_FILE=fixtures/${P_REFERENCE}")
    endif()

    add_test(NAME grid_${P_NAME}
        COMMAND ${CMAKE_COMMAND}
            -DTESTGRID=$<TARGET_FILE:testgrid>
            -DHEAD_FILE=${_head}
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

ast_read_test_data("${AST_FIXTURE_SOURCE_ROOT}/plot/cases.txt" _grid_tests 7)
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

    ast_add_test(testplotter
        SOURCES testplotter.c plplotutil.c
        NO_DEFAULT_GRF_LIBS
        EXTRA_LIBS $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_plplot> $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf3d_plplot> ${PLPLOT_LIBRARIES}
        TEST_ARGS "fixtures/wcsconv/inputs/lsst.fits-wcs" " " " " "testplotter.ps"
    )
    target_include_directories(testplotter PRIVATE ${PLPLOT_INCLUDE_DIRS})
    target_link_directories(testplotter PRIVATE ${PLPLOT_LIBRARY_DIRS})

    ast_add_test(testplot3d
        SOURCES testplot3d.c plplotutil.c
        NO_DEFAULT_GRF_LIBS
        EXTRA_LIBS $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf_plplot> $<LINK_LIBRARY:WHOLE_ARCHIVE,ast_grf3d_plplot> ${PLPLOT_LIBRARIES}
        TEST_ARGS "testplot3d.pdf"
    )
    target_include_directories(testplot3d PRIVATE ${PLPLOT_INCLUDE_DIRS})
    target_link_directories(testplot3d PRIVATE ${PLPLOT_LIBRARY_DIRS})

    # --- Plotter smoke tests ---
    # Run testplotter using the inline settings in fixtures/plot/cases.txt.
    # Pure smoke test: verifies astGrid does not crash on the given WCS.
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
        ast_require_fixture("plot/cases.txt" "${P_HEAD}")
        set(_head "fixtures/${P_HEAD}")

        set(_box_arg "")
        if(DEFINED P_BOX)
            set(_box_arg "${P_BOX}")
        endif()

        add_test(NAME plotter_${P_NAME}
            COMMAND ${CMAKE_COMMAND}
                -DTESTPLOTTER=$<TARGET_FILE:testplotter>
                -DHEAD_FILE=${_head}
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

    ast_read_test_data("${AST_FIXTURE_SOURCE_ROOT}/plot/cases.txt"
                       _plotter_tests 7)
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

