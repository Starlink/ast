# Driver script for ctest cases that run the simplify tool and diff its
# output against a reference file.  Invoked via `cmake -P` by
# add_simplify_test().
#
# Required cache variables (supplied via -D):
#   SIMPLIFY  : absolute path to the simplify executable
#   IN_FILE   : path to the input .map file
#   REF_FILE  : reference .simp file to diff against
#   OUT_FILE  : output file to produce

foreach(_req IN ITEMS SIMPLIFY IN_FILE REF_FILE OUT_FILE)
    if(NOT DEFINED ${_req})
        message(FATAL_ERROR "run_simplify_test.cmake: ${_req} not set")
    endif()
endforeach()

execute_process(
    COMMAND "${SIMPLIFY}" "${IN_FILE}" "${OUT_FILE}"
    RESULT_VARIABLE _rv
)
if(NOT _rv EQUAL 0)
    message(FATAL_ERROR "simplify exited with code ${_rv} "
                        "(cmd: ${SIMPLIFY} ${IN_FILE} ${OUT_FILE})")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files "${REF_FILE}" "${OUT_FILE}"
    RESULT_VARIABLE _rv
)
if(NOT _rv EQUAL 0)
    message(STATUS "Output mismatch; showing unified diff:")
    execute_process(
        COMMAND diff -u "${REF_FILE}" "${OUT_FILE}"
    )
    message(FATAL_ERROR
        "simplify output ${OUT_FILE} differs from reference ${REF_FILE}")
endif()
