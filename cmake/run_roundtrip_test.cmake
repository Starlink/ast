# Driver for the roundtrip ctests: read a fixture, dump it straight back, and
# require the output to be byte-identical to the input.
#
# Reading and re-writing a dump is the identity, so any difference is a
# serialisation defect.  This is the only family that can see one: both
# simplify families call astSimplify after reading, and astSimplify re-sets
# IsSimple on every node it visits, hiding a record a loader dropped.
#
# Required cache variables (via -D): ROUNDTRIP, IN_FILE, OUT_FILE

foreach(_req IN ITEMS ROUNDTRIP IN_FILE OUT_FILE)
    if(NOT DEFINED ${_req})
        message(FATAL_ERROR "run_roundtrip_test.cmake: ${_req} not set")
    endif()
endforeach()

execute_process(COMMAND "${ROUNDTRIP}" "${IN_FILE}" "${OUT_FILE}"
                RESULT_VARIABLE _rv)
if(NOT _rv EQUAL 0)
    message(FATAL_ERROR "roundtrip exited with code ${_rv} "
                        "(cmd: ${ROUNDTRIP} ${IN_FILE} ${OUT_FILE})")
endif()

execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
                        "${IN_FILE}" "${OUT_FILE}"
                RESULT_VARIABLE _rv)
if(NOT _rv EQUAL 0)
    message(STATUS "Round-trip mismatch; showing unified diff:")
    execute_process(COMMAND diff -u "${IN_FILE}" "${OUT_FILE}")
    message(FATAL_ERROR
        "reading ${IN_FILE} and writing it back did not reproduce it")
endif()
