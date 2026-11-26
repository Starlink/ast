# Driver script for ctest Plot3D logging tests.  Invoked via `cmake -P`
# by add_plot3dlog_test().
#
# Runs testplot3dlog with the logging grf3d module, then compares the
# grf3d call log byte for byte with REF_LOG and the Native dump of the
# Plot3D with REF_DUMP through compare_dumps.  When CAMERA is set an SVG
# is written as well and checked to be a complete SVG document.
#
# Required cache variables (supplied via -D):
#   PROG      : absolute path to the testplot3dlog executable
#   COMPARE   : absolute path to the compare_dumps executable
#   IN_FILE   : input fixture (resolved against CWD)
#   ATTR      : Plot3D attribute settings (may be empty)
#   BOX       : graphics box "x1 y1 z1 x2 y2 z2"
#   BBOX      : physical box "x1 y1 z1 x2 y2 z2"
#   MM        : millimetres per world unit (may be empty)
#   ACTIONS   : space-separated actions
#   CAMERA    : empty, "auto", or ten numbers
#   OUT_LOG, OUT_DUMP, OUT_SVG : output files to produce
#   REF_LOG, REF_DUMP          : reference files

foreach(_req IN ITEMS PROG COMPARE IN_FILE BOX BBOX ACTIONS
                      OUT_LOG OUT_DUMP OUT_SVG REF_LOG REF_DUMP)
    if(NOT DEFINED ${_req})
        message(FATAL_ERROR "run_plot3dlog_test.cmake: ${_req} not set")
    endif()
endforeach()

set(_cmd "${PROG}" --in "${IN_FILE}" --box "${BOX}" --bbox "${BBOX}"
         --log "${OUT_LOG}" --dump "${OUT_DUMP}")
if(DEFINED ATTR AND NOT ATTR STREQUAL "")
    list(APPEND _cmd --attr "${ATTR}")
endif()
if(DEFINED MM AND NOT MM STREQUAL "")
    list(APPEND _cmd --mm "${MM}")
endif()
if(DEFINED CAMERA AND NOT CAMERA STREQUAL "")
    list(APPEND _cmd --svg "${OUT_SVG}")
    if(NOT CAMERA STREQUAL "auto")
        list(APPEND _cmd --camera "${CAMERA}")
    endif()
endif()
separate_arguments(_actions UNIX_COMMAND "${ACTIONS}")
foreach(_a IN LISTS _actions)
    list(APPEND _cmd --do "${_a}")
endforeach()

execute_process(
    COMMAND ${_cmd}
    RESULT_VARIABLE _rv
)
if(NOT _rv EQUAL 0)
    message(FATAL_ERROR "testplot3dlog exited with code ${_rv} "
                        "(cmd: ${_cmd})")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files "${REF_LOG}" "${OUT_LOG}"
    RESULT_VARIABLE _rv
)
if(NOT _rv EQUAL 0)
    message(STATUS "grf3d log mismatch; showing diff:")
    execute_process(COMMAND diff -u "${REF_LOG}" "${OUT_LOG}")
    message(FATAL_ERROR
        "grf3d log ${OUT_LOG} differs from reference ${REF_LOG}")
endif()

execute_process(
    COMMAND "${COMPARE}" "${REF_DUMP}" "${OUT_DUMP}"
    RESULT_VARIABLE _rv
)
if(NOT _rv EQUAL 0)
    message(STATUS "Plot3D dump mismatch; showing diff:")
    execute_process(COMMAND diff -u "${REF_DUMP}" "${OUT_DUMP}")
    message(FATAL_ERROR
        "Plot3D dump ${OUT_DUMP} differs from reference ${REF_DUMP}")
endif()

if(DEFINED CAMERA AND NOT CAMERA STREQUAL "")
    file(READ "${OUT_SVG}" _svg)
    if(NOT _svg MATCHES "<svg " OR NOT _svg MATCHES "</svg>")
        message(FATAL_ERROR "${OUT_SVG} is not a complete SVG document")
    endif()
endif()
