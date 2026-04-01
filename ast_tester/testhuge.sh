#!/bin/sh
# Wrapper for testhuge: skip unless TEST_HUGE=1 is passed to make.
# testhuge takes a long time to run and significant memory resources;
# it is omitted from the default test run.
# Enable with: make check TEST_HUGE=1
if test "x${TEST_HUGE}" != "x1"; then
    echo "testhuge: skipped (re-run with TEST_HUGE=1 to enable)"
    exit 77
fi
exec ./testhuge "$@"
