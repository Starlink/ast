#!/bin/sh
# Verify fixture-status classification: a deliberate omission must configure
# cleanly, a missing manifest-listed fixture must fail, and a git checkout must
# keep failing on any missing fixture.
set -e
cd "$(dirname "$0")/.."
FIX=ast_tester/fixtures
STASH=$(mktemp -d)
fail=0

cleanup() {
    rm -f "$FIX/DIST_TIER"
    test -d "$STASH/expected" && mv "$STASH/expected" "$FIX/wcsconv/expected"
    test -f "$STASH/aitoff.svg" && \
        mv "$STASH/aitoff.svg" "$FIX/plot/expected/aitoff.svg"
    rm -rf "$STASH" _tiertest
    return 0
}
trap cleanup EXIT

say() { printf '%s\n' "$*"; }
check() {
    if [ "$1" = "$2" ]; then say "  ok   - $3"; else
        say "  FAIL - $3 (expected $2, got $1)"; fail=1; fi
}

configure() {
    rm -rf _tiertest
    if cmake -B _tiertest -DAST_C_STANDARD=11 > "$STASH/log" 2>&1; then
        echo pass
    else
        echo fail
    fi
}

say "case 1: git checkout, all fixtures present -> configure succeeds"
check "$(configure)" pass "baseline configures"

say "case 2: cutdown tarball, wcsconv/expected omitted -> configure succeeds"
mv "$FIX/wcsconv/expected" "$STASH/expected"
echo cutdown > "$FIX/DIST_TIER"
check "$(configure)" pass "deliberate omission is tolerated"

say "case 3: git checkout, wcsconv/expected missing -> configure FAILS"
rm -f "$FIX/DIST_TIER"
check "$(configure)" fail "missing fixture in a checkout is fatal"
mv "$STASH/expected" "$FIX/wcsconv/expected"

# A fixture that ast_require_fixture actually checks.  plot/expected/*.svg is
# manifest-listed, so its absence is a packaging defect rather than a deliberate
# omission, even inside a cut-down tarball.
say "case 4: cutdown tarball, a manifest-listed fixture missing -> FAILS"
mv "$FIX/plot/expected/aitoff.svg" "$STASH/aitoff.svg"
echo cutdown > "$FIX/DIST_TIER"
check "$(configure)" fail "missing manifest-listed fixture is fatal"
mv "$STASH/aitoff.svg" "$FIX/plot/expected/aitoff.svg"

exit $fail
