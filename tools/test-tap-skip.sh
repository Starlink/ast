#!/bin/sh
# A .tap producer whose fixtures are absent must emit a valid skipped plan and
# exit 0; test_keymap_mapsz.sh must exit 77, the code automake reads as a skip.
cd "$(dirname "$0")/.."
SRC=$PWD/ast_tester
BUILD=$PWD/build-dev/ast_tester
EMPTY=$(mktemp -d "${TMPDIR:-/tmp}/ast_tapskip.XXXXXX")
mkdir -p "$EMPTY/fixtures"
fail=0

check() {
    if [ "$1" = "$2" ]; then echo "  ok   - $3"; else
        echo "  FAIL - $3 (expected '$2', got '$1')"; fail=1; fi
}

cd "$BUILD"
for t in test_simplify test_wcsconverter test_roundtrip; do
    out=$(srcdir="$EMPTY" sh "$SRC/$t.tap" 2>&1); rc=$?
    check "$rc" 0 "$t.tap exits 0 with no fixtures"
    check "$(printf '%s' "$out" | head -1)" "1..0 # SKIP fixtures not distributed" \
          "$t.tap emits a skipped plan"
done

srcdir="$EMPTY" sh "$SRC/test_keymap_mapsz.sh" >/dev/null 2>&1; rc=$?
check "$rc" 77 "test_keymap_mapsz.sh exits 77 with no fixture"

rm -rf "$EMPTY"
exit $fail
