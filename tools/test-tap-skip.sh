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

# The shape a cut-down distribution actually has: the cases.txt manifests ship so
# that the suites can enumerate and report what they cannot run, but the fixtures
# those manifests name do not.  A producer must skip per row here, not try to run
# and fail 825 times.
CUT=$(mktemp -d "${TMPDIR:-/tmp}/ast_cut.XXXXXX")
mkdir -p "$CUT/fixtures/simplify"
cp "$SRC/fixtures/simplify/cases.txt" "$CUT/fixtures/simplify/"
( cd "$SRC/fixtures" && find programs -name '*.ast' -type f ) | while read -r f; do
    mkdir -p "$CUT/fixtures/$(dirname "$f")"
    cp "$SRC/fixtures/$f" "$CUT/fixtures/$f"
done

out=$(srcdir="$CUT" sh "$SRC/test_roundtrip.tap" 2>&1); rc=$?
check "$rc" 0 "roundtrip exits 0 when only the manifest shipped"
check "$(printf '%s\n' "$out" | grep -c '^not ok')" 0 \
      "roundtrip reports no failures for undistributed fixtures"
skips=$(printf '%s\n' "$out" | grep -c '# SKIP')
if [ "$skips" -gt 0 ]; then
    echo "  ok   - roundtrip reports $skips rows as skipped"
else
    echo "  FAIL - roundtrip reported no skips at all"; fail=1
fi
real=$(printf '%s\n' "$out" | grep '^ok' | grep -vc '# SKIP')
check "$real" 5 "roundtrip still runs the 5 shipped programs/*.ast fixtures"
rm -rf "$CUT"

rm -rf "$EMPTY"
exit $fail
