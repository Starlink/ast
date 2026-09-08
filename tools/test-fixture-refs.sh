#!/bin/sh
# Every fixture a test names by hand must be one the default distribution
# carries, or that test fails in the tarball while passing in a checkout.
#
# test_oracle_util.c loaded simplify/matrix_diagonal_to_zoom.map and so failed in
# the cut-down tarball; nothing local caught it because the path is relative to
# the fixture root rather than starting with "fixtures/", which is the form an
# earlier audit grepped for.  Both forms are checked here.
#
# Generators are exempt: they regenerate the corpus and only ever run in a full
# checkout.
cd "$(dirname "$0")/.."
AREAS='simplify|wcsconv|programs|serialisation|oracle|plot'
shipped=$(cd ast_tester && ./dist-fixtures.sh fixtures/DIST_MANIFEST) || exit 1
fail=0

for f in ast_tester/*.c ast_tester/*.f; do
    test -f "$f" || continue
    case ${f##*/} in gen_*) continue ;; esac
    refs=$(grep -ohE "\"(fixtures/)?($AREAS)/[A-Za-z0-9_./-]+\"" "$f" 2>/dev/null \
           | tr -d '"' | sed 's|^fixtures/||' | sort -u)
    for r in $refs; do
        case $r in */) continue ;; esac
        if ! printf '%s\n' "$shipped" | grep -qx "fixtures/$r"; then
            echo "  FAIL - $f names a fixture the default dist omits: $r"
            fail=1
        fi
    done
done
test "$fail" -eq 0 && echo "  ok   - every hand-named fixture is distributed"
exit $fail
