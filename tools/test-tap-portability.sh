#!/bin/sh
# Guard against a construct macOS's /bin/sh cannot parse.
#
# macOS ships bash 3.2, whose $( ) parser matches parentheses naively: the
# unbalanced ")" that ends a case pattern closes the substitution early, and the
# following ";;" is then a syntax error.  Linux shells lex it correctly, so the
# mistake passes every local test and fails only on macOS CI.
#
# The rule: no case statement inside a multi-line $( ) block.  Keep the loop at
# top level and collect its output in a file, as the other producers here do.
cd "$(dirname "$0")/.."
fail=0

for f in ast_tester/*.tap ast_tester/test_*.sh ast_tester/dist-fixtures.sh; do
    test -f "$f" || continue
    bad=$(awk '
        /\$\($/            { depth++; next }
        depth && /^[[:space:]]*\)[[:space:]]*$/ { depth--; next }
        depth && /(^|[[:space:]])case[[:space:]]/ { print FILENAME ":" FNR ": " $0 }
    ' "$f")
    if [ -n "$bad" ]; then
        echo "  FAIL - case inside a multi-line \$( ) block (breaks bash 3.2):"
        printf '%s\n' "$bad" | sed 's/^/         /'
        fail=1
    else
        echo "  ok   - $f"
    fi
done
exit $fail
