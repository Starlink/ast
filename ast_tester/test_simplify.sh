#!/bin/sh
# Test Mapping simplification using simplify.
# Each .map file in the source directory is simplified and the result
# diffed against the corresponding .simp reference file.
srcdir=${srcdir:-.}
ok=0

for map in "${srcdir}"/*.map; do
    base=$(basename "${map}" .map)
    ref=${srcdir}/${base}.simp
    out=${base}-new.simp
    ./simplify "${map}" "${out}" || { ok=1; continue; }
    if diff -c "${ref}" "${out}" > "${out}.diff" 2>&1; then
        rm -f "${out}" "${out}.diff"
    else
        echo "FAIL: ${base}.map: simplified output differs from reference"
        cat "${out}.diff"
        rm -f "${out}" "${out}.diff"
        ok=1
    fi
done

exit $ok
