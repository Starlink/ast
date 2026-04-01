#!/bin/sh
# Test Mapping simplification using simplify.
# Each .map file in the source directory is simplified and the result
# diffed against the corresponding .simp reference file.
srcdir=${srcdir:-.}
ok=0

python3="${PYTHON3:-python3}"
if ! "${python3}" -c '' 2>/dev/null; then
    echo "test_simplify: skipped (python3 not found)"
    exit 77
fi

for map in "${srcdir}"/*.map; do
    base=$(basename "${map}" .map)
    ref=${srcdir}/${base}.simp
    out=${base}-new.simp
    ./simplify "${map}" "${out}" || { ok=1; continue; }
    if "${python3}" "${srcdir}/numdiff.py" "${ref}" "${out}" > "${out}.diff" 2>&1; then
        rm -f "${out}" "${out}.diff"
    else
        echo "FAIL: ${base}.map: simplified output differs from reference"
        cat "${out}.diff"
        rm -f "${out}" "${out}.diff"
        ok=1
    fi
done

exit $ok
