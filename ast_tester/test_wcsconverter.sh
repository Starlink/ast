#!/bin/sh
# Test WCS encoding conversions using wcsconverter.
# Each case reads a FrameSet from an input file and writes it in a specified
# encoding, then diffs the result against a reference file.
srcdir=${srcdir:-.}
ok=0

run() {
    base=$1 in_suf=$2 enc=$3
    attrs="${4:+$4,}FitsDigits=8"
    in=${srcdir}/${base}.${in_suf}
    ref=${srcdir}/${base}.${enc}
    out=${base}-new.${enc}
    ./wcsconverter "${in}" "${enc}" "${out}" "${attrs}" || { ok=1; return; }
    if diff -c "${ref}" "${out}" > "${out}.diff" 2>&1; then
        rm -f "${out}" "${out}.diff"
    else
        echo "FAIL: ${base}.${in_suf} -> ${enc}: differs from reference"
        cat "${out}.diff"
        rm -f "${out}" "${out}.diff"
        ok=1
    fi
}

run timj        ast      fits-wcs  "cdmatrix=1"
run timj        ast      fits-iraf
run timj        ast      fits-aips
run timj        ast      fits-pc   "fitsrounding=0"
run timj        ast      native
run timj        ast      native
run a20070718_00010_02_cube ast fits-wcs
run dss         fits-dss ast
run dss         ast      dss
run dss         ast      fits-wcs  "cdmatrix=1"
run degen1      ast      fits-wcs  "cdmatrix=1"
run degen1      ast      fits-wcs  "cdmatrix=1"
run sip         head     fits-wcs  "cdmatrix=1,sipreplace=0"
run sip2        head     fits-wcs  "cdmatrix=1,sipreplace=0"
run lsst        ast      fits-wcs  "cdmatrix=1"
run longslit    fits-pc  fits-wcs  "cdmatrix=1"

exit $ok
