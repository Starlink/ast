#!/usr/bin/env python3
"""Numerically-tolerant file comparison.

Compares two text files.  Non-numeric content must match exactly;
floating-point values on the same line are compared with relative
and absolute tolerance.

Exits 0 if files are equivalent, 1 if they differ.

Usage: numdiff.py <reference> <got>
"""

import re
import sys
from difflib import SequenceMatcher

RTOL = 1e-10
ATOL = 1e-14

NUM = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?')


def close(a, b):
    return abs(a - b) <= ATOL + RTOL * max(abs(a), abs(b))


def near(a, b):
    if a == b:
        return True

    pa, pb = NUM.split(a), NUM.split(b)
    na, nb = NUM.findall(a), NUM.findall(b)

    if pa != pb or len(na) != len(nb):
        return False

    return all(close(float(x.replace('D', 'E').replace('d', 'e')),
                     float(y.replace('D', 'E').replace('d', 'e')))
               for x, y in zip(na, nb))


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: numdiff.py <reference> <got>\n")
        sys.exit(2)

    with open(sys.argv[1]) as f:
        want = f.readlines()
    with open(sys.argv[2]) as f:
        got = f.readlines()

    ok = True
    opcodes = SequenceMatcher(None, want, got).get_opcodes()
    for op, idx1, idx2, jdx1, jdx2 in opcodes:
        if op == 'equal':
            continue
        if (op == 'replace' and idx2 - idx1 == jdx2 - jdx1 and
                all(near(r.rstrip('\n'), g.rstrip('\n'))
                    for r, g in zip(want[idx1:idx2], got[jdx1:jdx2]))):
            continue
        for line in want[idx1:idx2]:
            print('-', line, end='')
        for line in got[jdx1:jdx2]:
            print('+', line, end='')
        ok = False

    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
