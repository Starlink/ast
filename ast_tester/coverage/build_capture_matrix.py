#!/usr/bin/env python3
"""Compute, per captured candidate, which of the target branches it covers.

Resumable: appends one JSONL record per candidate to MATRIX; on restart it
skips candidates already recorded. Designed to run detached in the
background so it survives interruptions. Run with the gcovr venv active and
from the repo root:

    source ~/pyenv/bin/activate
    python3 ast_tester/coverage/build_capture_matrix.py
"""
import glob
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.getcwd(), "ast_tester/coverage"))
import simplify_coverage_gaps as scg  # noqa: E402

REPO = os.getcwd()
ASTDIR = "build-cov/CMakeFiles/ast.dir"
SIMP = "build-cov/ast_tester/simplify"
GCOV = os.environ.get("GCOV", "xcrun llvm-cov gcov")
FILTER = (r"src/(mapping|box|cmpmap|dssmap|grismmap|interval|intramap|lutmap|"
          r"mathmap|matrixmap|normmap|nullregion|pcdmap|permmap|pointlist|"
          r"polymap|prism|ratemap|selectormap|shiftmap|slamap|specmap|sphmap|"
          r"splinemap|switchmap|timemap|tranmap|unitmap|unitnormmap|wcsmap|"
          r"winmap|xphmap|zoommap)\.c$")
CANDS = sorted(glob.glob("/tmp/cands/cand_*.map"))
TARGET = {tuple(k) for k in json.load(open("/tmp/target294.json"))}
MATRIX = "/tmp/matrix.jsonl"
TMPJSON = "/tmp/_cand_cov.json"


def reset_gcda():
    for x in glob.glob(f"{ASTDIR}/**/*.gcda", recursive=True):
        try:
            os.remove(x)
        except OSError:
            pass


def main():
    done = set()
    if os.path.exists(MATRIX):
        for line in open(MATRIX):
            try:
                done.add(json.loads(line)["cand"])
            except (ValueError, KeyError):
                pass
    out = open(MATRIX, "a")
    total = len(CANDS)
    for i, f in enumerate(CANDS):
        if f in done:
            continue
        reset_gcda()
        subprocess.run([SIMP, f, "/tmp/_o.simp"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(
            f'gcovr --root "{REPO}" "{ASTDIR}" --gcov-executable "{GCOV}" '
            f'--filter "{FILTER}" --json --output "{TMPJSON}"',
            shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        covers = []
        try:
            cov = scg.load_coverage(TMPJSON)
            covers = [list(k) for k in TARGET if cov["branches"].get(k, 0) > 0]
        except (OSError, ValueError):
            pass
        out.write(json.dumps({"cand": f, "covers": covers}) + "\n")
        out.flush()
        if i % 25 == 0:
            print(f"[{i}/{total}] {os.path.basename(f)} covers {len(covers)}",
                  flush=True)
    out.close()
    print("MATRIX_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
