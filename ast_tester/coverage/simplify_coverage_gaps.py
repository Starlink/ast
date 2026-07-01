#!/usr/bin/env python3
"""Diff two gcovr JSON coverage files into a simplify-pathway gap ledger.

Stdlib only. Reads gcovr 8.x JSON (the ``--json`` format): each file has
``lines`` (with ``line_number``, ``count``, ``branches:[{count}]``) and
``functions`` (with ``name``, ``lineno``). Classifies each uncovered branch
in a target function as ``differential`` (hit by the full suite but not by
simplify-only) or ``absolute-only`` (hit by neither).
"""
import argparse
import bisect
import json
import os
import re

# Functions in the astSimplify call graph. NOT included: Decompose and
# RemoveRegions -- those are reached only through the separate public APIs
# astDecompose and astRemoveRegions (e.g. from FitsChan/FrameSet/Region code),
# never from astSimplify, so a simplify fixture cannot exercise them. They are
# distinct efforts, out of scope for the simplify self-containment goal.
#
# The list covers both the per-class MapMerge/MapList/Simplify methods and the
# static merge/swap helper functions they call: the region-pair mergers
# (MergeBox/MergeInterval/MergeNullRegion/MergePointList), the PolyMap+ShiftMap
# mergers (MergeShift/MergeShifts), the UnitNormMap merger (MakeMergedMap,
# GetMappingType), the WinMap/MatrixMap swap executors (WinMat/MatWin/MatWin2/
# WinPerm/WinWcs) and the swap/merge feasibility helpers (CanSwap/CanMerge/
# PermGet). Without the helpers the ledger only saw the MapMerge method bodies
# and undercounted the merge engine.
ALLOWLIST = {"MapMerge", "MapList", "Simplify", "CombineMaps",
             "WinMat", "MatWin", "MatWin2", "WinPerm", "WinWcs",
             "MergeBox", "MergeInterval", "MergeNullRegion", "MergePointList",
             "MergeShift", "MergeShifts", "MakeMergedMap", "GetMappingType",
             "CanSwap", "CanMerge", "PermGet"}

# Region classes whose `Simplify` override performs geometric
# self-simplification rather than Mapping-pipeline merging. That is a separate
# effort with a different fixture methodology, so it is recorded as a deferred
# backlog rather than counted in the main merge-engine ledger. Their MapMerge
# (Region-as-Mapping) branches remain in scope.
REGION_SIMPLIFY_FILES = {"box.c", "interval.c", "prism.c",
                         "pointlist.c", "nullregion.c"}

INSCOPE = {
    "mapping.c", "box.c", "cmpmap.c", "dssmap.c", "grismmap.c", "interval.c",
    "intramap.c", "lutmap.c", "mathmap.c", "matrixmap.c", "normmap.c",
    "nullregion.c", "pcdmap.c", "permmap.c", "pointlist.c", "polymap.c",
    "prism.c", "ratemap.c", "selectormap.c", "shiftmap.c", "slamap.c",
    "specmap.c", "sphmap.c", "splinemap.c", "switchmap.c", "timemap.c",
    "tranmap.c", "unitmap.c", "unitnormmap.c", "wcsmap.c", "winmap.c",
    "xphmap.c", "zoommap.c",
}


def load_coverage_obj(obj):
    """Turn a parsed gcovr JSON object into a CovData dict.

    ``linefunc`` records gcovr's authoritative per-line ``function_name``
    when present; ``functions`` is the per-file function table used as a
    fallback to attribute a line to its enclosing function.
    """
    branches = {}
    functions = {}
    linefunc = {}
    for f in obj.get("files", []):
        path = f["file"]
        functions[path] = sorted(
            (fn["lineno"], fn["name"]) for fn in f.get("functions", []))
        for ln in f.get("lines", []):
            line = ln["line_number"]
            fname = ln.get("function_name")
            if fname:
                linefunc[(path, line)] = fname
            for idx, br in enumerate(ln.get("branches", [])):
                branches[(path, line, idx)] = br.get("count", 0)
    return {"branches": branches, "functions": functions, "linefunc": linefunc}


def load_coverage(path):
    with open(path) as fh:
        return load_coverage_obj(json.load(fh))


# A "pure" status guard is a line whose whole condition is just astOK, e.g.
# `if ( !astOK ) return;` or `if ( astOK ) {`. Its uncovered branch is the
# error direction, which never fires in a passing run, so it is filtered out
# of the gap ledger as defensive-code noise. Compound conditions such as
# `while ( astOK && cond )` are NOT filtered -- they carry real logic.
_GUARD = re.compile(r"^\s*(?:if|while)\s*\(\s*!?\s*astOK\s*\)")
_src_cache = {}


def _is_status_guard(file, line):
    lines = _src_cache.get(file)
    if lines is None:
        try:
            with open(file) as fh:
                lines = fh.read().split("\n")
        except OSError:
            lines = []
        _src_cache[file] = lines
    return 0 < line <= len(lines) and bool(_GUARD.match(lines[line - 1]))


def function_for(file, line, functions):
    table = functions.get(file, [])
    linenos = [t[0] for t in table]
    pos = bisect.bisect_right(linenos, line) - 1
    if pos < 0:
        return None
    return table[pos][1]


def classify_gaps(simp, full, allowlist=ALLOWLIST, inscope=INSCOPE,
                  filter_guards=True):
    gaps = []
    linefunc = full.get("linefunc", {})
    for key, full_count in full["branches"].items():
        file, line, idx = key
        if os.path.basename(file) not in inscope:
            continue
        fn = linefunc.get((file, line)) or function_for(file, line, full["functions"])
        if fn not in allowlist:
            continue
        if fn == "Simplify" and os.path.basename(file) in REGION_SIMPLIFY_FILES:
            continue  # geometric self-simplify; tracked in the deferred list
        if filter_guards and _is_status_guard(file, line):
            continue  # defensive astOK guard; error direction never fires
        simp_count = simp["branches"].get(key, 0)
        if simp_count > 0:
            continue  # already self-covered
        cls = "differential" if full_count > 0 else "absolute-only"
        gaps.append((file, line, idx, fn, cls))
    gaps.sort()
    return gaps


def region_simplify_gaps(simp, full):
    """Deferred backlog: uncovered branches in Region geometric `Simplify`."""
    gaps = []
    linefunc = full.get("linefunc", {})
    for key, full_count in full["branches"].items():
        file, line, idx = key
        if os.path.basename(file) not in REGION_SIMPLIFY_FILES:
            continue
        fn = linefunc.get((file, line)) or function_for(file, line, full["functions"])
        if fn != "Simplify":
            continue
        if simp["branches"].get(key, 0) > 0:
            continue
        cls = "differential" if full_count > 0 else "absolute-only"
        gaps.append((file, line, idx, fn, cls))
    gaps.sort()
    return gaps


_ROW = re.compile(
    r"\|\s*`([^:]+):(\d+)`\s*\|\s*(\d+)\s*\|\s*(\w+)\s*\|\s*[\w-]+\s*\|\s*([^|]*?)\s*\|")


def parse_prior(md_text):
    prior = {}
    for m in _ROW.finditer(md_text):
        file, line, idx, fn, status = m.groups()
        status = status.strip()
        if status and status != "open":
            prior[(file, fn, int(line), int(idx))] = status
    return prior


def render_ledger(gaps, prior, full, simp, deferred=(), n_guards=0):
    lines = ["# `astSimplify` Self-Contained Coverage Gap Ledger", "",
             "Regenerated by `run_simplify_coverage.sh`. Do not hand-edit",
             "rows except the Status column. See `README.md`.", "",
             f"_{n_guards} defensive astOK-guard error-direction branches "
             "filtered out as non-actionable noise._", ""]
    open_rows, closed_rows = [], []
    live_keys = set()
    for file, line, idx, fn, cls in gaps:
        key = (file, fn, line, idx)
        live_keys.add(key)
        status = prior.get(key, "open")
        open_rows.append(
            f"| `{file}:{line}` | {idx} | {fn} | {cls} | {status} |")
    # Closed: anything previously annotated that is no longer a gap.
    for (file, fn, line, idx), status in sorted(prior.items()):
        if (file, fn, line, idx) not in live_keys:
            closed_rows.append(
                f"| `{file}:{line}` | {idx} | {fn} | {status} |")
    diff_n = sum(1 for r in open_rows if " differential " in r)
    abs_n = sum(1 for r in open_rows if " absolute-only " in r)
    lines += [f"**Open differential gaps:** {diff_n}  ",
              f"**Open absolute-only gaps:** {abs_n}", "",
              "## Open gaps", "",
              "| Location | Branch | Function | Class | Status |",
              "| --- | --- | --- | --- | --- |"]
    lines += open_rows or ["| _none_ |  |  |  |  |"]
    lines += ["", "## Closed", "",
              "| Location | Branch | Function | Status |",
              "| --- | --- | --- | --- |"]
    lines += closed_rows or ["| _none_ |  |  |  |"]
    def_diff = sum(1 for gp in deferred if gp[4] == "differential")
    lines += ["", "## Deferred: Region geometric Simplify (out of scope)", "",
              "Self-simplification of Region geometry (Box/Interval/Prism/etc.).",
              "Tracked here as a known backlog; a separate effort with a",
              "different (geometric) fixture methodology — not part of the",
              "merge-engine self-containment goal.", "",
              f"**Deferred branches:** {len(deferred)} ({def_diff} differential)", "",
              "| Location | Branch | Function | Class |",
              "| --- | --- | --- | --- |"]
    lines += [f"| `{file}:{line}` | {idx} | {fn} | {cls} |"
              for file, line, idx, fn, cls in deferred] or ["| _none_ |  |  |  |"]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--simplify", required=True)
    ap.add_argument("--full", required=True)
    ap.add_argument("--ledger", required=True)
    args = ap.parse_args()
    simp = load_coverage(args.simplify)
    full = load_coverage(args.full)
    gaps = classify_gaps(simp, full)
    n_guards = len(classify_gaps(simp, full, filter_guards=False)) - len(gaps)
    deferred = region_simplify_gaps(simp, full)
    prior = {}
    if os.path.exists(args.ledger):
        with open(args.ledger) as fh:
            prior = parse_prior(fh.read())
    with open(args.ledger, "w") as fh:
        fh.write(render_ledger(gaps, prior, full, simp, deferred, n_guards))
    diff_n = sum(1 for gp in gaps if gp[4] == "differential")
    print(f"Wrote {args.ledger}: {len(gaps)} open gaps "
          f"({diff_n} differential); {len(deferred)} deferred region-Simplify; "
          f"{n_guards} astOK-guard branches filtered.")


if __name__ == "__main__":
    main()
