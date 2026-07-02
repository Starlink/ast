#!/usr/bin/env python3
"""Greedy set-cover over the capture matrix.

Reads /tmp/matrix.jsonl (one record per candidate: {cand, covers}) and
selects a near-minimal subset of candidates whose union covers the target
branch set. Writes the chosen candidate paths to /tmp/selected.json.
"""
import json

MATRIX = "/tmp/matrix.jsonl"
TARGET = "/tmp/target294.json"
OUT = "/tmp/selected.json"


def main():
    target = {tuple(k) for k in json.load(open(TARGET))}
    cand_cov = {}
    for line in open(MATRIX):
        rec = json.loads(line)
        cov = {tuple(k) for k in rec["covers"]} & target
        if cov:
            cand_cov[rec["cand"]] = cov

    uncovered = set(target)
    selected = []
    while uncovered:
        best, best_gain = None, 0
        for cand, cov in cand_cov.items():
            gain = len(cov & uncovered)
            if gain > best_gain:
                best, best_gain = cand, gain
        if not best:
            break
        selected.append({"cand": best, "new": best_gain})
        uncovered -= cand_cov[best]

    json.dump({"selected": selected,
               "n_selected": len(selected),
               "covered": len(target) - len(uncovered),
               "target": len(target),
               "uncovered": len(uncovered)},
              open(OUT, "w"), indent=2)
    print(f"selected {len(selected)} candidates covering "
          f"{len(target) - len(uncovered)}/{len(target)} target branches "
          f"({len(uncovered)} uncoverable)")


if __name__ == "__main__":
    main()
