import unittest

import simplify_coverage_gaps as g


# Minimal gcovr-8.x JSON: one file, two functions, branch lines.
def cov(branch_counts):
    # branch_counts: {line: [counts...]}
    lines = []
    for line, counts in branch_counts.items():
        lines.append({"line_number": line, "count": sum(counts) or 0,
                      "branches": [{"count": c} for c in counts]})
    return {"files": [{"file": "src/winmap.c",
                       "lines": lines,
                       "functions": [{"name": "MapMerge", "lineno": 900},
                                     {"name": "WinMat", "lineno": 2793}]}]}


ALLOW = {"MapMerge", "WinMat"}
INSCOPE = {"winmap.c"}


class ClassifyTest(unittest.TestCase):
    def test_differential_and_absolute(self):
        simp = g.load_coverage_obj(cov({905: [5, 0], 910: [0, 0]}))
        full = g.load_coverage_obj(cov({905: [5, 3], 910: [0, 0]}))
        gaps = g.classify_gaps(simp, full, ALLOW, INSCOPE)
        keys = {(gp[1], gp[2], gp[4]) for gp in gaps}  # (line, idx, cls)
        # line 905 branch idx 1: full-hit, simp-miss -> differential
        self.assertIn((905, 1, "differential"), keys)
        # line 910 both branches hit by nobody -> absolute-only
        self.assertIn((910, 0, "absolute-only"), keys)
        self.assertIn((910, 1, "absolute-only"), keys)
        # line 905 branch idx 0: simp-hit -> not a gap
        self.assertNotIn((905, 0, "differential"), keys)
        self.assertNotIn((905, 0, "absolute-only"), keys)

    def test_function_for(self):
        functions = {"src/winmap.c": [(900, "MapMerge"), (2793, "WinMat")]}
        self.assertEqual(g.function_for("src/winmap.c", 905, functions), "MapMerge")
        self.assertEqual(g.function_for("src/winmap.c", 2800, functions), "WinMat")
        self.assertIsNone(g.function_for("src/winmap.c", 10, functions))

    def test_prior_annotation_preserved(self):
        prior = g.parse_prior(
            "| `src/winmap.c:905` | 1 | MapMerge | differential | fixture=win_x |\n")
        self.assertEqual(prior[("src/winmap.c", "MapMerge", 905, 1)], "fixture=win_x")


def _cov(file, fname, branch_counts):
    lines = [{"line_number": l, "function_name": fname, "count": sum(c) or 0,
              "branches": [{"count": x} for x in c]}
             for l, c in branch_counts.items()]
    return {"files": [{"file": file, "lines": lines,
                       "functions": [{"name": fname, "lineno": min(branch_counts)}]}]}


class RegionSimplifyTest(unittest.TestCase):
    def test_region_simplify_excluded_from_main(self):
        simp = g.load_coverage_obj(_cov("src/interval.c", "Simplify", {100: [0]}))
        full = g.load_coverage_obj(_cov("src/interval.c", "Simplify", {100: [5]}))
        # Region-file Simplify must not appear in the main merge-engine ledger.
        self.assertEqual(g.classify_gaps(simp, full, {"Simplify"}, {"interval.c"}), [])
        deferred = g.region_simplify_gaps(simp, full)
        self.assertEqual([(d[0], d[1], d[4]) for d in deferred],
                         [("src/interval.c", 100, "differential")])

    def test_cmpmap_simplify_stays_in_main(self):
        simp = g.load_coverage_obj(_cov("src/cmpmap.c", "Simplify", {100: [0]}))
        full = g.load_coverage_obj(_cov("src/cmpmap.c", "Simplify", {100: [5]}))
        main = g.classify_gaps(simp, full, {"Simplify"}, {"cmpmap.c"})
        self.assertEqual([(m[0], m[4]) for m in main], [("src/cmpmap.c", "differential")])
        # ...and it is NOT double-counted in the deferred region list.
        self.assertEqual(g.region_simplify_gaps(simp, full), [])


if __name__ == "__main__":
    unittest.main()
