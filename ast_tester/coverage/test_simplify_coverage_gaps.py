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


if __name__ == "__main__":
    unittest.main()
