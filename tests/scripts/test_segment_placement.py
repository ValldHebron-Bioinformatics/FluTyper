"""Tests for scripts/segment_placement.py."""
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from segment_placement import placement  # noqa: E402


def tree_file(tmpdir, newick, name="tree.nwk"):
    path = Path(tmpdir) / name
    path.write_text(newick)
    return str(path)


class TestPlacement(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_focus_inside_an_hpai_clade_is_reported_with_its_support(self):
        tree = tree_file(self.tmp.name, "((focus,HPAI_a)95,(LPAI_b,LPAI_c)80);")
        out = placement(tree, "focus",
                        {"HPAI_a": "HPAI", "LPAI_b": "LPAI", "LPAI_c": "LPAI"})
        self.assertEqual(out["sister_pathotypes"], {"HPAI": 1})
        self.assertEqual(out["hpai_node_support"], 95.0)

    def test_focus_among_lpai_only_reports_no_hpai_node(self):
        tree = tree_file(self.tmp.name, "((focus,LPAI_b)90,LPAI_c);")
        out = placement(tree, "focus", {"LPAI_b": "LPAI", "LPAI_c": "LPAI"})
        self.assertIsNone(out["hpai_node_support"])
        self.assertEqual(out["sister_pathotypes"], {"LPAI": 1})

    def test_the_hpai_node_is_the_smallest_clade_joining_the_focus_to_any_hpai_tip(self):
        # focus sits with two LPAI tips; the HPAI tips join only at the outer node (70)
        tree = tree_file(self.tmp.name, "(((focus,LPAI_b)99,LPAI_c)88,(HPAI_a,HPAI_d)95)70;")
        out = placement(tree, "focus",
                        {"LPAI_b": "LPAI", "LPAI_c": "LPAI", "HPAI_a": "HPAI", "HPAI_d": "HPAI"})
        self.assertEqual(out["hpai_node_support"], 70.0)
        self.assertEqual(out["sister_pathotypes"], {"LPAI": 1})

    def test_sister_clade_tips_are_listed(self):
        tree = tree_file(self.tmp.name, "((focus,(HPAI_a,HPAI_d)99)95,LPAI_c);")
        out = placement(tree, "focus",
                        {"HPAI_a": "HPAI", "HPAI_d": "HPAI", "LPAI_c": "LPAI"})
        self.assertEqual(sorted(out["sister_tips"]), ["HPAI_a", "HPAI_d"])
        self.assertEqual(out["sister_pathotypes"], {"HPAI": 2})

    def test_a_tip_with_no_pathotype_is_counted_as_unknown_not_assumed(self):
        tree = tree_file(self.tmp.name, "((focus,MYSTERY)95,LPAI_c);")
        out = placement(tree, "focus", {"LPAI_c": "LPAI"})
        self.assertEqual(out["sister_pathotypes"], {"UNKNOWN": 1})
        self.assertIsNone(out["hpai_node_support"])

    def test_a_missing_focus_tip_is_reported_not_guessed(self):
        tree = tree_file(self.tmp.name, "((LPAI_b,HPAI_a)95,LPAI_c);")
        out = placement(tree, "focus", {"LPAI_b": "LPAI", "HPAI_a": "HPAI"})
        self.assertFalse(out["focus_present"])
        self.assertIsNone(out["hpai_node_support"])

    def test_support_absent_from_the_tree_is_none_not_zero(self):
        tree = tree_file(self.tmp.name, "((focus,HPAI_a),(LPAI_b,LPAI_c));")
        out = placement(tree, "focus",
                        {"HPAI_a": "HPAI", "LPAI_b": "LPAI", "LPAI_c": "LPAI"})
        self.assertEqual(out["sister_pathotypes"], {"HPAI": 1})
        self.assertIsNone(out["hpai_node_support"])

    def test_genotype_composition_of_the_sister_clade_is_reported(self):
        tree = tree_file(self.tmp.name, "((focus,(HPAI_a,HPAI_d)99)95,LPAI_c);")
        out = placement(tree, "focus",
                        {"HPAI_a": "HPAI", "HPAI_d": "HPAI", "LPAI_c": "LPAI"},
                        genotype_by_tip={"HPAI_a": "C.1", "HPAI_d": "B.24"})
        self.assertEqual(out["sister_genotypes"], {"B.24": 1, "C.1": 1})


if __name__ == "__main__":
    unittest.main()
