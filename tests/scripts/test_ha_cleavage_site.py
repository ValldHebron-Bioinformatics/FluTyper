"""Tests for scripts/ha_cleavage_site.py.

The motifs used here are the ones FluTyper's own protein output gives for this protocol's
samples: an H5 clade 2.3.4.4b HPAI virus ends HA1 with ...SPLREKRRKR and an H5 LPAI virus
with ...LRNVPQRETR, both followed by HA2 starting GLFGA.
"""
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from ha_cleavage_site import read_cleavage  # noqa: E402

HPAI_HA1 = "MEKIVLLLAIVSPLREKRRKR"
LPAI_HA1 = "MEKIVLLLAIVLRNVPQRETR"
HA2 = "GLFGAIAGFIE"


class TestReadCleavage(unittest.TestCase):
    def test_polybasic_motif_is_hpai(self):
        out = read_cleavage(HPAI_HA1, HA2)
        self.assertEqual(out["call"], "HPAI")
        self.assertTrue(out["motif"].endswith("/GLF"), out["motif"])
        self.assertIn("REKRRKR", out["motif"])

    def test_monobasic_motif_is_lpai(self):
        out = read_cleavage(LPAI_HA1, HA2)
        self.assertEqual(out["call"], "LPAI")
        self.assertIn("PQRETR", out["motif"])

    def test_a_furin_consensus_with_three_basics_is_hpai(self):
        # R-X-K/R-R immediately before the junction is the furin site, even with fewer basics
        out = read_cleavage("MEKIVLLLAIVPQRETRGKR", HA2)
        self.assertEqual(out["call"], "HPAI")

    def test_missing_ha1_is_indeterminate_not_lpai(self):
        out = read_cleavage("", HA2)
        self.assertEqual(out["call"], "indeterminate")
        self.assertEqual(out["motif"], "")

    def test_missing_ha2_anchor_is_indeterminate(self):
        out = read_cleavage(HPAI_HA1, "")
        self.assertEqual(out["call"], "indeterminate")
        self.assertIn("anchor", out["reason"])

    def test_wrong_ha2_anchor_is_indeterminate(self):
        out = read_cleavage(HPAI_HA1, "MKTIIALSYI")
        self.assertEqual(out["call"], "indeterminate")
        self.assertIn("anchor", out["reason"])

    def test_an_unresolved_residue_in_the_motif_is_indeterminate(self):
        out = read_cleavage("MEKIVLLLAIVSPLREKRXKR", HA2)
        self.assertEqual(out["call"], "indeterminate")
        self.assertIn("X", out["reason"])

    def test_a_gap_in_the_motif_is_indeterminate(self):
        out = read_cleavage("MEKIVLLLAIVSPLREK-RKR", HA2)
        self.assertEqual(out["call"], "indeterminate")

    def test_trailing_stop_and_gaps_are_stripped_before_reading(self):
        out = read_cleavage(HPAI_HA1 + "*--", HA2)
        self.assertEqual(out["call"], "HPAI")

    def test_an_ha1_too_short_to_hold_a_motif_is_indeterminate(self):
        out = read_cleavage("KR", HA2)
        self.assertEqual(out["call"], "indeterminate")

    def test_the_basic_count_is_reported(self):
        self.assertEqual(read_cleavage(HPAI_HA1, HA2)["basic_p1_p5"], 5)
        self.assertEqual(read_cleavage(LPAI_HA1, HA2)["basic_p1_p5"], 2)

    def test_a_non_basic_p1_is_indeterminate(self):
        # the cleavage P1 residue is an arginine in every readable HA; anything else means the
        # junction was not where the alignment put it
        out = read_cleavage("MEKIVLLLAIVLRNVPQRETQ", HA2)
        self.assertEqual(out["call"], "indeterminate")


if __name__ == "__main__":
    unittest.main()
