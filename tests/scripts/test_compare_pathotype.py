"""Tests for scripts/compare_pathotype.py."""
import sys
import unittest
from pathlib import Path

import pandas

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from compare_pathotype import compare  # noqa: E402

# The cleavage read is the sequence-derived call and the one agreement is measured against.
CLEAVAGE = pandas.DataFrame([
    {"ID": "EPI_ISL_1", "cleavage_motif": "PLREKRRKR/GLF", "pathotype_cleavage": "HPAI"},
    {"ID": "EPI_ISL_2", "cleavage_motif": "PQRETR/GLF", "pathotype_cleavage": "LPAI"},
    {"ID": "EPI_ISL_3", "cleavage_motif": "PQRETR/GLF", "pathotype_cleavage": "LPAI"},
    {"ID": "EPI_ISL_4", "cleavage_motif": "", "pathotype_cleavage": "indeterminate"},
    {"ID": "EPI_ISL_5", "cleavage_motif": "", "pathotype_cleavage": "indeterminate"},
])

NEXTCLADE = pandas.DataFrame([
    {"Sample_ID": "EPI_ISL_1", "inferred_subtype": "H5N1", "pathotype": "HPAI"},
    {"Sample_ID": "EPI_ISL_2", "inferred_subtype": "H5N1", "pathotype": "HPAI"},
    {"Sample_ID": "EPI_ISL_3", "inferred_subtype": "H5N3", "pathotype": "LPAI"},
    {"Sample_ID": "EPI_ISL_5", "inferred_subtype": "H13N2", "pathotype": ""},
])

GISAID = pandas.DataFrame([
    {"ID": "EPI_ISL_1", "PATHOGENICITY": "HPAI", "CLADE": "2.3.4.4b", "GENOTYPE": "C.1"},
    {"ID": "EPI_ISL_2", "PATHOGENICITY": "HPAI", "CLADE": "2.3.4.4b", "GENOTYPE": "C.1"},
    {"ID": "EPI_ISL_3", "PATHOGENICITY": "", "CLADE": "", "GENOTYPE": ""},
    {"ID": "EPI_ISL_4", "PATHOGENICITY": "LPAI", "CLADE": "", "GENOTYPE": ""},
    {"ID": "EPI_ISL_5", "PATHOGENICITY": "", "CLADE": "", "GENOTYPE": ""},
])


class TestCompare(unittest.TestCase):
    def setUp(self):
        self.out = compare(CLEAVAGE, NEXTCLADE, GISAID)
        self.by_id = self.out.set_index("ID")

    def test_columns_keep_the_three_pathotype_sources_apart(self):
        self.assertEqual(list(self.out.columns),
                         ["ID", "inferred_subtype", "cleavage_motif", "pathotype_cleavage",
                          "pathotype_nextclade", "pathotype_gisaid", "clade", "genotype",
                          "agreement"])

    def test_agreement_is_measured_against_the_cleavage_read(self):
        self.assertEqual(self.by_id["agreement"].to_dict(), {
            "EPI_ISL_1": "agree",          # cleavage HPAI, GISAID HPAI
            "EPI_ISL_2": "disagree",       # cleavage LPAI, GISAID HPAI
            "EPI_ISL_3": "only_cleavage",  # cleavage LPAI, GISAID blank
            "EPI_ISL_4": "only_gisaid",    # cleavage indeterminate, GISAID LPAI
            "EPI_ISL_5": "neither",        # cleavage indeterminate, GISAID blank
        })

    def test_an_indeterminate_cleavage_is_not_a_call(self):
        self.assertEqual(self.by_id.loc["EPI_ISL_4", "pathotype_cleavage"], "indeterminate")
        self.assertEqual(self.by_id.loc["EPI_ISL_4", "agreement"], "only_gisaid")

    def test_the_nextclade_label_is_carried_but_never_drives_agreement(self):
        # EPI_ISL_2: nextclade says HPAI, the sequence says LPAI; agreement follows the sequence
        self.assertEqual(self.by_id.loc["EPI_ISL_2", "pathotype_nextclade"], "HPAI")
        self.assertEqual(self.by_id.loc["EPI_ISL_2", "pathotype_cleavage"], "LPAI")
        self.assertEqual(self.by_id.loc["EPI_ISL_2", "agreement"], "disagree")

    def test_the_motif_is_carried_through(self):
        self.assertEqual(self.by_id.loc["EPI_ISL_1", "cleavage_motif"], "PLREKRRKR/GLF")

    def test_covers_the_union_of_the_tables(self):
        self.assertEqual(len(self.out), 5)

    def test_clade_and_genotype_come_from_the_gisaid_table(self):
        self.assertEqual(self.by_id.loc["EPI_ISL_1", "clade"], "2.3.4.4b")
        self.assertEqual(self.by_id.loc["EPI_ISL_1", "genotype"], "C.1")

    def test_a_sample_only_the_cleavage_table_knows_is_kept(self):
        extra = pandas.DataFrame([
            {"ID": "25CR6244", "cleavage_motif": "RNVPQRETR/GLF", "pathotype_cleavage": "LPAI"}])
        out = compare(pandas.concat([CLEAVAGE, extra]), NEXTCLADE, GISAID)
        row = out.set_index("ID").loc["25CR6244"]
        self.assertEqual(row["pathotype_cleavage"], "LPAI")
        self.assertEqual(row["agreement"], "only_cleavage")


if __name__ == "__main__":
    unittest.main()
