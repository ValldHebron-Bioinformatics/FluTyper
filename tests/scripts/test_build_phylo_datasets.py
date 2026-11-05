"""Tests for scripts/build_phylo_datasets.py."""
import sys
import tempfile
import unittest
from pathlib import Path

import pandas

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from build_phylo_datasets import build_all, classify  # noqa: E402

META_HEADER = ("ID,DATE,LOCATION,COUNTRY,HOST,SUBTYPE,PATHOGENICITY,CLADE,GENOTYPE,"
               "AGE GROUP,SEX,ORIGINATING_LAB")

# id, pathogenicity, clade — the call column lives in the pathotype table below
CONTEXT = [
    ("EPI_ISL_20060001", "HPAI", "2.3.4.4b"),
    ("EPI_ISL_20060002", "HPAI", "2.3.4.4b"),
    ("EPI_ISL_20060003", "LPAI", ""),
    ("EPI_ISL_20060004", "LPAI", ""),
    ("EPI_ISL_20060005", "", ""),          # blank field, resolved by the call
    ("EPI_ISL_20060006", "", ""),          # blank field, no call -> UNKNOWN
    ("EPI_ISL_20060007", "HPAI", "2.3.4.4b"),
    ("EPI_ISL_20060008", "HPAI", ""),      # field says HPAI, the call says LPAI
]
# the cleavage read per sample, as scripts/ha_cleavage_site.py writes it
CALLS = {"EPI_ISL_20060005": "LPAI", "EPI_ISL_20060008": "LPAI", "25CR6244": "LPAI"}
FOCUS = "25CR6244"
FOCUS_SEGMENTS = ["HA", "PB2", "PB1", "NP", "MP", "NS"]


def fasta_headers(path):
    with open(path) as fh:
        return [line[1:].strip() for line in fh if line.startswith(">")]


def write_fixture(tmp, context):
    """Write the four inputs build_all() reads into tmp; returns their paths."""
    meta = tmp / "gisaid_meta.csv"
    with open(meta, "w") as fh:
        fh.write(META_HEADER + "\n")
        for iid, patho, clade in context:
            fh.write(f"{iid},2025-07-01,Malmo,Sweden,Mallard,A / H5N3,{patho},{clade},,,,Some Lab\n")

    patho_csv = tmp / "ha_cleavage.csv"
    with open(patho_csv, "w") as fh:
        fh.write("ID,cleavage_motif,pathotype_cleavage,basic_p1_p5,cleavage_reason\n")
        for sid, call in CALLS.items():
            fh.write(f"{sid},PQRETR/GLF,{call},2,monobasic cleavage site\n")

    gisaid_fasta = tmp / "gisaid.fasta"
    with open(gisaid_fasta, "w") as fh:
        for iid, _, _ in context:
            for seg in ("HA", "NP"):
                fh.write(f">{iid}|{seg}|2025-07-01|A/Mallard/Sweden/{iid}/2025|Some_Lab\nACGT\n")

    own_fasta = tmp / "own.fasta"
    with open(own_fasta, "w") as fh:
        for seg in FOCUS_SEGMENTS:
            fh.write(f">{FOCUS}_{seg}_H5N3_Avian_A/Mallard/Netherlands/OM_17/2018_EPI_ISL_19613905\nACGT\n")

    own_meta = tmp / "own_meta.csv"
    with open(own_meta, "w") as fh:
        fh.write("ID, DATE, LOCATION\n25CR6244, 2025-11-26, Barcelona\n")

    return meta, patho_csv, gisaid_fasta, own_fasta, own_meta


class TestClassify(unittest.TestCase):
    """classify() follows the HA cleavage site read, which is a reading of the sequence."""

    def test_the_cleavage_read_wins_over_the_gisaid_field(self):
        row = {"ID": "EPI_ISL_1", "PATHOGENICITY": "HPAI", "CLADE": "2.3.4.4b"}
        self.assertEqual(classify(row, {"EPI_ISL_1": "LPAI"}), "LPAI")

    def test_blank_gisaid_field_is_resolved_by_the_cleavage_read(self):
        row = {"ID": "EPI_ISL_2", "PATHOGENICITY": "", "CLADE": ""}
        self.assertEqual(classify(row, {"EPI_ISL_2": "LPAI"}), "LPAI")

    def test_no_read_and_no_field_is_unknown_never_assumed(self):
        row = {"ID": "EPI_ISL_3", "PATHOGENICITY": "", "CLADE": ""}
        self.assertEqual(classify(row, {}), "UNKNOWN")

    def test_an_indeterminate_read_falls_back_to_the_gisaid_field(self):
        row = {"ID": "EPI_ISL_4", "PATHOGENICITY": "HPAI", "CLADE": "2.3.4.4b"}
        self.assertEqual(classify(row, {"EPI_ISL_4": "indeterminate"}), "HPAI")

    def test_an_indeterminate_read_with_no_field_is_unknown_not_lpai(self):
        row = {"ID": "EPI_ISL_5", "PATHOGENICITY": "", "CLADE": ""}
        self.assertEqual(classify(row, {"EPI_ISL_5": "indeterminate"}), "UNKNOWN")

    def test_a_blank_read_does_not_override_the_gisaid_field(self):
        row = {"ID": "EPI_ISL_6", "PATHOGENICITY": "HPAI", "CLADE": "2.3.4.4b"}
        self.assertEqual(classify(row, {"EPI_ISL_6": ""}), "HPAI")


class TestBuild(unittest.TestCase):
    """build_all() over a fixture where no isolate is both LPAI and clade 2.3.4.4b."""

    n_clade_2344b = 3

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls.tmp.name)
        cls.outdir = tmp / "out"
        paths = write_fixture(tmp, CONTEXT)
        cls.counts = build_all(*[str(p) for p in paths], FOCUS, str(cls.outdir))
        cls.report = (cls.outdir / "dataset_report.txt").read_text()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_focus_sample_is_in_every_dataset(self):
        for ds in ("lpai", "hpai", "combined"):
            ids = {h.split("_")[0] for h in fasta_headers(self.outdir / f"{ds}.fasta")}
            self.assertIn("25CR6244", ids)

    def test_focus_keeps_all_six_of_its_segments_in_every_dataset(self):
        for ds in ("lpai", "hpai", "combined"):
            segs = [h.split("_")[1] for h in fasta_headers(self.outdir / f"{ds}.fasta")
                    if h.startswith("25CR6244_")]
            self.assertEqual(sorted(segs), sorted(FOCUS_SEGMENTS))

    def test_hpai_dataset_is_every_2344b_isolate(self):
        self.assertEqual(self.counts["hpai"]["context"], self.n_clade_2344b)

    def test_lpai_dataset_uses_the_call_over_the_gisaid_field(self):
        ids = {h.split("|")[0] for h in fasta_headers(self.outdir / "lpai.fasta")
               if h.startswith("EPI_ISL_")}
        self.assertEqual(ids, {"EPI_ISL_20060003", "EPI_ISL_20060004",
                               "EPI_ISL_20060005", "EPI_ISL_20060008"})

    def test_combined_is_the_union_of_both_contexts(self):
        self.assertEqual(self.counts["combined"]["context"],
                         self.counts["lpai"]["context"] + self.counts["hpai"]["context"])

    def test_unknown_isolates_are_listed_not_silently_dropped(self):
        self.assertIn("UNKNOWN", self.report)
        self.assertIn("EPI_ISL_20060006", self.report)

    def test_unknown_isolate_is_in_no_dataset(self):
        for ds in ("lpai", "hpai", "combined"):
            headers = " ".join(fasta_headers(self.outdir / f"{ds}.fasta"))
            self.assertNotIn("EPI_ISL_20060006", headers)

    def test_focus_metadata_row_carries_its_own_called_pathotype(self):
        meta = pandas.read_csv(self.outdir / "lpai_metadata.csv", keep_default_na=False)
        row = meta[meta["ID"] == FOCUS].iloc[0]
        self.assertEqual(row["PATHOGENICITY"], "LPAI")
        self.assertEqual(row["DATE"], "2025-11-26")
        self.assertEqual(row["LOCATION"], "Barcelona")

    def test_every_dataset_metadata_has_the_pipeline_header(self):
        for ds in ("lpai", "hpai", "combined"):
            meta = pandas.read_csv(self.outdir / f"{ds}_metadata.csv")
            self.assertEqual(list(meta.columns), META_HEADER.split(","))

    def test_metadata_and_fasta_cover_the_same_samples(self):
        for ds in ("lpai", "hpai", "combined"):
            meta_ids = set(pandas.read_csv(self.outdir / f"{ds}_metadata.csv")["ID"])
            fasta_ids = {h.split("|")[0] if h.startswith("EPI_ISL_") else h.split("_")[0]
                         for h in fasta_headers(self.outdir / f"{ds}.fasta")}
            self.assertEqual(meta_ids, fasta_ids)


class TestBuildWithOverlap(unittest.TestCase):
    """An isolate that is both LPAI and clade 2.3.4.4b belongs to both contexts."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls.tmp.name)
        cls.outdir = tmp / "out"
        context = CONTEXT + [("EPI_ISL_20060009", "", "2.3.4.4b")]
        calls = dict(CALLS, **{"EPI_ISL_20060009": "LPAI"})
        paths = list(write_fixture(tmp, context))
        with open(paths[1], "w") as fh:   # rewrite the cleavage table with the extra read
            fh.write("ID,cleavage_motif,pathotype_cleavage,basic_p1_p5,cleavage_reason\n")
            for sid, call in calls.items():
                fh.write(f"{sid},PQRETR/GLF,{call},2,monobasic cleavage site\n")
        cls.counts = build_all(*[str(p) for p in paths], FOCUS, str(cls.outdir))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_overlapping_isolate_is_in_both_contexts(self):
        self.assertEqual(self.counts["lpai"]["context"], 5)
        self.assertEqual(self.counts["hpai"]["context"], 4)

    def test_combined_counts_it_once(self):
        self.assertEqual(self.counts["combined"]["context"], 8)
        self.assertEqual(self.counts["overlap"], 1)


if __name__ == "__main__":
    unittest.main()
