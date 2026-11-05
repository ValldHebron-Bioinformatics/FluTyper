"""Tests for scripts/gisaid_to_metadata.py.

Geocoding is network I/O, so every test replaces gisaid_to_metadata._nominatim with a
stub: any place name resolves except the ones in NOT_FOUND, which stand for a location
GISAID spells in a way Nominatim cannot resolve.
"""
import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

import pandas

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import gisaid_to_metadata as g2m  # noqa: E402  (needs the sys.path line above)
from gisaid_to_metadata import isolate_key, read_export  # noqa: E402

FIXTURE = REPO / "tests" / "data" / "gisaid_avian_sample.csv"

# Every level of "Europe / Atlantis / Nowhereland / Nowhereville": no level geocodes, so the
# row exercises the fall-through where the script has nothing to put in LOCATION.
NOT_FOUND = {"Atlantis", "Nowhereland", "Nowhereville"}


def fake_nominatim(params):
    name = params["q"].split(",")[0].strip() if "q" in params else params.get("country", "")
    if name in NOT_FOUND:
        return None
    return {"address": {"country_code": name[:2].lower()}, "lat": "50.0", "lon": "4.0"}


class TestIsolateKey(unittest.TestCase):
    def test_accepts_both_spellings(self):
        self.assertEqual(isolate_key("EPI_ISL_20298047"), "EPI_ISL_20298047")
        self.assertEqual(isolate_key("20298047"), "EPI_ISL_20298047")


class TestReadExport(unittest.TestCase):
    def test_reads_a_csv_export(self):
        self.assertEqual(len(read_export(str(FIXTURE))), 12)


class TestMain(unittest.TestCase):
    """main() over the fixture, with the FASTA headers carrying the EPI_ISL_ prefix."""

    @classmethod
    def setUpClass(cls):
        g2m._nominatim = fake_nominatim
        g2m._cache.clear()
        g2m._country_codes.clear()

        cls.tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls.tmp.name)
        rows = pandas.read_csv(FIXTURE, dtype=str, keep_default_na=False)
        fasta = tmp / "sequences.fasta"
        with open(fasta, "w") as fh:
            for _, r in rows.iterrows():
                fh.write(f">{r['Isolate_Id']}|HA|{r['Collection_Date']}|"
                         f"{r['Isolate_Name']}|Some_Reference_Lab\nACGT\n")

        cls.out_csv = tmp / "meta.csv"
        cls.out_tsv = tmp / "coords.tsv"
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            g2m.main(str(FIXTURE), str(fasta), str(cls.out_csv), str(cls.out_tsv))
        cls.stderr = err.getvalue()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_keeps_rows_whose_fasta_ids_carry_the_epi_isl_prefix(self):
        meta = pandas.read_csv(self.out_csv)
        self.assertEqual(len(meta), 12)
        self.assertEqual(list(meta.columns), ["ID", "DATE", "LOCATION", "COUNTRY", "HOST",
                                              "SUBTYPE", "PATHOGENICITY", "CLADE", "GENOTYPE",
                                              "AGE GROUP", "SEX", "ORIGINATING_LAB"])

    def test_blank_pathogenicity_becomes_empty_not_dropped(self):
        meta = pandas.read_csv(self.out_csv, keep_default_na=False)
        self.assertEqual(len(meta[meta["PATHOGENICITY"] == ""]), 1)

    def test_ungeocodable_location_is_kept_as_written_and_warned(self):
        meta = pandas.read_csv(self.out_csv)
        self.assertIn("Nowhereville", list(meta["LOCATION"]))
        self.assertIn("could not geocode", self.stderr)

    def test_country_is_emitted_for_a_location_with_no_city_level(self):
        meta = pandas.read_csv(self.out_csv).set_index("ID")
        self.assertEqual(meta.loc["EPI_ISL_20060006", "COUNTRY"], "Netherlands")


if __name__ == "__main__":
    unittest.main()
