"""Tests for scripts/gisaid_to_metadata.py.

Geocoding is network I/O, so every test replaces gisaid_to_metadata._nominatim with a
stub: any place name resolves except the ones in NOT_FOUND, which stand for a location
GISAID spells in a way Nominatim cannot resolve.
"""
import contextlib
import io
import json
import sys
import urllib.error
import urllib.request
from unittest import mock
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

# Replacing a module function is global, so every test that does it restores this.
ORIGINAL_NOMINATIM = g2m._nominatim


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
        g2m._nominatim = ORIGINAL_NOMINATIM
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


class TestNominatimResilience(unittest.TestCase):
    """A 213-location conversion must survive a transient Nominatim hiccup.

    An empty result list means 'this place is not in Nominatim' and is a normal answer; an
    HTTP or network error means the service is unavailable and is not an answer at all. The
    two must not be collapsed, or an outage would silently blank every location.
    """

    @staticmethod
    def _ok_body():
        return io.BytesIO(b'[{"lat":"1.0","lon":"2.0","address":{"country_code":"no"}}]')

    def setUp(self):
        g2m._cache.clear()
        g2m._country_codes.clear()
        self.sleeps = []

    def test_retries_a_transient_503_then_succeeds(self):
        calls = []

        def flaky(req, timeout=None):
            calls.append(1)
            if len(calls) < 3:
                raise urllib.error.HTTPError(req.full_url, 503, "Service Unavailable", {}, None)
            return self._ok_body()

        with mock.patch.object(urllib.request, "urlopen", flaky), \
             mock.patch.object(g2m.time, "sleep", self.sleeps.append):
            hit = g2m._nominatim({"q": "Vardo"})
        self.assertEqual(len(calls), 3)
        self.assertEqual(hit["lat"], "1.0")

    def test_backoff_waits_longer_each_retry(self):
        def always_503(req, timeout=None):
            raise urllib.error.HTTPError(req.full_url, 503, "Service Unavailable", {}, None)

        with mock.patch.object(urllib.request, "urlopen", always_503), \
             mock.patch.object(g2m.time, "sleep", self.sleeps.append):
            with self.assertRaises(g2m.GeocodeUnavailable):
                g2m._nominatim({"q": "Vardo"})
        waits = [s for s in self.sleeps if s > 1.1]
        self.assertEqual(waits, sorted(waits))
        self.assertGreater(len(waits), 1)

    def test_a_persistent_outage_aborts_naming_the_status(self):
        def always_503(req, timeout=None):
            raise urllib.error.HTTPError(req.full_url, 503, "Service Unavailable", {}, None)

        with mock.patch.object(urllib.request, "urlopen", always_503), \
             mock.patch.object(g2m.time, "sleep", self.sleeps.append):
            with self.assertRaises(g2m.GeocodeUnavailable) as ctx:
                g2m._nominatim({"q": "Vardo"})
        self.assertIn("503", str(ctx.exception))

    def test_an_empty_result_is_not_found_not_an_outage(self):
        with mock.patch.object(urllib.request, "urlopen",
                               lambda req, timeout=None: io.BytesIO(b"[]")), \
             mock.patch.object(g2m.time, "sleep", self.sleeps.append):
            self.assertIsNone(g2m._nominatim({"q": "Nowhereville"}))


class TestGeocodeCache(unittest.TestCase):
    """Lookups are cached on disk: an interrupted run must not have to pay for them twice."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        tmp = Path(self.tmp.name)
        rows = pandas.read_csv(FIXTURE, dtype=str, keep_default_na=False)
        self.fasta = tmp / "sequences.fasta"
        with open(self.fasta, "w") as fh:
            for _, r in rows.iterrows():
                fh.write(f">{r['Isolate_Id']}|HA|{r['Collection_Date']}|"
                         f"{r['Isolate_Name']}|Some_Reference_Lab\nACGT\n")
        self.out_csv = tmp / "meta.csv"
        self.out_tsv = tmp / "coords.tsv"
        self.cache = tmp / "cache.json"
        g2m._nominatim = fake_nominatim
        g2m._cache.clear()
        g2m._country_codes.clear()

    def tearDown(self):
        g2m._nominatim = ORIGINAL_NOMINATIM
        self.tmp.cleanup()

    def _run(self, nominatim):
        g2m._nominatim = nominatim
        g2m._cache.clear()
        g2m._country_codes.clear()
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            g2m.main(str(FIXTURE), str(self.fasta), str(self.out_csv), str(self.out_tsv),
                     str(self.cache))
        return err.getvalue()

    def test_a_second_run_needs_no_lookups_at_all(self):
        self._run(fake_nominatim)
        self.assertTrue(self.cache.exists())

        def forbidden(params):
            raise AssertionError(f"looked up {params} again instead of using the cache")

        self._run(forbidden)
        self.assertEqual(len(pandas.read_csv(self.out_csv)), 12)

    def test_an_outage_still_leaves_the_lookups_it_had_already_resolved(self):
        resolved = []

        def dies_after_three(params):
            if len(resolved) >= 3:
                raise g2m.GeocodeUnavailable("HTTP Error 503: Service Unavailable")
            resolved.append(params)
            return fake_nominatim(params)

        with self.assertRaises(SystemExit):
            self._run(dies_after_three)
        self.assertTrue(self.cache.exists())
        self.assertGreater(len(json.loads(self.cache.read_text())["places"]), 0)


class TestOfflineMode(unittest.TestCase):
    """--no-geocode: build the metadata without Nominatim at all.

    PhyloGeography resolves country, division and coordinates offline from GeoNames using
    LOCATION and COUNTRY, so the Auspice map never needs Nominatim. Offline mode keeps the
    deepest level GISAID wrote and leaves the coordinates table empty rather than guessing.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        tmp = Path(self.tmp.name)
        rows = pandas.read_csv(FIXTURE, dtype=str, keep_default_na=False)
        self.fasta = tmp / "sequences.fasta"
        with open(self.fasta, "w") as fh:
            for _, r in rows.iterrows():
                fh.write(f">{r['Isolate_Id']}|HA|{r['Collection_Date']}|"
                         f"{r['Isolate_Name']}|Some_Reference_Lab\nACGT\n")
        self.out_csv = tmp / "meta.csv"
        self.out_tsv = tmp / "coords.tsv"

        def forbidden(params):
            raise AssertionError(f"offline mode must not query Nominatim (got {params})")

        g2m._nominatim = forbidden
        g2m._cache.clear()
        g2m._country_codes.clear()
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            g2m.main(str(FIXTURE), str(self.fasta), str(self.out_csv), str(self.out_tsv),
                     None, False)
        self.stderr = err.getvalue()
        self.meta = pandas.read_csv(self.out_csv, keep_default_na=False)

    def tearDown(self):
        g2m._nominatim = ORIGINAL_NOMINATIM
        self.tmp.cleanup()

    def test_every_row_survives_without_a_single_lookup(self):
        self.assertEqual(len(self.meta), 12)

    def test_location_is_the_deepest_level_gisaid_wrote(self):
        by_id = self.meta.set_index("ID")
        self.assertEqual(by_id.loc["EPI_ISL_20060001", "LOCATION"], "Vardo")
        self.assertEqual(by_id.loc["EPI_ISL_20060008", "LOCATION"], "Nowhereville")

    def test_a_location_with_no_city_level_keeps_its_country(self):
        by_id = self.meta.set_index("ID")
        self.assertEqual(by_id.loc["EPI_ISL_20060006", "LOCATION"], "Netherlands")
        self.assertEqual(by_id.loc["EPI_ISL_20060006", "COUNTRY"], "Netherlands")

    def test_country_is_the_country_level(self):
        by_id = self.meta.set_index("ID")
        self.assertEqual(by_id.loc["EPI_ISL_20060001", "COUNTRY"], "Norway")
        self.assertEqual(by_id.loc["EPI_ISL_20060009", "COUNTRY"], "Spain")

    def test_no_coordinate_row_is_invented(self):
        coords = pandas.read_csv(self.out_tsv, sep="\t")
        self.assertEqual(len(coords), 0)
        self.assertEqual(list(coords.columns),
                         ["Comunidad", "Provincia", "Población", "Latitud", "Longitud", "Altitud"])

    def test_it_says_coordinates_were_not_resolved(self):
        self.assertIn("no-geocode", self.stderr)
