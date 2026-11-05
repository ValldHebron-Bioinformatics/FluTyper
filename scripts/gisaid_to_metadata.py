#!/usr/bin/env python3
"""Convert a GISAID EpiFlu isolates export (.xls, .xlsx or .csv) into FluTyper inputs:

  out.csv  --metadata     ID,DATE,LOCATION,COUNTRY,HOST,SUBTYPE,PATHOGENICITY,CLADE,
                          GENOTYPE,AGE GROUP,SEX,ORIGINATING_LAB
  out.tsv  --coordinates  same columns as RESOURCES/coordenades_cat.tsv, covering
                          every location, region and lab of the metadata

HOST, SUBTYPE, PATHOGENICITY, CLADE and GENOTYPE are copied from the export when it has
them (an avian download does, a human one may not) and left empty otherwise. PATHOGENICITY
is GISAID's own call: FluTyper's pathotype is computed from the sequence, and the two are
compared downstream rather than merged here.

ORIGINATING_LAB is not in the GISAID table, so it is taken from the last field
of the FASTA headers (>ID|SEGMENT|DATE|NAME|LAB) of the same download.
Places are geocoded with OpenStreetMap Nominatim (place names only, 1 query/s).

Every Nominatim lookup is cached on disk (out.csv.geocache.json by default), so a run that
is interrupted - or hits a Nominatim outage - does not pay for the same lookups twice.

Usage: gisaid_to_metadata.py gisaid_epiflu_isolates.xls sequences.fasta out.csv out.tsv \
                             [geocache.json]
"""
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

# GISAID spellings that do not geocode or would split one region in two
FIXES = {"Basqsue Country": "Basque Country", "Catalunya": "Catalonia",
         # Greek prefectures come in the genitive with a "Nomos" prefix
         "Nomos Attica": "Attica", "Nomos Ioanninon": "Ioannina", "Nomos Lasithiou": "Lasithi"}

_cache = {}
_country_codes = {}

# Nominatim statuses worth another try: rate limiting and the transient server-side ones
RETRY_STATUS = {429, 500, 502, 503, 504}
MAX_ATTEMPTS = 5


class GeocodeUnavailable(RuntimeError):
    """Nominatim could not be reached or kept failing: no answer, as opposed to 'not found'."""

# Columns copied straight from the export when present: metadata column -> export column
PASSTHROUGH = {"HOST": "Host", "SUBTYPE": "Subtype", "PATHOGENICITY": "Pathogenicity",
               "CLADE": "Clade", "GENOTYPE": "Genotype"}


def isolate_key(raw_id):
    """'20298047' and 'EPI_ISL_20298047' are the same isolate: key both as the prefixed form.

    GISAID writes the bare number in the isolates table but the prefixed form in the FASTA
    headers of some downloads, so the two sides only join once both go through here."""
    text = str(raw_id).strip()
    return text if text.startswith("EPI_ISL_") else f"EPI_ISL_{text}"


def read_export(path):
    """Read a GISAID isolates export: .csv as text, anything else through pandas' Excel readers."""
    if str(path).lower().endswith(".csv"):
        return pd.read_csv(path, dtype={"Isolate_Id": str})
    return pd.read_excel(path, dtype={"Isolate_Id": str})


def _nominatim(params):
    """One Nominatim query, retried with backoff.

    An empty result list is an answer - the place is not in Nominatim - and returns None.
    A transport or server error is not an answer: after MAX_ATTEMPTS it raises
    GeocodeUnavailable rather than letting the caller record the place as 'not found'."""
    q = urllib.parse.urlencode({**params, "format": "json", "limit": 1, "addressdetails": 1})
    req = urllib.request.Request(f"https://nominatim.openstreetmap.org/search?{q}",
                                 headers={"User-Agent": "FluTyper-metadata/1.0"})
    last = None
    for attempt in range(MAX_ATTEMPTS):
        try:
            hits = json.load(urllib.request.urlopen(req, timeout=30))
            time.sleep(1.1)  # Nominatim usage policy
            return hits[0] if hits else None
        except urllib.error.HTTPError as exc:
            last = f"HTTP Error {exc.code}: {exc.reason}"
            if exc.code not in RETRY_STATUS:
                raise GeocodeUnavailable(last)
        except (urllib.error.URLError, TimeoutError, ValueError) as exc:
            last = f"{type(exc).__name__}: {exc}"
        if attempt < MAX_ATTEMPTS - 1:
            time.sleep(2 ** (attempt + 1))
    raise GeocodeUnavailable(f"{last}, still failing after {MAX_ATTEMPTS} attempts")


def load_cache(path):
    """Warm the in-memory lookup caches from a previous run's cache file, if there is one."""
    if not path or not os.path.exists(path):
        return
    with open(path) as fh:
        data = json.load(fh)
    for key, value in data.get("places", {}).items():
        name, _, country = key.partition("\t")
        _cache[(name, country)] = tuple(value) if value else None
    _country_codes.update(data.get("countries", {}))


def save_cache(path):
    """Persist every lookup made so far, including the ones that came back 'not found'."""
    if not path:
        return
    data = {"places": {f"{name}\t{country}": (list(v) if v else None)
                       for (name, country), v in _cache.items()},
            "countries": _country_codes}
    with open(path, "w") as fh:
        json.dump(data, fh, indent=0, sort_keys=True)


def country_code(country):
    """ISO code of a country, so place lookups stay inside it ('Denmark' alone resolves to Denmark, SC, USA)."""
    if country not in _country_codes:
        hit = _nominatim({"country": country})
        _country_codes[country] = hit["address"].get("country_code") if hit else None
    return _country_codes[country]


def geocode(name, country):
    key = (name, country)
    if key not in _cache:
        cc = country_code(country)
        if name == country:
            hit = _nominatim({"country": country})
        else:
            hit = _nominatim({"q": name, **({"countrycodes": cc} if cc else {"q": f"{name}, {country}"})})
        _cache[key] = (float(hit["lat"]), float(hit["lon"])) if hit else None
    return _cache[key]


def parse_location(loc):
    """'Europe / Spain / Galicia / Pontevedra / Vigo' -> (place, region, country, geocoded).

    place is the deepest level that geocodes. When no level geocodes, place is the deepest
    level exactly as GISAID wrote it and geocoded is False: the sample keeps the location it
    was reported with (never a guess, never blank) and the caller warns about it instead.
    lower-case levels ('swine') are GISAID noise."""
    levels = [FIXES.get(x, x) for x in (s.strip() for s in str(loc).split(" / "))]
    # lower-case ('swine') and numeric ('0') levels are GISAID noise
    levels = [x.title() if x.isupper() else x for x in levels if x and not x.islower() and not x.isdigit()]
    country = levels[1] if len(levels) > 1 else levels[0]
    region = levels[2] if len(levels) > 2 else country
    for place in reversed(levels[1:]):
        if geocode(place, country):
            return place, region, country, True
    return (levels[-1] if levels else ""), region, country, False


def age_group(age, unit):
    if pd.isna(age):
        return ""
    years = age / 12 if unit == "M" else age
    if years < 3:
        return "0-2"
    if years < 5:
        return "3-4"
    if years < 15:
        return "5-14"
    if years <= 65:
        return "15-65"
    return ">65"


def main(xls, fasta, out_csv, out_tsv, cache_path=None):
    if cache_path is None:
        cache_path = f"{out_csv}.geocache.json"
    load_cache(cache_path)
    try:
        _convert(xls, fasta, out_csv, out_tsv)
    except GeocodeUnavailable as exc:
        save_cache(cache_path)
        sys.exit(f"gisaid_to_metadata: Nominatim is unavailable ({exc}). The lookups already "
                 f"resolved are cached in {cache_path}; rerun the same command to continue "
                 f"from there.")
    save_cache(cache_path)


def _convert(xls, fasta, out_csv, out_tsv):
    labs = {}
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                fields = line[1:].rstrip().split("|")
                labs[isolate_key(fields[0])] = fields[-1].replace("_", " ")

    g = read_export(xls)
    g["Isolate_Id"] = g["Isolate_Id"].map(isolate_key)
    g = g[g["Isolate_Id"].isin(labs)]
    missing = set(labs) - set(g["Isolate_Id"])
    if missing:
        print(f"WARNING: {len(missing)} FASTA IDs not in the GISAID table", file=sys.stderr)

    places = {loc: parse_location(loc) for loc in g["Location"].dropna().unique()}
    geo = g["Location"].map(places)

    # A date column read from a real .xls arrives as a Timestamp; write it back as YYYY-MM-DD.
    # A partial GISAID date ('2025-06') is left exactly as written rather than completed.
    dates = g["Collection_Date"]
    if pd.api.types.is_datetime64_any_dtype(dates):
        dates = dates.dt.strftime("%Y-%m-%d")

    meta = pd.DataFrame({
        "ID": g["Isolate_Id"],
        "DATE": dates,
        "LOCATION": geo.str[0],
        "COUNTRY": geo.str[2],
        **{out: (g[src].fillna("").astype(str) if src in g.columns else "")
           for out, src in PASSTHROUGH.items()},
        "AGE GROUP": [age_group(a, u) for a, u in zip(g["Host_Age"], g["Host_Age_Unit"])],
        "SEX": g["Host_Gender"].map({"Male": "M", "Female": "F"}).fillna(""),
        "ORIGINATING_LAB": g["Isolate_Id"].map(labs),
    })
    meta.to_csv(out_csv, index=False)

    # Coordinates table: one row per place, per region (the "Province" map level), per country (the "Country"
    # map level: Población == Comunidad) and per lab. Población must stay unique, so a country already present as a
    # place/region row (same name) is not added again.
    rows = {}
    for place, region, country, _ok in places.values():
        for name in (place, region):
            if name and name not in rows and geocode(name, country):
                rows[name] = (country, region, name, *geocode(name, country))
    for _, _, country, _ok in places.values():
        if country not in rows and geocode(country, country):
            rows[country] = (country, "", country, *geocode(country, country))
    lab_place = (pd.DataFrame({"lab": meta["ORIGINATING_LAB"], "geo": geo})
                 .dropna().query("lab != ''").groupby("lab")["geo"].agg(lambda s: s.mode()[0]))
    for lab, (place, region, country, ok) in lab_place.items():
        if place and ok:
            if lab not in rows:
                rows[lab] = (country, region, lab, *geocode(place, country))
    coords = pd.DataFrame(rows.values(), columns=["Comunidad", "Provincia", "Población", "Latitud", "Longitud"])
    coords["Altitud"] = ""
    coords.to_csv(out_tsv, sep="\t", index=False)

    unresolved = [loc for loc, (_, _, _, ok) in places.items() if not ok]
    if unresolved:
        print(f"WARNING: could not geocode {unresolved}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    main(*sys.argv[1:])
