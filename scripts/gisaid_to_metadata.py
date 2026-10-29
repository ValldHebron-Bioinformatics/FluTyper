#!/usr/bin/env python3
"""Convert a GISAID EpiFlu isolates export (.xls) into FluTyper inputs:

  out.csv  --metadata     ID,DATE,LOCATION,AGE GROUP,SEX,ORIGINATING_LAB
  out.tsv  --coordinates  same columns as RESOURCES/coordenades_cat.tsv, covering
                          every location, region and lab of the metadata

ORIGINATING_LAB is not in the GISAID table, so it is taken from the last field
of the FASTA headers (>ID|SEGMENT|DATE|NAME|LAB) of the same download.
Places are geocoded with OpenStreetMap Nominatim (place names only, 1 query/s).

Usage: gisaid_to_metadata.py gisaid_epiflu_isolates.xls sequences.fasta out.csv out.tsv
"""
import json
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

# GISAID spellings that do not geocode or would split one region in two
FIXES = {"Basqsue Country": "Basque Country", "Catalunya": "Catalonia",
         # Greek prefectures come in the genitive with a "Nomos" prefix
         "Nomos Attica": "Attica", "Nomos Ioanninon": "Ioannina", "Nomos Lasithiou": "Lasithi"}

_cache = {}
_country_codes = {}


def _nominatim(params):
    q = urllib.parse.urlencode({**params, "format": "json", "limit": 1, "addressdetails": 1})
    req = urllib.request.Request(f"https://nominatim.openstreetmap.org/search?{q}",
                                 headers={"User-Agent": "FluTyper-metadata/1.0"})
    hits = json.load(urllib.request.urlopen(req, timeout=30))
    time.sleep(1.1)  # Nominatim usage policy
    return hits[0] if hits else None


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
    """'Europe / Spain / Galicia / Pontevedra / Vigo' -> (place, region, country).

    place is the deepest level that geocodes; lower-case levels ('swine') are GISAID noise."""
    levels = [FIXES.get(x, x) for x in (s.strip() for s in str(loc).split(" / "))]
    # lower-case ('swine') and numeric ('0') levels are GISAID noise
    levels = [x.title() if x.isupper() else x for x in levels if x and not x.islower() and not x.isdigit()]
    country = levels[1] if len(levels) > 1 else levels[0]
    region = levels[2] if len(levels) > 2 else country
    for place in reversed(levels[1:]):
        if geocode(place, country):
            return place, region, country
    return "", region, country


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


def main(xls, fasta, out_csv, out_tsv):
    labs = {}
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                fields = line[1:].rstrip().split("|")
                labs[fields[0]] = fields[-1].replace("_", " ")

    g = pd.read_excel(xls, dtype={"Isolate_Id": str})
    g["Isolate_Id"] = g["Isolate_Id"].str.replace("EPI_ISL_", "", regex=False)
    g = g[g["Isolate_Id"].isin(labs)]
    missing = set(labs) - set(g["Isolate_Id"])
    if missing:
        print(f"WARNING: {len(missing)} FASTA IDs not in the GISAID table", file=sys.stderr)

    places = {loc: parse_location(loc) for loc in g["Location"].dropna().unique()}
    geo = g["Location"].map(places)

    meta = pd.DataFrame({
        "ID": g["Isolate_Id"],
        "DATE": g["Collection_Date"],
        "LOCATION": geo.str[0],
        "AGE GROUP": [age_group(a, u) for a, u in zip(g["Host_Age"], g["Host_Age_Unit"])],
        "SEX": g["Host_Gender"].map({"Male": "M", "Female": "F"}).fillna(""),
        "ORIGINATING_LAB": g["Isolate_Id"].map(labs),
    })
    meta.to_csv(out_csv, index=False)

    # Coordinates table: one row per place, per region (the "Province" map level), per country (the "Country"
    # map level: Población == Comunidad) and per lab. Población must stay unique, so a country already present as a
    # place/region row (same name) is not added again.
    rows = {}
    for place, region, country in places.values():
        for name in (place, region):
            if name and name not in rows and geocode(name, country):
                rows[name] = (country, region, name, *geocode(name, country))
    for _, _, country in places.values():
        if country not in rows and geocode(country, country):
            rows[country] = (country, "", country, *geocode(country, country))
    lab_place = (pd.DataFrame({"lab": meta["ORIGINATING_LAB"], "geo": geo})
                 .dropna().query("lab != ''").groupby("lab")["geo"].agg(lambda s: s.mode()[0]))
    for lab, (place, region, country) in lab_place.items():
        if place:
            if lab not in rows:
                rows[lab] = (country, region, lab, *geocode(place, country))
    coords = pd.DataFrame(rows.values(), columns=["Comunidad", "Provincia", "Población", "Latitud", "Longitud"])
    coords["Altitud"] = ""
    coords.to_csv(out_tsv, sep="\t", index=False)

    unresolved = [loc for loc, (place, _, _) in places.items() if not place]
    if unresolved:
        print(f"WARNING: could not geocode {unresolved}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
