#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloGeography {
    // Geography for the Auspice export (plan D12), resolved ONCE per run and shared by every tree: loading GeoNames
    // takes ~3-5 s and ~0.5 GB, so it is not repeated per tree. Input is the city only (metadata LOCATION, no GISAID
    // location string parsing); country, division (admin-1) and coordinates are derived offline, first match wins:
    //  1. User tables. params.coordinates (coordenades_cat.tsv: Spanish/Catalan towns, no country column -> Spain; its
    //     Comunidad names are mapped to the GeoNames admin-1 name, "Catalunya" -> "Catalonia", so table and GeoNames
    //     cities share one division). The table is only used when Spain is in scope (--auspiceContinent, or the
    //     sample's COUNTRY). --auspiceLatLongs (augur format, level<TAB>name<TAB>lat<TAB>long) overrides every
    //     coordinate (location, division and country points) and can place a city nothing else knows.
    //  2. GeoNames (geonamescache 3.0.2, cities with population >= 500, CC-BY 4.0) restricted to the continent scope
    //     (--auspiceContinent, default EU), or to the country of an optional metadata COUNTRY column, which overrides
    //     the scope. City names are ambiguous worldwide (Valencia ES/VE, Córdoba ES/AR), hence the tiers: primary-name
    //     matches rank before alternate-name matches (Wien/Vienna, Lisboa/Lisbon; but "Porto" is also an alternate
    //     name of Bordeaux, so a pure most-populous rule is wrong), exact before accent-folded, then most populous.
    //     Division names come from the vendored GeoNames admin1CodesASCII.txt (geonamescache ships codes only).
    //  3. Unresolved -> WARNING, country/division/region stay NA and the city gets no map point.
    // Ambiguity is never silent: >1 candidate in the winning tier logs the choice and the alternatives -- a WARNING
    // when a namesake is a close contender (>= 10% of the chosen city's population or >= 50,000 people, a tie, or no
    // population; D16g: Frankfurt (Oder), 57,107 = 8.8% of am Main, still warns),
    // else INFO -- and every resolution, with ALL alternatives, is written to lat_longs_resolved.tsv for auditing. Metadata is never imputed: the metadata LOCATION
    // is kept as the Auspice location label (Wien stays Wien), a missing LOCATION stays NA, and no coordinate is ever
    // attached to a sample that has none. Only when one city name resolves to two places in the same run (e.g.
    // Córdoba ES and, through COUNTRY, Córdoba AR) is the label disambiguated as "<city> (<country>)", because augur
    // matches lat_longs on the (lower-cased) name and the two would otherwise share a map point.
    // Country and division rows of lat_longs.tsv are representative points for the aggregate map levels (a country at
    // its GeoNames capital, a division at its most populous GeoNames city, unless --auspiceLatLongs gives one); each
    // is logged. augur's bundled lat_longs.tsv is not relied on (no location level; "Catalunya" vs GeoNames
    // "Catalonia"), so every point used is written here.
    // Runs in the Nextstrain env (geonamescache is not in the FluTyper env), per D14 through moduleDir.
    conda params.auspiceCondaEnv ?: "${moduleDir}/../../../envs/nextstrain.yaml"   // --auspiceCondaEnv: pre-built prefix (D21)
    errorStrategy 'ignore'   // a failure (bad --auspiceContinent / --auspiceLatLongs / coordinates table) exits 1 with
                             // a "PhyloGeography: ..." message on stderr and emits nothing; the caller reports it

    input:
    // User files are staged into their own folders: a user file called e.g. lat_longs.tsv or geo.tsv would otherwise be
    // a symlink with an output's name, and writing that output would overwrite the user's file through the link.
    path(metadata, stageAs: "in_metadata/*")          // metadata CSV (LOCATION = city, optional COUNTRY = name or ISO code), or []
    path(coordinates, stageAs: "in_coordinates/*")    // params.coordinates (Comunidad Provincia Población Latitud Longitud Altitud), or []
    path(user_lat_longs, stageAs: "in_lat_longs/*")   // --auspiceLatLongs TSV (level<TAB>name<TAB>lat<TAB>long, no header), or []
    path(admin1)        // RESOURCES/geonames/admin1CodesASCII.txt (GeoNames admin-1 code -> name)
    val(continent)      // --auspiceContinent: comma-separated GeoNames continent codes (AF AN AS EU NA OC SA)
    path(env_report, stageAs: "in_env/*")   // AuspiceEnvCheck's nextstrain_env.txt: this task never starts without a
                                            // verified env (plan D21); its lines head geography.log, so every phylo.log

    output:
    path("lat_longs.tsv"),          emit: lat_longs   // augur export --lat-longs: country/division/location rows
    path("geo.tsv"),                emit: geo         // (location, country_input) -> auspice_location, country, division, region
    path("lat_longs_resolved.tsv"), emit: resolved    // audit table, one row per distinct (LOCATION, COUNTRY)
    path("geography.log"),          emit: log         // INFO/WARNING lines for the phylo.log of every tree

    script:
    """
    python3 - <<'PYEOF'
META = "${metadata}".strip()
COORDS = "${coordinates}".strip()
USER_LL = "${user_lat_longs}".strip()
ADMIN1 = "${admin1}"
CONTINENT = "${continent}"
import csv, re, sys, unicodedata
from importlib.metadata import version
import pandas as pd
import geonamescache

P = "PhyloGeography"
MIN_POP = 500   # geonamescache cities500 (plan D12): small Catalan towns such as Aiguafreda only exist there
MAX_ALTERNATIVES = 5
NA = "NA"

def fail(msg):
    sys.stderr.write(f"{P}: {msg}\\n")
    sys.exit(1)

log_lines = [line.rstrip("\\n") for line in open("${env_report}", encoding="utf-8") if line.strip()]
def say(level, msg):
    log_lines.append(f"{P}: {level} {msg}")

def fold(s):
    # case + accent folding: 'Córdoba' == 'cordoba'; used only after an exact (case-folded) match failed
    return "".join(ch for ch in unicodedata.normalize("NFKD", s.casefold()) if not unicodedata.combining(ch))

def fold_town(s):
    # Same normalisation as GeographicReport.nf's normalize_str, so a town matches the Catalan table exactly as it
    # does in the geographic report ("L'Hospitalet de Llobregat" == "Hospitalet de Llobregat")
    s = fold(str(s).strip())
    s = s.replace("l'", " ").replace("d'", " ").replace("l’", " ").replace("d’", " ")
    s = re.sub(r"\\b(?:el|la|els|les)\\b", " ", s)
    s = re.sub(r"[^a-z0-9]", " ", s)
    return " ".join(s.split())

def coord(x):
    # as written in the source (float repr): no rounding, so table coordinates stay identical to the table
    return x if isinstance(x, str) else repr(float(x))

# ---- GeoNames: countries, continents, admin-1 names ----------------------------------------------------------------
gc = geonamescache.GeonamesCache(min_city_population=MIN_POP)
countries = gc.get_countries()
continents = {code: c.get("name", code) for code, c in gc.get_continents().items()}
continent_codes = [c.strip().upper() for c in CONTINENT.split(",") if c.strip()]
unknown = [c for c in continent_codes if c not in continents]
if unknown or not continent_codes:
    fail(f"unknown --auspiceContinent code(s) {','.join(unknown) or '(empty)'}; valid GeoNames continent codes: "
         f"{', '.join(f'{k} ({v})' for k, v in sorted(continents.items()))}")
scope_cc = {cc for cc, c in countries.items() if c["continentcode"] in continent_codes}

admin1 = {}
with open(ADMIN1, encoding="utf-8") as f:
    for line in f:
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) >= 2:
            admin1[parts[0]] = parts[1]

country_by_key = {}
for cc, c in countries.items():
    for key in (cc, c["iso3"], c["name"]):
        country_by_key[fold(key)] = cc

def country_name(cc):
    return countries[cc]["name"] if cc else NA

def region_name(cc):
    return continents.get(countries[cc]["continentcode"], NA) if cc else NA

def division_of(cc, a1):
    return admin1.get(f"{cc}.{a1}", NA)

# ---- Metadata: distinct (LOCATION, COUNTRY) pairs ------------------------------------------------------------------
keys = {}   # (city or None, country_input or None) -> number of samples
if not META:
    say("WARNING", "no --metadata: no sample has a location, the Auspice map will be empty.")
else:
    # keep_default_na=False: COUNTRY "NA" is Namibia's ISO code, not a missing value; blanks are missing (value())
    df = pd.read_csv(META, dtype=str, skipinitialspace=True, keep_default_na=False)
    df.columns = [str(c).strip().upper() for c in df.columns]
    if "LOCATION" not in df.columns:
        say("WARNING", f"the metadata has no LOCATION column (columns: {', '.join(df.columns)}): no sample has a "
                       "location, the Auspice map will be empty.")
    else:
        has_country = "COUNTRY" in df.columns
        def value(row, col):
            v = row.get(col)
            if v is None or pd.isna(v):
                return None
            v = str(v).strip()
            return v or None
        no_location = 0
        for _, row in df.iterrows():
            city = value(row, "LOCATION")
            given = value(row, "COUNTRY") if has_country else None
            if city is None:
                no_location += 1
                if given is None:
                    continue
            keys[(city, given)] = keys.get((city, given), 0) + 1
        n_cities = len({c for c, _ in keys if c is not None})
        say("INFO", f"{len(df)} samples, {n_cities} distinct LOCATION values"
                    f"{', COUNTRY column present (it overrides the lookup)' if has_country else ''}; continent scope "
                    f"{','.join(continent_codes)} ({len(scope_cc)} countries); GeoNames cities with population >= "
                    f"{MIN_POP} via geonamescache {version('geonamescache')}.")
        if no_location:
            say("INFO", f"{no_location} sample(s) have no LOCATION: their location and division stay NA"
                        " (country too, unless the COUNTRY column gives it).")

# ---- User tables ---------------------------------------------------------------------------------------------------
# params.coordinates (coordenades_cat.tsv: Comunidad Provincia Población Latitud Longitud Altitud) has no country
# column: its rows are Spanish (plan D12). Its Comunidad names are mapped to the GeoNames admin-1 names so that
# table cities and GeoNames cities share one division ("Catalunya" -> "Catalonia", admin1CodesASCII.txt ES.56).
ES_COMUNIDAD = {"catalunya": "56", "cataluna": "56", "catalonia": "56", "comunitat valenciana": "60",
                "comunidad valenciana": "60", "illes balears": "07", "islas baleares": "07", "pais vasco": "59",
                "euskadi": "59", "andalucia": "51", "aragon": "52", "castilla y leon": "55", "castilla la mancha": "54",
                "comunidad de madrid": "29", "madrid": "29", "region de murcia": "31", "murcia": "31",
                "canarias": "53", "principado de asturias": "34", "asturias": "34", "navarra": "32",
                "comunidad foral de navarra": "32", "galicia": "58", "cantabria": "39", "la rioja": "27",
                "extremadura": "57", "ceuta": "CE", "melilla": "ML"}
ES_COMUNIDAD = {fold_town(k): v for k, v in ES_COMUNIDAD.items()}   # same normalisation as the lookup ("la" dropped)
table = {}   # fold_town(town) -> [(town, division, admin1 code, lat, long)]
TABLE_CC = "ES"
if COORDS:
    ct = pd.read_csv(COORDS, sep="\\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
    blank = lambda x: pd.isna(x) or not str(x).strip()
    need = ["Comunidad", "Población", "Latitud", "Longitud"]
    missing = [c for c in need if c not in ct.columns]
    if missing:
        fail(f"--coordinates {COORDS} lacks column(s) {', '.join(missing)} (expected {', '.join(need)})")
    unmapped = set()
    for _, r in ct.iterrows():
        if blank(r["Población"]) or blank(r["Latitud"]) or blank(r["Longitud"]):
            continue
        comunidad = "" if blank(r["Comunidad"]) else str(r["Comunidad"]).strip()
        code = ES_COMUNIDAD.get(fold_town(comunidad))
        if code is None:
            unmapped.add(comunidad)
        division = division_of(TABLE_CC, code) if code else (comunidad or NA)
        try:
            float(r["Latitud"]), float(r["Longitud"])
        except ValueError:
            fail(f"--coordinates {COORDS}: non-numeric coordinates for {r['Población']}: {r['Latitud']}, {r['Longitud']}")
        table.setdefault(fold_town(r["Población"]), []).append((str(r["Población"]).strip(), division, code,
                                                                 r["Latitud"].strip(), r["Longitud"].strip()))
    for c in sorted(unmapped):
        say("WARNING", f"coordinates table Comunidad '{c}' has no GeoNames admin-1 name: kept as written, so it may "
                       "not merge with the same division resolved through GeoNames.")
    if TABLE_CC not in scope_cc:
        say("INFO", f"coordinates table skipped (its places are in {country_name(TABLE_CC)}, outside continent scope "
                    f"{','.join(continent_codes)}) unless a sample's COUNTRY column says {country_name(TABLE_CC)}.")

user_ll = {}   # (level, name.lower()) -> (lat, long, name) -- augur matches lat_longs on lower-case (level, name)
if USER_LL:
    with open(USER_LL, encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) != 4 or parts[0].strip().lower() not in ("region", "country", "division", "location"):
                fail(f"--auspiceLatLongs line {i}: expected level<TAB>name<TAB>lat<TAB>long with level "
                     f"region/country/division/location, got {line.strip()!r}")
            try:
                lat, lon = float(parts[2]), float(parts[3])
            except ValueError:
                fail(f"--auspiceLatLongs line {i}: latitude/longitude are not numbers: {line.strip()!r}")
            if not (-90 <= lat <= 90 and -180 <= lon <= 180):
                fail(f"--auspiceLatLongs line {i}: coordinates out of range: {line.strip()!r}")
            user_ll[(parts[0].strip().lower(), parts[1].strip().lower())] = (parts[2].strip(), parts[3].strip(), parts[1].strip())
    say("INFO", f"--auspiceLatLongs: {len(user_ll)} user coordinate(s), they override every other source.")

# ---- GeoNames city index, restricted to the countries that can be asked for ----------------------------------------
given_cc = {}
for _, given in keys:
    if given is not None and given not in given_cc:
        given_cc[given] = country_by_key.get(fold(given))
        if given_cc[given] is None:
            say("WARNING", f"COUNTRY '{given}' is not a GeoNames country name or ISO code: kept as written, no "
                           "region, and its cities are not looked up.")
index_cc = scope_cc | {cc for cc in given_cc.values() if cc}
cities = [c for c in gc.get_cities().values() if c["countrycode"] in index_cc]
idx = {"name": {}, "alternate": {}, "name-folded": {}, "alternate-folded": {}}
for c in cities:
    name = c["name"]
    idx["name"].setdefault(name.casefold(), []).append(c)
    idx["name-folded"].setdefault(fold(name), []).append(c)
    alts = {a for a in c.get("alternatenames", []) if a and a != name}
    for a in {a.casefold() for a in alts} - {name.casefold()}:
        idx["alternate"].setdefault(a, []).append(c)
    for a in {fold(a) for a in alts} - {fold(name)}:
        idx["alternate-folded"].setdefault(a, []).append(c)

def geonames_lookup(city, allowed_cc):
    # D12: primary-name matches rank before alternate-name matches (else "Porto" -> Bordeaux, which lists Porto as
    # an alternate name); exact (case-folded) before accent-folded; within the first non-empty tier, most populous.
    for tier, key in (("name", city.casefold()), ("alternate", city.casefold()),
                      ("name-folded", fold(city)), ("alternate-folded", fold(city))):
        cand = [c for c in idx[tier].get(key, []) if c["countrycode"] in allowed_cc]
        if cand:
            cand.sort(key=lambda c: (-c["population"], c["geonameid"]))
            return tier, cand
    return None, []

# D16g (user ✔ 2026-09-29): a namesake is a close contender if it has >= 10% of the chosen city's population OR is
# itself a sizeable city (>= 50,000 people): Frankfurt (Oder), 57,107 = 8.8% of Frankfurt am Main, must still warn.
CLOSE_PERCENT = 10
CLOSE_ABSOLUTE = 50000

def is_close(pick, other, close_percent=CLOSE_PERCENT, close_absolute=CLOSE_ABSOLUTE):
    # Integer-exact (100 * 652 >= 10 * 6516 is close, 100 * 628 < 10 * 6307 is not). A missing population (key absent
    # / None) counts as close. GeoNames population 0 means "not recorded" but is compared as a number: it is only close
    # to a chosen city that also has 0 (a tie), otherwise e.g. Porto ES (0) would make Porto PT (252,687) a WARNING.
    a, b = pick.get("population"), other.get("population")
    if a is None or b is None:
        return True
    return 100 * b >= close_percent * a or b >= close_absolute

def describe(c):
    return f"{c['name']} ({c['countrycode']}, {division_of(c['countrycode'], c['admin1code'])}, pop {c['population']})"

# ---- Resolve every distinct (LOCATION, COUNTRY) --------------------------------------------------------------------
resolved = {}
for (city, given), n in sorted(keys.items(), key=lambda kv: (kv[0][0] or "", kv[0][1] or "")):
    r = {"city": city, "given": given, "n": n, "cc": None, "country": NA, "division": NA, "a1": None,
         "lat": None, "long": None, "source": "unresolved", "n_candidates": 0, "alternatives": ""}
    gcc = given_cc.get(given) if given is not None else None
    if given is not None:
        r["cc"], r["country"] = gcc, (country_name(gcc) if gcc else given)
    if city is None:
        r["source"] = "COUNTRY column" if given is not None else "unresolved"
        resolved[(city, given)] = r
        continue
    allowed = {gcc} if gcc else (set() if given is not None else scope_cc)
    hits = table.get(fold_town(city), []) if TABLE_CC in allowed else []
    if hits:
        town, division, code, lat, lon = hits[0]
        r.update(cc=TABLE_CC, country=country_name(TABLE_CC), division=division, a1=code, lat=lat, long=lon,
                 source="coordinates", n_candidates=len(hits))
        if len(hits) > 1:
            # the table has no population: always a WARNING (D16g, missing population = close contender)
            r["alternatives"] = "; ".join(f"{h[0]} ({h[1]}, {h[3]}, {h[4]})" for h in hits[1:])
            say("WARNING", f"'{city}' is ambiguous in the coordinates table ({len(hits)} rows): using the first row, "
                           f"{town} ({lat}, {lon}); others: {r['alternatives']}.")
        # M4 (Phase 6): the table only knows Spain, so a hit must not hide in-scope namesakes elsewhere (Breda NL, Pau FR).
        # D16g rule against the Spanish town's GeoNames population (none known -> close): WARNING if close, else INFO.
        _, others = geonames_lookup(city, allowed - {TABLE_CC})
        if others:
            _, spanish = geonames_lookup(city, {TABLE_CC})
            ref = spanish[0] if spanish else {}
            r["alternatives"] = "; ".join(x for x in (r["alternatives"], "; ".join(describe(c) for c in others)) if x)
            r["n_candidates"] += len(others)
            shown = "; ".join(describe(c) for c in others[:MAX_ALTERNATIVES])
            more = f" and {len(others) - MAX_ALTERNATIVES} more" if len(others) > MAX_ALTERNATIVES else ""
            if any(is_close(ref, c) for c in others):
                say("WARNING", f"'{city}' was placed with the coordinates table ({town}, {division}, Spain), but GeoNames "
                               f"has a close namesake elsewhere in scope (>= {CLOSE_PERCENT}% of the Spanish town's population, "
                               f">= {CLOSE_ABSOLUTE:,} people, or no population to compare): {shown}{more}. Add a COUNTRY "
                               "column (or an --auspiceLatLongs location row) if this sample is not from Spain.")
            else:
                say("INFO", f"'{city}' was placed with the coordinates table ({town}, {division}, Spain); smaller GeoNames "
                            f"namesake(s) elsewhere in scope: {shown}{more}.")
    else:
        tier, cand = geonames_lookup(city, allowed)
        if cand:
            pick = cand[0]
            r.update(cc=pick["countrycode"], country=country_name(pick["countrycode"]),
                     division=division_of(pick["countrycode"], pick["admin1code"]), a1=pick["admin1code"],
                     lat=pick["latitude"], long=pick["longitude"], source=f"geonames:{tier}", n_candidates=len(cand))
            if tier != "name":
                say("INFO", f"'{city}' matched a GeoNames {tier.replace('-', ' ')} name of {pick['name']} ({pick['countrycode']}).")
            if len(cand) > 1:
                # D16g: the audit table keeps EVERY alternative (and n_candidates the full count); the log lists up to
                # MAX_ALTERNATIVES. WARNING only when a namesake is a close contender (see is_close: >= 10% of the
                # chosen population or >= 50,000 people, a tie, or a population missing); smaller namesakes (Porto
                # IT/ES next to Porto PT) are INFO, so real ambiguities are not drowned in noise on Europe-wide data.
                others = cand[1:]
                r["alternatives"] = "; ".join(describe(c) for c in others)
                shown = "; ".join(describe(c) for c in others[:MAX_ALTERNATIVES])
                more = f" and {len(others) - MAX_ALTERNATIVES} more" if len(others) > MAX_ALTERNATIVES else ""
                if any(is_close(pick, c) for c in others):
                    say("WARNING", f"'{city}' is ambiguous ({len(cand)} GeoNames {tier} matches in scope, a namesake "
                                   f"has >= {CLOSE_PERCENT}% of the chosen population or >= {CLOSE_ABSOLUTE:,} people): chose the most populous, "
                                   f"{describe(pick)}; alternatives: {shown}{more}. Add a COUNTRY column or an "
                                   "--auspiceLatLongs location row if this is wrong.")
                else:
                    say("INFO", f"'{city}': chose the most populous of {len(cand)} GeoNames {tier} matches in scope, "
                                f"{describe(pick)}; smaller namesake(s) (< {CLOSE_PERCENT}% of its population and < {CLOSE_ABSOLUTE:,} people): "
                                f"{shown}{more}.")
            if r["division"] == NA:
                say("INFO", f"'{city}' ({describe(pick)}) has no GeoNames admin-1 name: division stays NA.")
        elif given is not None:
            say("WARNING", f"'{city}' not found in {r['country']} (COUNTRY column): country kept, division NA, no map "
                           "point at location level.")
        elif ("location", city.lower()) in user_ll:
            say("WARNING", f"'{city}' not found (coordinates table, GeoNames in continent scope "
                           f"{','.join(continent_codes)}): country, division and region stay NA; its map point comes "
                           "from --auspiceLatLongs only. Add a COUNTRY column to give it a country.")
        else:
            say("WARNING", f"'{city}' not found (coordinates table, GeoNames in continent scope "
                           f"{','.join(continent_codes)}): country, division and region stay NA and it gets no map "
                           "point. Add a COUNTRY column or an --auspiceLatLongs location row.")
    resolved[(city, given)] = r

# ---- Location labels: one Auspice location name must mean one place --------------------------------------------------
# augur matches lat_longs on the lower-cased name, so the same city name resolved to two places (e.g. Córdoba ES and
# Córdoba AR through a COUNTRY column) would share one map point: each place is labelled "<city> (<country>)" instead.
places_by_name = {}
for r in resolved.values():
    if r["city"] is not None and r["lat"] is not None:
        places_by_name.setdefault(r["city"].casefold(), set()).add((r["country"], coord(r["lat"]), coord(r["long"])))
for r in resolved.values():
    r["label"] = r["city"]
    if r["city"] is not None and len(places_by_name.get(r["city"].casefold(), ())) > 1:
        r["label"] = f"{r['city']} ({r['country']})"
for name, places in sorted(places_by_name.items()):
    if len(places) > 1:
        labels = sorted({r["label"] for r in resolved.values() if r["city"] and r["city"].casefold() == name and r["lat"] is not None})
        say("WARNING", f"'{next(r['city'] for r in resolved.values() if r['city'] and r['city'].casefold() == name)}' "
                       f"names {len(places)} different places in this run: shown in Auspice as {', '.join(labels)}.")

# User location coordinates override (and may place a city nothing else knows)
for r in resolved.values():
    if r["city"] is None:
        continue
    u = user_ll.get(("location", r["label"].lower()))
    if u:
        r["lat"], r["long"] = u[0], u[1]
        r["source"] = "auspiceLatLongs" if r["source"] == "unresolved" else f"{r['source']}+auspiceLatLongs"
        say("INFO", f"'{r['label']}' placed at {u[0]}, {u[1]} from --auspiceLatLongs.")

# ---- Representative points for the country and division levels ----------------------------------------------------
# Auspice needs one point per country/division shown on the map. Order: --auspiceLatLongs, then GeoNames: a country
# at its capital (countries.json "capital", primary-name match inside the country), a division at its most populous
# GeoNames city. These place aggregate levels on the map; no sample gets these coordinates as its own.
top_city, top_admin1, capital_city = {}, {}, {}
rank = lambda c: (c["population"], -c["geonameid"])
for c in cities:
    cc, a1 = c["countrycode"], (c["countrycode"], c["admin1code"])
    if cc not in top_city or rank(c) > rank(top_city[cc]):
        top_city[cc] = c
    if a1 not in top_admin1 or rank(c) > rank(top_admin1[a1]):
        top_admin1[a1] = c
    if c["name"].casefold() == countries[cc].get("capital", "").casefold():
        if cc not in capital_city or rank(c) > rank(capital_city[cc]):
            capital_city[cc] = c

latlong_rows = {}
def add_row(level, name, lat, lon):
    latlong_rows.setdefault((level, name.lower()), (level, name, lat, lon))

countries_used = {(r["country"], r["cc"]) for r in resolved.values() if r["country"] != NA}
for name, cc in sorted(countries_used, key=lambda x: x[0]):
    u = user_ll.get(("country", name.lower()))
    if u:
        add_row("country", name, u[0], u[1])
        say("INFO", f"country '{name}' drawn at {u[0]}, {u[1]} from --auspiceLatLongs.")
        continue
    if not cc:
        say("WARNING", f"country '{name}' has no map point (not a GeoNames country; add an --auspiceLatLongs country row).")
        continue
    pick = capital_city.get(cc) or top_city.get(cc)
    if pick:
        add_row("country", name, coord(pick["latitude"]), coord(pick["longitude"]))
        say("INFO", f"country '{name}' drawn at {'its capital ' if cc in capital_city else 'its most populous city '}"
                    f"{pick['name']} ({pick['latitude']:.3f}, {pick['longitude']:.3f}).")
    else:
        say("WARNING", f"country '{name}' has no map point (no GeoNames city; add an --auspiceLatLongs country row).")

divisions_used = {(r["division"], r["cc"], r["a1"]) for r in resolved.values() if r["division"] != NA}
done = set()
for name, cc, a1 in sorted(divisions_used, key=lambda x: (x[0], x[1] or "", x[2] or "")):
    if name.lower() in done:
        continue
    done.add(name.lower())
    u = user_ll.get(("division", name.lower()))
    if u:
        add_row("division", name, u[0], u[1])
        say("INFO", f"division '{name}' drawn at {u[0]}, {u[1]} from --auspiceLatLongs.")
        continue
    pick = top_admin1.get((cc, a1)) if cc and a1 else None
    if pick:
        add_row("division", name, coord(pick["latitude"]), coord(pick["longitude"]))
        say("INFO", f"division '{name}' drawn at its most populous GeoNames city {pick['name']} "
                    f"({pick['latitude']:.3f}, {pick['longitude']:.3f}).")
    else:
        say("WARNING", f"division '{name}' has no map point (add an --auspiceLatLongs division row).")
same_division = {}
for name, cc, a1 in divisions_used:
    same_division.setdefault(name.lower(), set()).add((cc, a1))
for name, codes in sorted(same_division.items()):
    if len(codes) > 1:
        say("WARNING", f"division name '{name}' belongs to {len(codes)} admin-1 areas ({', '.join(sorted(f'{c}.{a}' for c, a in codes))}): "
                       "Auspice shows them as one division with one map point.")

for r in sorted(resolved.values(), key=lambda r: r["label"] or ""):
    if r["city"] is not None and r["lat"] is not None:
        add_row("location", r["label"], coord(r["lat"]), coord(r["long"]))

# ---- Outputs -------------------------------------------------------------------------------------------------------
order = {"country": 0, "division": 1, "location": 2}
with open("lat_longs.tsv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    for row in sorted(latlong_rows.values(), key=lambda x: (order[x[0]], x[1].lower())):
        w.writerow(row)

rows = sorted(resolved.values(), key=lambda r: (r["city"] or "", r["given"] or ""))
with open("geo.tsv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["location", "country_input", "auspice_location", "country", "division", "region"])
    for r in rows:
        # country_input: the COUNTRY as written, "" when there was none ("NA" would be Namibia's ISO code)
        w.writerow([r["city"] or NA, r["given"] or "", r["label"] or NA, r["country"], r["division"], region_name(r["cc"])])

with open("lat_longs_resolved.tsv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["city", "country", "division", "lat", "long", "source", "n_candidates",
                "region", "country_input", "auspice_location", "n_samples", "alternatives"])
    for r in rows:
        if r["city"] is None:
            continue
        lat = NA if r["lat"] is None else coord(r["lat"])
        lon = NA if r["long"] is None else coord(r["long"])
        w.writerow([r["city"], r["country"], r["division"], lat, lon, r["source"], r["n_candidates"],
                    region_name(r["cc"]), r["given"] or "", r["label"], r["n"], r["alternatives"] or NA])

n_city = sum(1 for r in rows if r["city"] is not None)
n_point = sum(1 for r in rows if r["city"] is not None and r["lat"] is not None)
by_source = {}
for r in rows:
    if r["city"] is not None:
        by_source[r["source"]] = by_source.get(r["source"], 0) + 1
say("INFO", f"{n_point} of {n_city} distinct (LOCATION, COUNTRY) values have a map point "
            f"({', '.join(f'{k}: {v}' for k, v in sorted(by_source.items()))}).")
with open("geography.log", "w", encoding="utf-8") as f:
    f.write("\\n".join(log_lines) + ("\\n" if log_lines else ""))
PYEOF
    """
}
