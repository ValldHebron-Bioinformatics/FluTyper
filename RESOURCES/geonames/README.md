# GeoNames admin-1 names (vendored)

`admin1CodesASCII.txt` maps GeoNames first-level administrative division codes to names
(`<CC>.<admin1 code>\t<name>\t<ASCII name>\t<geonameid>`, no header; e.g. `ES.56	Catalonia	Catalonia	3336901`).
PhyloGeography uses it to turn the `admin1code` of a GeoNames city (geonamescache 3.0.2 ships codes only, no names)
into the Auspice `division` (plan decision D12).

- Source: https://download.geonames.org/export/dump/admin1CodesASCII.txt
- Downloaded: 2026-09-29 (server `Last-Modified: Tue, 29 Sep 2026 01:58:15 GMT`, 151,583 bytes, 3,865 rows)
- SHA-256: `1da92a6323a5fec3176f3f743bf4cf4040fd56a876da55e46fbca23c863aa60a`
- Vendored unmodified (the whole file is small, so no filtering).
- Snapshot mismatch: the city table comes from geonamescache 3.0.2's bundled `cities500.json`, an older GeoNames
  snapshot. 93 of its `CC.admin1` codes are absent from this file, nearly all empty or `00` codes ("no admin-1");
  such cities get division `NA` (logged), never a guessed name.

To refresh: download the same URL, replace the file, update the date, size and checksum above.

## Licence and attribution
GeoNames data, © GeoNames (https://www.geonames.org), licensed under
Creative Commons Attribution 4.0 (https://creativecommons.org/licenses/by/4.0/). Credit GeoNames when publishing
maps or tables derived from it. The city coordinates used by PhyloGeography come from the same source through the
`geonamescache` Python package (MIT-licensed code, GeoNames data under CC-BY 4.0).
