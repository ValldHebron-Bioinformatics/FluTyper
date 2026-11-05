#!/usr/bin/env python3
"""Put the three pathotype sources side by side, isolate by isolate, without merging them.

  cleavage_motif        the motif read at the HA1/HA2 junction, e.g. PLREKRRKR/GLF
  pathotype_cleavage    LPAI / HPAI / indeterminate, read off that motif
                        (scripts/ha_cleavage_site.py). This is the sequence-derived call.
  pathotype_nextclade   the HPAI/LPAI suffix of the NEAREST REFERENCE found by
                        `nextclade sort` (SubtypeDetection). A nearest-neighbour label, not a
                        reading of the sequence, and only set for H5 and H7. Carried for
                        comparison; it never drives the agreement.
  pathotype_gisaid      whatever the submitting laboratory recorded in the GISAID export.

agreement compares pathotype_cleavage with pathotype_gisaid and is one of:
agree, disagree, only_cleavage, only_gisaid, neither. 'indeterminate' is not a call.

Usage: compare_pathotype.py <ha_cleavage.csv> <inferred_subtypes.csv> <gisaid_meta.csv> <out.csv>
"""
import sys

import pandas as pd

AGREEMENT_COLUMNS = ["ID", "inferred_subtype", "cleavage_motif", "pathotype_cleavage",
                     "pathotype_nextclade", "pathotype_gisaid", "clade", "genotype",
                     "agreement"]
CALLS = ("LPAI", "HPAI")


def _id_column(df, where):
    for c in df.columns:
        if str(c).lower().replace("_", "").replace(" ", "") in ("id", "sampleid"):
            return c
    sys.exit(f"{where}: no sample ID column ('ID', 'Sample_ID' or 'SampleID'). "
             f"Columns: {list(df.columns)}")


def _agreement(cleavage, gisaid):
    # 'indeterminate' and '' are both "no call": neither can agree or disagree with anything
    cleavage = cleavage if cleavage in CALLS else ""
    gisaid = gisaid if gisaid in CALLS else ""
    if cleavage and gisaid:
        return "agree" if cleavage == gisaid else "disagree"
    if cleavage:
        return "only_cleavage"
    if gisaid:
        return "only_gisaid"
    return "neither"


def _text(df, column):
    if column in df.columns:
        return df[column].fillna("").astype(str).str.strip()
    return pd.Series([""] * len(df), index=df.index)


def compare(cleavage, nextclade, gisaid):
    """One row per isolate of the union of the three tables, with the three pathotypes."""
    left = pd.DataFrame({
        "ID": _text(cleavage, _id_column(cleavage, "the cleavage table")),
        "cleavage_motif": _text(cleavage, "cleavage_motif"),
        "pathotype_cleavage": _text(cleavage, "pathotype_cleavage"),
    })
    middle = pd.DataFrame({
        "ID": _text(nextclade, _id_column(nextclade, "the inferred subtypes table")),
        "inferred_subtype": _text(nextclade, "inferred_subtype"),
        "pathotype_nextclade": _text(nextclade, "pathotype").str.upper(),
    })
    right = pd.DataFrame({
        "ID": _text(gisaid, _id_column(gisaid, "the GISAID metadata")),
        "pathotype_gisaid": _text(gisaid, "PATHOGENICITY").str.upper(),
        "clade": _text(gisaid, "CLADE"),
        "genotype": _text(gisaid, "GENOTYPE"),
    })
    merged = (left.drop_duplicates("ID")
                  .merge(middle.drop_duplicates("ID"), on="ID", how="outer")
                  .merge(right.drop_duplicates("ID"), on="ID", how="outer")
                  .fillna(""))
    merged["agreement"] = [_agreement(c, g) for c, g in
                           zip(merged["pathotype_cleavage"].str.upper().replace(
                               {"INDETERMINATE": "indeterminate"}),
                               merged["pathotype_gisaid"])]
    return merged[AGREEMENT_COLUMNS].sort_values("ID").reset_index(drop=True)


def main(cleavage_csv, nextclade_csv, gisaid_csv, out_csv):
    tables = []
    for path in (cleavage_csv, nextclade_csv, gisaid_csv):
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
        df.columns = [str(c).strip() for c in df.columns]
        tables.append(df)
    out = compare(*tables)
    out.to_csv(out_csv, index=False)

    counts = out["agreement"].value_counts()
    print(f"{len(out)} isolates compared -> {out_csv}")
    for name in ("agree", "disagree", "only_cleavage", "only_gisaid", "neither"):
        print(f"  {name}: {int(counts.get(name, 0))}")
    disagree = out[out["agreement"] == "disagree"]
    if not disagree.empty:
        print("\ndisagreements (HA cleavage site read vs GISAID Pathogenicity):")
        for _, r in disagree.iterrows():
            print(f"  {r['ID']}  {r['inferred_subtype']}  motif {r['cleavage_motif'] or '-'}  "
                  f"cleavage {r['pathotype_cleavage']} / GISAID {r['pathotype_gisaid']}  "
                  f"clade {r['clade'] or '-'}")
    return out


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
