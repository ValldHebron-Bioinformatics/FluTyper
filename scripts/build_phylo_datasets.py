#!/usr/bin/env python3
"""Build the H5 context datasets for a reassortment screen.

Three datasets, each a FASTA plus a FluTyper --metadata CSV, each one containing the focus
sample so that its position can be read against both contexts:

  lpai      every isolate classified LPAI
  hpai      every isolate whose GISAID CLADE is 2.3.4.4b
  combined  the union of the two, which is the dataset that can actually show a focus
            segment sitting inside the HPAI clade

The pathotype of an isolate is the HA cleavage site read (scripts/ha_cleavage_site.py) when
that read produced a call, and GISAID's Pathogenicity field only when it did not. A cleavage
site that could not be read is 'indeterminate', which is not a call. An isolate with neither
is UNKNOWN: it joins no dataset and is listed in dataset_report.txt, never assigned by
assumption.

Usage: build_phylo_datasets.py <gisaid_meta.csv> <ha_cleavage.csv> <gisaid.fasta> \
                               <own.fasta> <own_meta.csv> <focus_id> <outdir>
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gisaid_to_metadata import isolate_key  # noqa: E402

HPAI_CLADE = "2.3.4.4b"
DATASETS = ("lpai", "hpai", "combined")
META_COLUMNS = ["ID", "DATE", "LOCATION", "COUNTRY", "HOST", "SUBTYPE", "PATHOGENICITY",
                "CLADE", "GENOTYPE", "AGE GROUP", "SEX", "ORIGINATING_LAB"]


def _key(raw):
    """Join key for an ID: a bare accession ('20060005') becomes 'EPI_ISL_20060005'; an own
    sample ID ('25CR6244') is left alone."""
    text = str(raw).strip()
    return f"EPI_ISL_{text}" if text.isdigit() else text


def classify_with_basis(row, cleavage_by_id):
    """(call, basis) for one metadata row: call is 'LPAI', 'HPAI' or 'UNKNOWN'; basis is
    'cleavage read', 'GISAID field' or 'UNKNOWN'.

    The HA cleavage site read wins: it is read off this very sequence, while GISAID's
    Pathogenicity is whatever the submitter recorded. 'indeterminate' is not a call, so the
    GISAID field is used for those, and an isolate with neither stays UNKNOWN."""
    call = str(cleavage_by_id.get(_key(row.get("ID", "")), "") or "").strip().upper()
    if call in ("LPAI", "HPAI"):
        return call, "cleavage read"
    field = str(row.get("PATHOGENICITY", "") or "").strip().upper()
    if field in ("LPAI", "HPAI"):
        return field, "GISAID field"
    return "UNKNOWN", "UNKNOWN"


def classify(row, cleavage_by_id):
    """'LPAI', 'HPAI' or 'UNKNOWN' for one metadata row (see classify_with_basis)."""
    return classify_with_basis(row, cleavage_by_id)[0]


def _id_column(df, where):
    for c in df.columns:
        if c.lower().replace("_", "").replace(" ", "") in ("id", "sampleid"):
            return c
    sys.exit(f"{where}: no sample ID column ('ID', 'Sample_ID' or 'SampleID'). "
             f"Columns: {list(df.columns)}")


def read_table(path):
    """Read a CSV whose header and values may carry stray spaces, as text throughout."""
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    df.columns = [str(c).strip() for c in df.columns]
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
    return df


# Pathotype column names accepted in the cleavage table, best first
PATHOTYPE_COLUMNS = ("pathotype_cleavage", "pathotype", "pathogenicity")


def read_calls(path):
    """The HA cleavage site pathotype per sample (scripts/ha_cleavage_site.py output)."""
    df = read_table(path)
    id_col = _id_column(df, path)
    lower = {c.lower(): c for c in df.columns}
    patho = next((lower[name] for name in PATHOTYPE_COLUMNS if name in lower), None)
    if patho is None:
        sys.exit(f"{path}: no pathotype column (looked for "
                 f"{', '.join(PATHOTYPE_COLUMNS)}). Columns: {list(df.columns)}")
    return {_key(i): c for i, c in zip(df[id_col], df[patho])}


def gisaid_record_id(header):
    return isolate_key(header.split("|")[0])


def own_record_id(header):
    return header.split("|")[0].split("_")[0]


def read_fasta_by_id(path, key_fn, wanted=None):
    """Group a FASTA's records by sample, keeping each record's header and sequence verbatim."""
    records, current = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                key = key_fn(line[1:].strip())
                current = None if (wanted is not None and key not in wanted) else \
                    records.setdefault(key, [])
                if current is not None:
                    current.append(line.rstrip("\n"))
            elif current is not None:
                current.append(line.rstrip("\n"))
    return records


def build_all(gisaid_meta, cleavage_csv, gisaid_fasta, own_fasta, own_meta, focus_id, outdir):
    """Write the three datasets into outdir and return their counts."""
    meta = read_table(gisaid_meta)
    meta_id = _id_column(meta, gisaid_meta)
    if meta_id != "ID":
        meta = meta.rename(columns={meta_id: "ID"})
    for col in META_COLUMNS:
        if col not in meta.columns:
            meta[col] = ""
    calls = read_calls(cleavage_csv)

    groups = {name: [] for name in ("LPAI", "HPAI", "UNKNOWN")}
    called = {}
    basis_count = {"cleavage read": 0, "GISAID field": 0, "UNKNOWN": 0}
    for row in meta.to_dict("records"):
        call, basis = classify_with_basis(row, calls)
        groups[call].append(row["ID"])
        if row["ID"] not in called:
            called[row["ID"]] = call
            basis_count[basis] += 1
    hpai_ids = [row["ID"] for row in meta.to_dict("records")
                if row.get("CLADE", "").strip() == HPAI_CLADE]

    members = {
        "lpai": list(dict.fromkeys(groups["LPAI"])),
        "hpai": list(dict.fromkeys(hpai_ids)),
    }
    members["combined"] = list(dict.fromkeys(members["lpai"] + members["hpai"]))
    overlap = len(set(members["lpai"]) & set(members["hpai"]))

    # Focus sample: its own records and its own metadata row, with the pathotype FluTyper
    # called for it. It is the sample under investigation, so it is in every dataset.
    focus_records = read_fasta_by_id(own_fasta, own_record_id, {focus_id}).get(focus_id)
    if not focus_records:
        sys.exit(f"{own_fasta}: no records for the focus sample {focus_id}.")
    own = read_table(own_meta)
    own = own.rename(columns={_id_column(own, own_meta): "ID"})
    focus_rows = own[own["ID"] == focus_id]
    if focus_rows.empty:
        sys.exit(f"{own_meta}: no metadata row for the focus sample {focus_id}.")
    focus_meta = {col: "" for col in META_COLUMNS}
    for col in META_COLUMNS:
        if col in focus_rows.columns:
            focus_meta[col] = focus_rows.iloc[0][col]
    focus_meta["ID"] = focus_id
    focus_meta["PATHOGENICITY"] = str(calls.get(_key(focus_id), "")).strip()

    wanted = set(members["combined"])
    context_records = read_fasta_by_id(gisaid_fasta, gisaid_record_id, wanted)
    no_sequence = sorted(wanted - set(context_records))

    os.makedirs(outdir, exist_ok=True)
    counts = {"overlap": overlap, "unknown": sorted(set(groups["UNKNOWN"])),
              "no_sequence": no_sequence,
              "unknown_in_hpai": sorted(set(groups["UNKNOWN"]) & set(members["hpai"]))}
    indexed = meta.drop_duplicates("ID").set_index("ID")
    for name in DATASETS:
        ids = [i for i in members[name] if i in context_records]
        with open(os.path.join(outdir, f"{name}.fasta"), "w") as fh:
            for i in ids:
                fh.write("\n".join(context_records[i]) + "\n")
            fh.write("\n".join(focus_records) + "\n")
        # PATHOGENICITY is the call the membership was built on, not the submitter's field
        rows = [dict(indexed.loc[i], ID=i, PATHOGENICITY=called[i]) for i in ids] + [focus_meta]
        pd.DataFrame(rows)[META_COLUMNS].to_csv(
            os.path.join(outdir, f"{name}_metadata.csv"), index=False)
        counts[name] = {"context": len(ids), "total": len(ids) + 1}

    with open(os.path.join(outdir, "dataset_report.txt"), "w") as fh:
        fh.write(f"Focus sample: {focus_id} "
                 f"(cleavage pathotype: {focus_meta['PATHOGENICITY'] or 'none'}, "
                 f"{len(focus_records) // 2} segments)\n")
        fh.write(f"Classified: LPAI {len(set(groups['LPAI']))}, "
                 f"clade {HPAI_CLADE} {len(set(hpai_ids))}, "
                 f"UNKNOWN {len(set(groups['UNKNOWN']))}\n")
        fh.write(f"Pathotype call based on: cleavage read {basis_count['cleavage read']}, "
                 f"GISAID field {basis_count['GISAID field']}, UNKNOWN {basis_count['UNKNOWN']}\n")
        fh.write(f"LPAI and clade {HPAI_CLADE} both: {overlap} "
                 f"(counted once in combined)\n\n")
        for name in DATASETS:
            fh.write(f"{name}: {counts[name]['context']} context + 1 focus "
                     f"= {counts[name]['total']}\n")
        unknown_in_hpai = sorted(set(counts["unknown"]) & set(members["hpai"]))
        fh.write(f"\nUNKNOWN (no pathotype call: neither an HA cleavage site read nor a "
                 f"GISAID Pathogenicity): {len(counts['unknown'])}\n")
        fh.write(f"  They are not LPAI context, but HPAI membership is by clade, so "
                 f"{len(unknown_in_hpai)} of them carry clade {HPAI_CLADE} and ARE in the "
                 f"HPAI and combined datasets.\n")
        for i in counts["unknown"]:
            mark = f"  (in HPAI context: clade {HPAI_CLADE})" if i in unknown_in_hpai else ""
            fh.write(f"  {i}{mark}\n")
        if no_sequence:
            fh.write(f"\nSelected but absent from the FASTA: {len(no_sequence)}\n")
            for i in no_sequence:
                fh.write(f"  {i}\n")
    return counts


if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    result = build_all(*sys.argv[1:])
    for dataset in DATASETS:
        print(f"{dataset}: {result[dataset]['total']} samples")
