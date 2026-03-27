#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import sys

SUBTYPE_COLUMNS = {
    "H3": "reference_site(H3_numbering)",
    "H1": "reference_H1_site(H1_numbering)",
    "H5": "mature_H5_site(no_signal_peptide)",
    "H7": "H7_NUMBERING",
    "H9": "H9_NUMBERING"
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--subtype", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--dictionary", required=True)
    parser.add_argument("--base", default="H5")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, dtype=str)
    
    if df.empty:
        df['POSITION_SUBTYPE'] = None
        df.to_csv(args.output, index=False)
        return

    m_dict = pd.read_csv(args.dictionary, dtype=str)
    m_dict.columns = m_dict.columns.str.strip()
    m_dict['region'] = m_dict['region'].fillna("HA1").str.strip().str.upper()
    # Assume any region containing "HA2" is HA2, otherwise HA1
    m_dict['region'] = m_dict['region'].apply(lambda x: "HA2" if "HA2" in x else "HA1")

    start_col = SUBTYPE_COLUMNS.get(args.base.upper(), SUBTYPE_COLUMNS["H5"])
    target_col = SUBTYPE_COLUMNS.get(args.subtype.upper(), SUBTYPE_COLUMNS["H5"])
    # Create a lookup dictionary for (region, start_col) -> target_col
    lookup = m_dict.set_index(['region', start_col])[target_col].to_dict()

    fname = os.path.basename(args.input).upper()
    file_prot = "HA2" if "HA2" in fname else "HA1"

    # Determine the protein for each row based on the 'PROTEIN' column
    prots = df['PROTEIN'].str.strip().str.upper()
    
    df['POSITION_SUBTYPE'] = [lookup.get((p, str(pos).strip()), str(pos).strip()) for p, pos in zip(prots, df['POSITION'])]

    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()