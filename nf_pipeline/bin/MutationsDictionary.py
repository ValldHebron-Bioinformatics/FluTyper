#!/usr/bin/env python3
import argparse
import os
import pandas as pd

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

    # Determine protein context from filename if column is missing
    filename = os.path.basename(args.input).upper()
    file_prot = "HA2" if filename.startswith("HA2") else "HA1"

    # Load dictionary and create the lookup map
    master_dict = pd.read_csv(args.dictionary, dtype=str)
    start_col = SUBTYPE_COLUMNS.get(args.base)
    target_col = SUBTYPE_COLUMNS.get(args.subtype)
    
    master_dict['region'] = master_dict['region'].fillna("HA1").str.strip().str.upper()
    lookup = master_dict.set_index(['region', start_col])[target_col].to_dict()

    # Load input and preserve original data
    df = pd.read_csv(args.input, dtype=str)
    prots = df['PROTEIN'].str.strip().str.upper() if 'PROTEIN' in df.columns else [file_prot] * len(df)
    
    # Create the new column while keeping the original POSITION
    new_col_name = f"POSITION_{args.subtype}"
    df[new_col_name] = [lookup.get((p, pos), "-") for p, pos in zip(prots, df['POSITION'])]
    
    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()