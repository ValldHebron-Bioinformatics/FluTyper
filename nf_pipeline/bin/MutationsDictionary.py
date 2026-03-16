#!/usr/bin/env python3
import argparse, re
import pandas as pd

COL = {
    "H3": "reference_site(H3_numbering)",
    "H1": "reference_H1_site(H1_numbering)",
    "H5": "mature_H5_site(no_signal_peptide)",
    "H7": "H7_NUMBERING",
    "H9": "H9_NUMBERING"
}

def main():
    parser = argparse.ArgumentParser(description="Translate HA markers and save an edited CSV file")
    parser.add_argument("--subtype", required=True, choices=COL.keys())
    parser.add_argument("--markers", required=True)
    parser.add_argument("--dictionary", required=True)
    parser.add_argument("--base", default="H5", choices=COL.keys())
    parser.add_argument("--output", default="EDITED_MARKERS.csv")
    args = parser.parse_args()

    # Load dictionary and create a fast mapping dictionary
    sites = pd.read_csv(args.dictionary, dtype=str).dropna(subset=[COL[args.base], COL[args.subtype]])
    mapping = dict(zip(sites[COL[args.base]].str.strip(), sites[COL[args.subtype]].str.strip()))

    # Load the CSV file
    df = pd.read_csv(args.markers, dtype=str)

    # Translate POSITION column using mapping
    df['TRANSLATED_POSITION'] = df['POSITION'].apply(lambda pos: mapping.get(pos.strip(), '-'))

    # Optionally, combine TRANSLATED_POSITION and AA if needed
    df['TRANSLATED_MARKER'] = df.apply(lambda row: f"{row['TRANSLATED_POSITION']}{row['AA']}" if row['TRANSLATED_POSITION'] != '-' else '-', axis=1)

    # Save the edited DataFrame to CSV
    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()