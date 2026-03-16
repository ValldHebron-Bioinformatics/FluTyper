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
    parser = argparse.ArgumentParser(description="Translate HA markers and save an edited Excel file")
    parser.add_argument("--subtype", required=True, choices=COL.keys())
    parser.add_argument("--markers", required=True)
    parser.add_argument("--dictionary", required=True)
    parser.add_argument("--base", default="H5", choices=COL.keys())
    parser.add_argument("--output", default="EDITED_MARKERS.xlsx")
    args = parser.parse_args()

    # Load dictionary and create a fast mapping dictionary
    sites = pd.read_csv(args.dictionary, dtype=str).dropna(subset=[COL[args.base], COL[args.subtype]])
    mapping = dict(zip(sites[COL[args.base]].str.strip(), sites[COL[args.subtype]].str.strip()))

    # Load the Excel file
    df = pd.read_excel(args.markers, sheet_name="Markers", header=None)
    rx = re.compile(r"^\s*(-?\d+[A-Za-z]?)\s*([A-Za-z])\s*$")

    # Translation function applied to each cell
    def translate_cell(val):
        match = rx.match(str(val))
        if match:
            pos, aa = match.groups()
            target_pos = mapping.get(pos, "")
            return f"{target_pos}{aa.upper()}" if target_pos and target_pos != "-" else "-"
        return val

    # Find HA columns and modify them in place starting from row 2
    for col in df.columns:
        if str(df.iat[1, col]).strip().upper() == "HA":
            df.iloc[2:, col] = df.iloc[2:, col].apply(translate_cell)

    # Save the entirely edited sheet back to Excel
    df.to_excel(args.output, index=False, header=False)
    print(f"Modifications saved successfully to {args.output}")

if __name__ == "__main__":
    main()