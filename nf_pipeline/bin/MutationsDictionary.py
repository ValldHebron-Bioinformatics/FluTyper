#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--subtype", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--dictionary", required=True)
    parser.add_argument("--base", default="H5")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Prevents turning "NA" protein into a NaN missing value
    df = pd.read_csv(args.input, dtype=str, keep_default_na=False)
    
    if df.empty:
        df['POSITION_SUBTYPE'] = None
        df.to_csv(args.output, index=False)
        return

    # Read the new mutations dictionary
    m_dict = pd.read_csv(args.dictionary, dtype=str, keep_default_na=False)
    m_dict.columns = m_dict.columns.str.strip()
    
    # Clean the PROTEIN column
    if 'PROTEIN' in m_dict.columns:
        m_dict['PROTEIN'] = m_dict['PROTEIN'].str.strip().str.upper()

    base_prefix = f"{args.base.upper()}_"
    subtype_prefix = f"{args.subtype.upper()}_"

    start_col = next((col for col in m_dict.columns if col.startswith(base_prefix)), None)
    target_col = next((col for col in m_dict.columns if col.startswith(subtype_prefix)), None)

    if not start_col:
        print(f"ERROR: No column starting with '{base_prefix}' found in the mutations dictionary.")
        sys.exit(1)

    if not target_col:
        df['POSITION_SUBTYPE'] = df['POSITION']
        df.to_csv(args.output, index=False)
        return

    lookup = m_dict.set_index(['PROTEIN', start_col])[target_col].to_dict()

    prots = df['PROTEIN'].str.strip().str.upper()
    
    df['POSITION_SUBTYPE'] = [lookup.get((p, str(pos).strip()), str(pos).strip()) for p, pos in zip(prots, df['POSITION'])]

    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()