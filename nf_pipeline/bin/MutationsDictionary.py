#!/usr/bin/env python3
import argparse
import pandas as pd

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

    # Read the new mutations dictionary
    m_dict = pd.read_csv(args.dictionary, dtype=str)
    m_dict.columns = m_dict.columns.str.strip()
    
    # Clean the PROTEIN column
    m_dict['PROTEIN'] = m_dict['PROTEIN'].str.strip().str.upper()

    # Generate column names dynamically based on the new format
    start_col = f"{args.base.upper()}_numbering"
    target_col = f"{args.subtype.upper()}_numbering"

    # Check if the target subtype exists in the dictionary
    if target_col not in m_dict.columns:
        print(f"Warning: {target_col} not found in the dictionary. Original position will be kept.")
        df['POSITION_SUBTYPE'] = df['POSITION']
        df.to_csv(args.output, index=False)
        return

    # Create the lookup dictionary: (PROTEIN, base_pos) -> target_pos
    lookup = m_dict.set_index(['PROTEIN', start_col])[target_col].to_dict()

    # Determine the protein for each row based on the 'PROTEIN' column in the input
    prots = df['PROTEIN'].str.strip().str.upper()
    
    # Apply conversion with fallback to the original position if unmapped
    df['POSITION_SUBTYPE'] = [lookup.get((p, str(pos).strip()), str(pos).strip()) for p, pos in zip(prots, df['POSITION'])]

    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()