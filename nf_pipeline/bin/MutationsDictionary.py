#!/usr/bin/env python3
import argparse
import pandas as pd

# These are the exact column names expected in the master dictionary file
SUBTYPE_COLUMNS = {
    "H3": "reference_site(H3_numbering)",
    "H1": "reference_H1_site(H1_numbering)",
    "H5": "mature_H5_site(no_signal_peptide)",
    "H7": "H7_NUMBERING",
    "H9": "H9_NUMBERING"
}

def main():
    # Set up the instructions for anyone running the script from the command line
    parser = argparse.ArgumentParser(description="Translate HA marker positions and save them to a new CSV.")
    parser.add_argument("--subtype", required=True, choices=SUBTYPE_COLUMNS.keys(), help="The target subtype you want to translate to.")
    parser.add_argument("--input", required=True, help="Your input CSV file containing a 'POSITION' column.")
    parser.add_argument("--dictionary", required=True, help="The master CSV file that maps all the subtypes together.")
    parser.add_argument("--base", default="H5", choices=SUBTYPE_COLUMNS.keys(), help="The starting subtype of your markers (defaults to H5).")
    parser.add_argument("--output", default="EDITED_MARKERS.csv", help="What to name the final saved file.")
    
    args = parser.parse_args()

    # Find the exact column names for the starting and target subtypes
    starting_column = SUBTYPE_COLUMNS[args.base]
    target_column = SUBTYPE_COLUMNS[args.subtype]

    # Load the master dictionary, ignoring any rows that are missing the required numbers
    master_dictionary = pd.read_csv(args.dictionary, dtype=str) # dtype=str ensures we read everything as text, so no NaNs from empty cells
    valid_rows = master_dictionary.dropna(subset=[starting_column, target_column])

    # Create a dictionary for the translation lookup
    translation_lookup = {}
    for index, row in valid_rows.iterrows():
        start_position = row[starting_column].strip()
        target_position = row[target_column].strip()
        translation_lookup[start_position] = target_position

    # Load the user's marker data
    input_data = pd.read_csv(args.input, dtype=str)

    # Define a quick helper function to translate a single position
    def translate_position(current_position):
        clean_position = current_position.strip()
        return translation_lookup.get(clean_position, "-")

    # Apply the translation to the entire POSITION column
    input_data['POSITION'] = input_data['POSITION'].apply(translate_position)

    # Save the newly translated data
    input_data.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()