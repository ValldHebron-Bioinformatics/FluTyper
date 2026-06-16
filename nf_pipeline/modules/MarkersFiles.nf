#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MarkersFiles {
    // This process generates marker files for each protein based on the protocol (AVIAN or HUMAN).
    // For AVIAN, it queries a SQLite database to extract mutation information and generates CSV files.
    // For HUMAN, it reads existing CSV files from a specified directory and combines them.
    // It also integrates extra markers from a user-provided CSV file, ensuring that the required columns 
    // are present and merging them with the existing data.
    errorStrategy 'ignore'
    debug true

    input:
    path input_data

    output:
    path "*_markers.csv", emit: results

    script:
    // Resolve the absolute path outside the script block for reliability
    def extra_markers_abs = params.extraMarkers ? file(params.extraMarkers).toAbsolutePath().toString() : ""

    """#!/usr/bin/env python3
import sqlite3
import pandas as pd
import os
import glob

mutations_dataframe = pd.DataFrame()
# For AVIAN protocol, connect to the SQLite database and extract mutation information
if "${params.protocol}" == "AVIAN":
    # Connect directly to the physical database file provided by the previous process
    db_connection = sqlite3.connect('${input_data}')

    mutations_query = '''
    SELECT 
        mm.marker_id AS MARKER_ID, -- Retrieve the marker ID to correctly associate mutations combinations with markers
        m.protein_name, -- Retrieve the protein name for grouping
        m.name AS mutation_name, -- Grabs the mutation name (e.g., "M1:N30D") 
        me.effect_name AS EFFECT, -- Pulls the effect description from the markers_effects table
        me.subtype AS FOUND_IN, -- Extract the subtype where the mutation was found
        me.paper_id AS REFERENCE -- Pulls the reference paper ID for the mutation effect
    FROM mutations m
    JOIN markers_mutations mm ON m.name = mm.mutation_name
    JOIN markers_effects me ON mm.marker_id = me.marker_id
    '''
    mutations_dataframe = pd.read_sql_query(mutations_query, db_connection)
    db_connection.close()

    # Extract the position and amino acid change from the mutation name
    mutations_dataframe[['POSITION', 'AA']] = mutations_dataframe['mutation_name'].str.extract(r'(\\d+)([A-Z]+)')
    mutations_dataframe = mutations_dataframe.dropna(subset=['POSITION', 'AA'])

else:
    # HUMAN protocol: read from the provided directory of CSVs
    csv_files = glob.glob('${input_data}/*_markers.csv')
    df_list = []
    for f in csv_files:
        prot_name = os.path.basename(f).replace('_markers.csv', '').upper()
        tmp_df = pd.read_csv(f)
        tmp_df['protein_name'] = prot_name
        df_list.append(tmp_df)
        
    if df_list:
        mutations_dataframe = pd.concat(df_list, ignore_index=True)


# Extra Markers Integration
target_prots = ["HA1", "HA2", "M1", "M2", "NA", "NP", "NS-1", "NS-2", "PA", "PA-X"," "PB1", "PB1-F2", "PB2"]
extra_file_path = "${extra_markers_abs}"

# Check if the path points to a file
if extra_file_path and os.path.isfile(extra_file_path):
    try:
        extra_df = pd.read_csv(extra_file_path)

        # Standardize column names to uppercase
        extra_df.columns = [c.upper().strip() for c in extra_df.columns]
        
        # Required columns for the extra markers CSV
        required_columns = {'MARKER_ID', 'POSITION', 'AA', 'PROTEIN', 'EFFECT', 'FOUND_IN', 'REFERENCE'}
        actual_columns = set(extra_df.columns)

        # Validate that all required columns are present in the extra markers CSV
        if required_columns.issubset(actual_columns):
            extra_df = extra_df.rename(columns={'PROTEIN': 'protein_name'})
            if not mutations_dataframe.empty:
                mutations_dataframe = pd.concat([mutations_dataframe, extra_df], ignore_index=True)
            else:
                mutations_dataframe = extra_df
        else:
            missing_cols = required_columns - actual_columns
            print(f"WARNING: Skipping extra markers file. Missing required columns: {', '.join(missing_cols)}")
    
    except Exception as e:
        print(f"WARNING: Could not process extra markers file. Error: {e}")

if not mutations_dataframe.empty:
    mutations_dataframe['POSITION'] = pd.to_numeric(mutations_dataframe['POSITION'], errors='coerce')
    mutations_dataframe = mutations_dataframe.dropna(subset=['POSITION'])

    # Remove duplicates based on the combination of marker ID, protein name, position, amino acid change, effect, found in, and reference
    mutations_dataframe = mutations_dataframe.drop_duplicates(subset=['MARKER_ID', 'protein_name', 'POSITION', 'AA', 'EFFECT', 'FOUND_IN', 'REFERENCE'])
    mutations_dataframe['POSITION'] = mutations_dataframe['POSITION'].astype(int)
    mutations_dataframe = mutations_dataframe.sort_values(['MARKER_ID', 'POSITION'])

    
    name_map = {
        "NS-1": "NS1",
        "NS-2": "NS2"
    }

    for protein_id, protein_specific_dataframe in mutations_dataframe.groupby('protein_name'):
        if protein_id in target_prots:
            # Use the map to get the new name, or default to the original protein_id
            export_name = name_map.get(protein_id, protein_id)
            
            protein_specific_dataframe[['MARKER_ID', 'POSITION', 'AA', 'EFFECT', 'FOUND_IN', 'REFERENCE']].to_csv(f"{export_name}_markers.csv", index=False)
    """
}
