process MarkersFiles {
    errorStrategy 'ignore'

    input:
    path 'flumut_db.sqlite'

    output:
    path "*_markers.csv", emit: results

    script:
    // Resolve the absolute path outside the script block for reliability
    def extra_markers_abs = params.extraMarkers ? file(params.extraMarkers).toAbsolutePath().toString() : ""

    """#!/usr/bin/env python3
import sqlite3
import pandas as pd
import os

# Connect directly to the physical database file provided by the previous process
db_connection = sqlite3.connect('flumut_db.sqlite')

mutations_query = '''
SELECT 
    m.protein_name, -- Retrieve the protein name for grouping
    m.name AS mutation_name, -- Grabs the mutation name (e.g., "M1:N30D") 
    me.effect_name AS EFFECT, -- Pulls the effect description from the markers_effects table
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

# Extra Markers Integration
target_prots = ["HA1", "HA2", "M1", "M2", "NA", "NP", "NS-1", "NS-2", "PA", "PB1", "PB1-F2", "PB2"]
extra_file_path = "${extra_markers_abs}"

# Check if the path points to a file
if extra_file_path and os.path.isfile(extra_file_path):
    extra_df = pd.read_csv(extra_file_path)
    
    # Standardize column names to uppercase
    extra_df.columns = [c.upper() for c in extra_df.columns]
    
    # Align the PROTEIN column name with the database query output
    if 'PROTEIN' in extra_df.columns:
        extra_df = extra_df.rename(columns={'PROTEIN': 'protein_name'})
        
    # Combine the extra markers with the database mutations
    mutations_dataframe = pd.concat([mutations_dataframe, extra_df], ignore_index=True)

# Final cleanup: ensure POSITION is numeric, deduplicate, and sort
mutations_dataframe['POSITION'] = pd.to_numeric(mutations_dataframe['POSITION'], errors='coerce')
mutations_dataframe = mutations_dataframe.dropna(subset=['POSITION'])
mutations_dataframe = mutations_dataframe.drop_duplicates(subset=['protein_name', 'POSITION', 'AA', 'EFFECT'])
mutations_dataframe['POSITION'] = mutations_dataframe['POSITION'].astype(int)
mutations_dataframe = mutations_dataframe.sort_values('POSITION')

# Export grouped files
for protein_id, protein_specific_dataframe in mutations_dataframe.groupby('protein_name'):
    if protein_id in target_prots:
        protein_specific_dataframe[['POSITION', 'AA', 'EFFECT', 'REFERENCE']].to_csv(f"{protein_id}_markers.csv", index=False)
    """
}