process MarkersFiles {
    errorStrategy 'ignore'
    debug true

    input:
    path 'flumut_db.sqlite'

    output:
    path "*_markers.csv"

    script:
    """#!/usr/bin/env python3
import sqlite3
import pandas as pd

# Connect directly to the physical database file provided by the previous process
db_connection = sqlite3.connect('flumut_db.sqlite')

mutations_query = '''
SELECT 
    m.protein_name, -- Retrieve the protein name for grouping
    m.name AS mutation_name, -- Grabs the mutation name (e.g., "M1:N30D") 
    me.effect_name AS EFFECT, -- Pulls the effect description from the markers_effects table
    me.paper_id AS REFERENCE -- Pulls the reference paper ID for the mutation effect
FROM mutations m -- Set the mutations table as the primary source of mutation data
JOIN markers_mutations mm ON m.name = mm.mutation_name -- Join to link mutations to their associated markers
JOIN markers_effects me ON mm.marker_id = me.marker_id -- Join to get the effect details for each mutation
'''
mutations_dataframe = pd.read_sql_query(mutations_query, db_connection)
db_connection.close()

# Extract the position and amino acid change from the mutation name, then clean and sort the dataframe
mutations_dataframe[['POSITION', 'AA']] = mutations_dataframe['mutation_name'].str.extract(r'(\\d+)([A-Z]+)')
mutations_dataframe = mutations_dataframe.dropna(subset=['POSITION', 'AA'])
mutations_dataframe['POSITION'] = mutations_dataframe['POSITION'].astype(int)
mutations_dataframe = mutations_dataframe.sort_values('POSITION')

for protein_id, protein_specific_dataframe in mutations_dataframe.groupby('protein_name'):
    protein_specific_dataframe[['POSITION', 'AA', 'EFFECT', 'REFERENCE']].to_csv(f"{protein_id}_markers.csv", index=False)
    """
}