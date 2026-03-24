#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsCompiler {
    errorStrategy 'ignore'
    debug true
    input:
    path csv_files

    output:
    tuple path("final_mutations_report.xlsx"), path("relevant_mutations.xlsx"), emit: results
    
    script:
    """#!/usr/bin/env python3
import pandas as pd
import sys

csv_list = "${csv_files}".split()
dataframes = [pd.read_csv(f) for f in csv_list if f.endswith('.csv')]

if not dataframes:
    print("ERROR: No CSV files found to compile.")
    sys.exit(1)

# Concatenate all dataframes into one master dataframe
master_df = pd.concat(dataframes, ignore_index=True)

if 'PROTEIN' not in master_df.columns:
    print("ERROR: 'PROTEIN' column not found in the data. Cannot compile report.")
    sys.exit(1)

# Open both Excel files simultaneously to avoid reading data back from the disk
with pd.ExcelWriter("final_mutations_report.xlsx", engine='openpyxl') as writer_all, \\
     pd.ExcelWriter("relevant_mutations.xlsx", engine='openpyxl') as writer_rel:
    
    # Group the dataframe by the 'PROTEIN' column
    # When you call .groupby('PROTEIN'), you aren't just looking at a list,
    # but rather you are creating a DataFrameGroupBy object. 
    # This object acts like a dictionary where the "keys" are the protein names and the "values" 
    # are the actual rows of data belonging to those names.
    for protein, protein_df in master_df.dropna(subset=['PROTEIN']).groupby('PROTEIN'):
        
        # Write to the final mutations report
        if 'SAMPLE_ID' in protein_df.columns and 'POSITION' in protein_df.columns:
            protein_df = protein_df.sort_values(by=['SAMPLE_ID', 'POSITION'])
            
        protein_df.to_excel(writer_all, sheet_name=str(protein), index=False)
        
        # Filter and write to the relevant mutations report, nunique() returns the number of unique sample IDs
        threshold = protein_df['SAMPLE_ID'].nunique() * 0.1 # 10% of the samples have the mutations, subject to change ASK ALEJANDRA
        
        counts = protein_df['POSITION'].value_counts()
        relevant_positions = counts[counts > threshold].index
        # A mutation is considered relevant if it occurs in more than 10% of the samples for that protein, or if it is marked as a marker mutation
        relevants = protein_df['POSITION'].isin(relevant_positions) | (protein_df['MARKER'] == True)
        combined_df = protein_df[relevants]
        
        if not combined_df.empty:
            combined_df.to_excel(writer_rel, sheet_name=str(protein), index=False)
"""
}