#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsCompiler {
    errorStrategy 'ignore'
    debug true
    input:
    path csv_files

    output:
    path "final_mutations_report.xlsx", emit: results

   script:
    """#!/usr/bin/env python3
import pandas as pd
import sys

csv_list = "${csv_files}".split()
dataframes = [pd.read_csv(f) for f in csv_list if f.endswith('.csv')]

if not dataframes:
    print("ERROR: No CSV files found to compile.")
    sys.exit(1)
# Concatenate all dataframes into one master dataframe ignoring the index to avoid duplicate indices from individual files 
master_df = pd.concat(dataframes, ignore_index=True)

if 'PROTEIN' not in master_df.columns:
    print("ERROR: 'PROTEIN' column not found in the data. Cannot compile report.")
    sys.exit(1)

# Write the master dataframe to an Excel file with separate sheets for each unique protein
with pd.ExcelWriter("final_mutations_report.xlsx", engine='openpyxl') as writer:
    unique_proteins = master_df['PROTEIN'].dropna().unique()
    
    for protein in unique_proteins:
        protein_df = master_df[master_df['PROTEIN'] == protein]
        
        # Sort the data for better readability
        if 'SAMPLE_ID' in protein_df.columns and 'POSITION' in protein_df.columns:
            protein_df = protein_df.sort_values(by=['SAMPLE_ID', 'POSITION'])
            
        protein_df.to_excel(writer, sheet_name=str(protein), index=False)
"""
}