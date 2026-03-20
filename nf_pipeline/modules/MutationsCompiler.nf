#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsCompiler {
    errorStrategy 'ignore'
    debug true
    input:
    path csv_files

    output:
    tuple path("final_mutations_report.xlsx"), path("relevant_mutations.xlsx") , emit: results
    

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
# Additionally, create a separate Excel file for relevant mutations
# Open an ExcelWriter object to handle creating multiple sheets properly
with pd.ExcelWriter("relevant_mutations.xlsx") as writer:
    
    for protein in unique_proteins:
        # Read the specific sheet into a DataFrame directly
        df = pd.read_excel("final_mutations_report.xlsx", sheet_name=str(protein))
        
        # Calculate unique samples correctly from the isolated Series
        n = len(df['SAMPLE_ID'].dropna().unique())
        
        # If a certain position has mutations in more than 10% of the samples, consider it relevant
        threshold = n * 0.1
        relevant_mutations = df['POSITION'].value_counts()[df['POSITION'].value_counts() > threshold].index.tolist()
        relevant_df = df[df['POSITION'].isin(relevant_mutations)]
        
        # Also include all MARKERS=TRUE in the relevant mutations file
        markers_df = df[df['MARKER'] == True]
        
        # Combine both dataframes and drop duplicates to prevent overlapping rows
        combined_df = pd.concat([relevant_df, markers_df]).drop_duplicates()
        
        # Write the combined data to its specific sheet if it is not empty
        if not combined_df.empty:
            combined_df.to_excel(writer, sheet_name=str(protein), index=False)
"""
}