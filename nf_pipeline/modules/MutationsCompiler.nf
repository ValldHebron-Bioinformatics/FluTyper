#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsCompiler {
    errorStrategy 'ignore'

    input:
    path csv_files

    output:
    tuple path("final_mutations_report.xlsx"), path("relevant_mutations.xlsx"), emit: results

    
    script:
    """#!/usr/bin/env python3
import pandas as pd
import os

csv_list = "${csv_files}".split()

if not csv_list:
    print("No CSV files found to process.")
    exit(1)

all_data = []

for csv_file in csv_list:
    df = pd.read_csv(csv_file)
    all_data.append(df)

# Combine all individual CSV dataframes into one master dataframe
master_df = pd.concat(all_data, ignore_index=True)

# Calculate total unique samples PER PROTEIN to establish a dynamic denominator
samples_per_protein = master_df.groupby('PROTEIN')['SAMPLE_ID'].transform('nunique')
threshold = samples_per_protein * ${params.threshold}  # Using the threshold from params, default is 0.25 (25%)

# Calculate mutation frequency per protein and position
sample_counts_per_pos = master_df.groupby(['PROTEIN', 'POSITION'])['SAMPLE_ID'].transform('nunique')

# Filter for rows that are marked as markers OR appear in >25% (default) of total samples
relevant_df = master_df[(master_df['MARKER'] == "Yes") | (sample_counts_per_pos > threshold)]

# Write the Full Report with a sheet for each protein
with pd.ExcelWriter("final_mutations_report.xlsx") as writer:
    for prot_name, group in master_df.groupby("PROTEIN"):
        group.to_excel(writer, sheet_name=str(prot_name).strip("()',"), index=False)

# Write the Filtered Report with a sheet for each protein
with pd.ExcelWriter("relevant_mutations.xlsx") as writer:
    for prot_name, group in relevant_df.groupby("PROTEIN"):
        group.to_excel(writer, sheet_name=str(prot_name).strip("()',"), index=False)
"""
}