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
    df = pd.read_csv(csv_file, keep_default_na=False)
    all_data.append(df)

# Combine all individual CSV dataframes into one master dataframe
master_df = pd.concat(all_data, ignore_index=True)

# Calculate total unique samples per protein
samples_per_protein = master_df.groupby('PROTEIN')['SAMPLE_ID'].transform('nunique')
threshold = samples_per_protein * ${params.threshold}  # Using the threshold from params

# Calculate mutation frequency per protein, position, AND specific amino acid mutation
sample_counts_per_mutation = master_df.groupby(['PROTEIN', 'POSITION_REF', 'AA_MUTATION'])['SAMPLE_ID'].transform('nunique')

# Filter for rows that are marked as markers OR appear in more than the threshold of their specific protein's samples
relevant_df = master_df[(master_df['MARKER'] == "Yes") | (sample_counts_per_mutation > threshold)]

# Write the Full Report with a sheet for each protein, plus a combined sheet
with pd.ExcelWriter("final_mutations_report.xlsx") as writer:
    master_df.to_excel(writer, sheet_name="All_Proteins", index=False)
    for prot_name, group in master_df.groupby("PROTEIN"):
        group.to_excel(writer, sheet_name=str(prot_name).strip("()',"), index=False)

# Write the Filtered Report with a sheet for each protein, plus a combined sheet
with pd.ExcelWriter("relevant_mutations.xlsx") as writer:
    relevant_df.to_excel(writer, sheet_name="All_Proteins", index=False)
    for prot_name, group in relevant_df.groupby("PROTEIN"):
        group.to_excel(writer, sheet_name=str(prot_name).strip("()',"), index=False)
"""
}