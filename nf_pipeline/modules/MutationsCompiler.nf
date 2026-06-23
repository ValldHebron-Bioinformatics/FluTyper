#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsCompiler {
    // This process compiles all individual mutation CSV files into a single Excel report
    // with separate sheets for each protein and a combined sheet for all proteins.
    errorStrategy 'ignore'

    input:
    path csv_files

    output:
    path("final_mutations_report.xlsx"), emit: results

    
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
master_df['PROTEIN'] = master_df['PROTEIN'].replace('NA', "NA ")  # Ensure 'NA' is treated as a string, not as NaN

# Write the Full Report with a sheet for each protein, plus a combined sheet
with pd.ExcelWriter("final_mutations_report.xlsx") as writer:
    master_df.to_excel(writer, sheet_name="All_Proteins", index=False)
    for prot_name, group in master_df.groupby("PROTEIN"):
        group.to_excel(writer, sheet_name=str(prot_name).strip("()',"), index=False)
"""
}
