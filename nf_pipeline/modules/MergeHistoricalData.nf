#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MergeHistoricalData {
    errorStrategy 'ignore'
    
    input:
    path "new_sub.csv"
    path "new_geno.csv"
    path "new_mut.xlsx"
    val meta_path
    path append_dir

    output:
    path "inferred_subtypes.csv",        emit: subtypes
    path "final_genotyping_results.csv", emit: genotyping
    path "final_mutations_report.xlsx",  emit: mutations
    path "metadata.csv",                 emit: metadata

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os

    # Define historical file paths
    old_sub_path = os.path.join("${append_dir}", "inferred_subtypes.csv")
    old_geno_path = os.path.join("${append_dir}", "final_genotyping_results.csv")
    old_mut_path = os.path.join("${append_dir}", "final_mutations_report.xlsx")
    old_meta_path = os.path.join("${append_dir}", "metadata.csv")

    # Merge Subtypes
    df_new_sub = pd.read_csv("new_sub.csv")
    if os.path.exists(old_sub_path):
        df_old_sub = pd.read_csv(old_sub_path)
        df_old_sub = df_old_sub[~df_old_sub['seqName'].isin(df_new_sub['seqName'])]
        df_final_sub = pd.concat([df_old_sub, df_new_sub], ignore_index=True)
    else:
        df_final_sub = df_new_sub
    df_final_sub.to_csv("inferred_subtypes.csv", index=False)

    # Merge Genotyping
    df_new_geno = pd.read_csv("new_geno.csv")
    if os.path.exists(old_geno_path):
        df_old_geno = pd.read_csv(old_geno_path)
        if 'SampleID' in df_old_geno.columns and 'SampleID' in df_new_geno.columns:
            df_old_geno = df_old_geno[~df_old_geno['SampleID'].isin(df_new_geno['SampleID'])]
        df_final_geno = pd.concat([df_old_geno, df_new_geno], ignore_index=True)
    else:
        df_final_geno = df_new_geno
    df_final_geno.to_csv("final_genotyping_results.csv", index=False)

    # Merge Mutations
    new_mut_sheets = pd.read_excel("new_mut.xlsx", sheet_name=None, keep_default_na=False)
    
    with pd.ExcelWriter("final_mutations_report.xlsx") as writer:
        if os.path.exists(old_mut_path):
            old_mut_sheets = pd.read_excel(old_mut_path, sheet_name=None, keep_default_na=False)
            
            for sheet_name, df_new in new_mut_sheets.items():
                if sheet_name in old_mut_sheets:
                    df_old = old_mut_sheets[sheet_name]
                    if 'SAMPLE_ID' in df_old.columns and 'SAMPLE_ID' in df_new.columns:
                        df_old = df_old[~df_old['SAMPLE_ID'].isin(df_new['SAMPLE_ID'])]
                    
                    df_final = pd.concat([df_old, df_new], ignore_index=True)
                    df_final.to_excel(writer, sheet_name=sheet_name, index=False)
                else:
                    df_new.to_excel(writer, sheet_name=sheet_name, index=False)
        else:
            for sheet_name, df_new in new_mut_sheets.items():
                df_new.to_excel(writer, sheet_name=sheet_name, index=False)

    # Merge Metadata
    new_meta_path = "${meta_path}"
    df_meta_list = []
    
    if os.path.exists(old_meta_path):
        try:
            df_meta_list.append(pd.read_csv(old_meta_path))
        except Exception:
            pass

    if os.path.exists(new_meta_path) and os.path.getsize(new_meta_path) > 0:
        try:
            df_meta_list.append(pd.read_csv(new_meta_path))
        except Exception:
            pass

    if df_meta_list:
        df_final_meta = pd.concat(df_meta_list, ignore_index=True)
        if 'ID' in df_final_meta.columns:
            df_final_meta = df_final_meta.drop_duplicates(subset=['ID'], keep='last')
        df_final_meta.to_csv("metadata.csv", index=False)
    else:
        pd.DataFrame(columns=['ID', 'DATE']).to_csv("metadata.csv", index=False)
    """
}