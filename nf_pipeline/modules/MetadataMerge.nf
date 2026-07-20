#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MetadataMerge {
    errorStrategy 'ignore'

    input:
    path(subtype_csv)
    path(genotyping_csv)
    path(mutations_xlsx)
    path(metadata_csv)

    output:
    path("inferred_subtypes.csv"),        emit: subtypes
    path("final_genotyping_results.csv"), emit: genotyping
    path("final_mutations_report.xlsx"),  emit: mutations

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    def find_id_col(df):
        df.columns = df.columns.str.strip()
        for col in df.columns:
            if col.lower().replace('_', '').replace(' ', '') in ('id', 'sampleid'):
                return col
        return None

    def get_season(d):
        if pd.isna(d):
            return "Unknown Season"
        week = d.isocalendar()[1]
        y = d.year if week >= 40 else d.year - 1
        return "Season " + str(y) + "-" + str(y + 1)

    # Load and prepare metadata
    meta = pd.read_csv("${metadata_csv}")
    meta_id_col = find_id_col(meta)
    meta = meta.rename(columns={meta_id_col: '_MERGE_KEY'})
    meta['_MERGE_KEY'] = meta['_MERGE_KEY'].astype(str).str.strip()

    if 'DATE' in meta.columns:
        meta['Season'] = pd.to_datetime(meta['DATE'], errors='coerce').apply(get_season)
    else:
        meta['Season'] = "Unknown Season"

    meta = meta[['_MERGE_KEY', 'Season'] + [c for c in ['AGE GROUP', 'SEX', 'LOCATION', 'ORIGINATING_LAB'] if c in meta.columns]]

    def merge_with_meta(df):
        original_id_col = find_id_col(df)
        df = df.rename(columns={original_id_col: '_MERGE_KEY'})
        df['_MERGE_KEY'] = df['_MERGE_KEY'].astype(str).str.strip()
        merged = df.merge(meta, on='_MERGE_KEY', how='left')
        return merged.rename(columns={'_MERGE_KEY': original_id_col})

    merge_with_meta(pd.read_csv("${subtype_csv}")).to_csv("inferred_subtypes.csv", index=False)
    merge_with_meta(pd.read_csv("${genotyping_csv}")).to_csv("final_genotyping_results.csv", index=False)

    sheets = pd.read_excel("${mutations_xlsx}", sheet_name=None)
    with pd.ExcelWriter("final_mutations_report.xlsx") as w:
        for name, df in sheets.items():
            merge_with_meta(df).to_excel(w, sheet_name=name, index=False)
    """
}