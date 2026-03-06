#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GenotypingNextclade {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(sample_fasta)
    path(inferred_subtypes)

    output:
    tuple val(sample_id),
            path("seqid_clade_${sample_id}.csv")

    script:
    """
    # Filtra només capçaleres HA
    awk '/^>/ {f=(index(\$0, "|HA|") > 0 || index(\$0, "_HA_") > 0)} f' "${sample_fasta}" > filtered_HA.fasta

    # Dataset Nextclade: usa local i, si no existeix, el descarrega
    DATASET_NAME='community/moncla-lab/iav-h5/ha/2.3.4.4'
    LOCAL_DATASET="${params.workDir}/../docs/nextclade_dataset"
    [ -f "\${LOCAL_DATASET}/pathogen.json" ] || nextclade dataset get --name "\$DATASET_NAME" --output-dir "\$LOCAL_DATASET"
    cp -r "\${LOCAL_DATASET}" ./nextclade_dataset

    # Executa anàlisi Nextclade
    nextclade run \
        --input-dataset nextclade_dataset \
        --output-csv nextclade_results_${sample_id}.csv \
        filtered_HA.fasta

    ### python ${params.workDir}/${params.programs.extractClades} \
        nextclade_results_${sample_id}.csv > clade_check_${sample_id}.csv

    python ${params.workDir}/${params.programs.createCSV} \
        nextclade_results_${sample_id}.csv seqid_clade_${sample_id}.csv 
    """
}



