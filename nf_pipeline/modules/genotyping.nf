#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GenotypingNextclade {
    errorStrategy 'ignore'

    input:
    tuple val(params.sample), path(params.dirSample)

    output:
    tuple val(params.sample),
          path("seqid_clade_${params.sample}.csv"),
          path("clade_check_${params.sample}.csv")

    script:
    """
    # Filtra només capçaleres HA
    awk '/^>/ {f=(index(\$0, "|HA|") > 0 || index(\$0, "_HA_") > 0)} f' "${params.dirSample}/${params.sample}" > filtered_HA.fasta

    # Dataset Nextclade: usa local i, si no existeix, el descarrega
    DATASET_NAME='community/moncla-lab/iav-h5/ha/2.3.4.4'
    LOCAL_DATASET="${params.workDir}/../docs/nextclade_dataset"
    [ -f "\${LOCAL_DATASET}/pathogen.json" ] || nextclade dataset get --name "\$DATASET_NAME" --output-dir "\$LOCAL_DATASET"
    cp -r "\${LOCAL_DATASET}" ./nextclade_dataset

    # Executa anàlisi Nextclade
    nextclade run \
        --input-dataset nextclade_dataset \
        --output-csv nextclade_results_${params.sample}.csv \
        filtered_HA.fasta

    python ${params.workDir}/${params.programs.extractClades} \
        nextclade_results_${params.sample}.csv > clade_check_${params.sample}.csv

    python ${params.workDir}/${params.programs.createCSV} \
        nextclade_results_${params.sample}.csv seqid_clade_${params.sample}.csv 
    """
}



