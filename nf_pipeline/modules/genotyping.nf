#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GenotypingNextclade {
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(dirSample)

    output:
    tuple val(sample), path("seqid_clade_${sample}.tsv")

    script:
    """
    nextclade dataset get --name 'community/moncla-lab/iav-h5/ha/all-clades' --output-dir nextclade_dataset
    
    nextclade run \\
        --input-dataset nextclade_dataset \\
        --output-json nextclade_results_${sample}.json \\
        ${dirSample}/${sample}

    python /home/vhir/Desktop/FluTyper/nf_pipeline/bin/extract_clades.py \\
        nextclade_results_${sample}.json > seqid_clade_${sample}.tsv
    """
}