#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GenotypingNextclade {
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(ha_fasta), val(h_tag), val(n_tag), val(pathotype), val(dataset_dir)
    output:
    path("nextclade_results_${sample_id}.csv")
    script:
    def logDir = file(params.outDir)
    """
    # Genotyping using Nextclade with the appropriate dataset based on the H subtype
    if [[ ${h_tag} == "H5" ]]; then
        dataset_dir="${dataset_dir}/H5/nextclade_H5_dataset"
    elif [[ ${h_tag} == "H7" ]]; then
        ##dataset_dir="${dataset_dir}/H7/nextclade_H7_dataset"
        touch nextclade_results_${sample_id}.csv
        exit 0
    elif [[ ${h_tag} == "H9" ]]; then
        ##dataset_dir="${dataset_dir}/H9/nextclade_H9_dataset"
        touch nextclade_results_${sample_id}.csv
        exit 0
    else
        echo "No valid H subtype found for genotyping: ${h_tag}" >> "${logDir}/errors.log"
        touch nextclade_results_${sample_id}.csv
        exit 0 ## Theoretically, this should not happen because we only get the H tags from the inferred subtypes, but we add this check just in case.
    fi
    nextclade run \
        --input-dataset "${dataset_dir}" \
        --output-csv nextclade_results_${sample_id}.csv \
        "${ha_fasta}"
    """
}