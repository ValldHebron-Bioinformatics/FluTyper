#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GenotypingNextclade {
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(ha_fasta), val(h_tag), val(n_tag), val(pathotype), val(dataset_dir)
    
    output:
    tuple val(sample_id), path("nextclade_results_${sample_id}.csv"), emit: results
    tuple val(sample_id), path("GNerrors.log"), optional: true, emit: errors
    
    script:
    """
    # Genotyping using Nextclade with the appropriate dataset based on the H subtype
    if [ "${params.protocol}" = "HUMAN" ]; then
        if [[ ${h_tag} != "H1" && ${h_tag} != "H3" ]]; then
            echo "No valid H subtype found for HUMAN genotyping: ${h_tag}" >> GNerrors.log
            touch nextclade_results_${sample_id}.csv
            exit 0
        fi
    else
        if [[ ${h_tag} == "H7" || ${h_tag} == "H9" ]]; then
            touch nextclade_results_${sample_id}.csv
            exit 0
        elif [[ ${h_tag} != "H5" ]]; then
            echo "No valid H subtype found for AVIAN genotyping: ${h_tag}" >> GNerrors.log
            touch nextclade_results_${sample_id}.csv
            exit 0
        fi
    fi

    nextclade run \
        --input-dataset "${dataset_dir}" \
        --output-csv nextclade_results_${sample_id}.csv \
        "${ha_fasta}"
    """
}