#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GenotypingNextclade {
    // This process performs genotyping of influenza samples using Nextclade and
    // if clade 2.3.4.4b is detected, it also runs genin2 for further genotyping.
    errorStrategy 'ignore'
    debug true

    input:
    tuple val(sample_id), path(ha_fasta), val(h_tag), val(n_tag), val(pathotype), val(dataset_dir), path(sample_dir)

    output:
    tuple val(sample_id), path("nextclade_results_${sample_id}.csv"), emit: results
    tuple val(sample_id), path("genin_results_${sample_id}.tsv"), optional: true, emit: genin
    tuple val(sample_id), path("GNerrors.log"), optional: true, emit: errors
    // Copy of the CSV inside the sample's own published folder (samples/<id>/nextclade_results.csv), so a later
    // --append run (with --phylogenetics) can reuse this sample's subclade instead of rerunning Nextclade for it
    path("samples/${sample_id}/nextclade_results.csv"), optional: true, emit: sample_csv

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

    mkdir -p "samples/${sample_id}"
    cp "nextclade_results_${sample_id}.csv" "samples/${sample_id}/nextclade_results.csv"

    # Genotyping with genin2 if clade 2.3.4.4b is detected in the Nextclade results
    if grep -q "2.3.4.4b" "nextclade_results_${sample_id}.csv"; then
        if [ -s "${sample_dir}/${sample_id}.fasta" ]; then
            # Adapt the header for genin2 input (replace '|' with '_' and keep only the first two fields)
            # Written to this task's own directory, not into the staged sample_dir that belongs to OrganizeBySample
            cat "${sample_dir}/${sample_id}.fasta" | tr "|" "_" | cut -d "_" -f1,2 > "./${sample_id}_genin_input.fasta"
            
            genin2 -o "genin_results_${sample_id}.tsv" "./${sample_id}_genin_input.fasta"
            
            if [ ! -f "genin_results_${sample_id}.tsv" ]; then
                echo "GenotypingNextclade: Genin2 failed to produce output for ${sample_id}" >> GNerrors.log
            fi
        else
            echo "GenotypingNextclade: Input FASTA file not found or empty for ${sample_id}" >> GNerrors.log
        fi
    fi
    """
}
