#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SubtypeDetection {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(sample_fasta)

    output:
    tuple val(sample_id), path("inferred_subtypes.tsv")

    script:
    """
    input_fasta="${sample_fasta}"

    if [[ "${params.protocol}" == "AVIAN" ]]; then
        echo "Subtype detection for AVIAN protocol"
        minimizer_index="${params.protocols.AVIAN.resources}/Avian_minimizers.json"
    elif [[ "${params.protocol}" == "SWINE" ]]; then
        echo "Subtype detection for SWINE protocol"
        minimizer_index="${params.protocols.SWINE.resources}/Swine_minimizers.json"
    else
        echo "No valid protocol specified for subtype detection: ${params.protocol}"
        : > minimizers_results.tsv
        printf 'id\tha_match\tha_score\tna_match\tna_score\tinferred_subtype\n' > inferred_subtypes.tsv
        exit 0
    fi

    nextclade sort -m "\${minimizer_index}" -r minimizers_results.tsv "\${input_fasta}"

    python3 ${params.workDir}/${params.programs.subtypeInference}

        

    """
}