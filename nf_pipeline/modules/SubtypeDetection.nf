#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SubtypeDetection {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(ha_fasta), path(na_fasta)

    output:
    tuple val(sample_id), path("inferred_subtypes_${sample_id}.tsv")

    script:
    """
    input_fasta="${sample_id}_HA_NA.fasta"
    cat ${ha_fasta} ${na_fasta} > "\${input_fasta}"

    if [[ "${params.protocol}" == "AVIAN" ]]; then
        minimizer_index="${params.protocols.AVIAN.resources}/Avian_minimizers.json"
    elif [[ "${params.protocol}" == "SWINE" ]]; then
        minimizer_index="${params.protocols.SWINE.resources}/Swine_minimizers.json"
    else
        echo "No valid protocol specified for subtype detection: ${params.protocol}"
        : > minimizers_results.tsv
        printf '%s\tIncomplete\n' "${sample_id}" > inferred_subtypes_${sample_id}.tsv
        exit 0
    fi

    nextclade sort -m "\${minimizer_index}" -r minimizers_results.tsv "\${input_fasta}"

    ha_line=\$(grep -E '^0[[:space:]]' minimizers_results.tsv | head -n 1 || true)
    na_line=\$(grep -E '^1[[:space:]]' minimizers_results.tsv | head -n 1 || true)

    if [[ -z "\${ha_line}" && -z "\${na_line}" ]]; then
        printf '%s\tIncomplete\n' "${sample_id}" > inferred_subtypes_${sample_id}.tsv
        exit 0
    fi

    ha_match=\$(printf '%s' "\${ha_line}" | cut -f3)
    na_match=\$(printf '%s' "\${na_line}" | cut -f3)

    h_tag=\$(printf '%s' "\${ha_match}" | grep -oE 'H[0-9]+' | head -n 1 || true)
    n_tag=\$(printf '%s' "\${na_match}" | grep -oE 'N[0-9]+' | head -n 1 || true)
    if [[ -n "\${h_tag}" && -n "\${n_tag}" ]]; then
        subtype="\${h_tag}\${n_tag}"
    else
        subtype="Incomplete"
    fi

    printf '%s\t%s\n' "${sample_id}" "\${subtype}" > inferred_subtypes_${sample_id}.tsv

    """
}