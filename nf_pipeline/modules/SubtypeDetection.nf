#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SubtypeDetection {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(ha_fasta), path(na_fasta)

    output:
    tuple val(sample_id), path("inferred_subtypes_${sample_id}.csv")
    
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
        printf '%s,%s,%s\n' "${sample_id}" "Incomplete" "" > inferred_subtypes_${sample_id}.csv
        exit 0
    fi

    nextclade sort -m "\${minimizer_index}" -r minimizers_results.tsv "\${input_fasta}"

    
    h_tag=\$(grep -E '^0\t' minimizers_results.tsv | head -n 1 | cut -f3 | grep -oE 'H[0-9]+' | head -n 1 || true)
    n_tag=\$(grep -E '^1\t' minimizers_results.tsv | head -n 1 | cut -f3 | grep -oE 'N[0-9]+' | head -n 1 || true)
    pathotype="" # Default pathotype for non-H5/H7, will be updated later if needed
    if [[ "\$h_tag" == "H5" || "\$h_tag" == "H7" ]]; then
        pathotype=\$(grep -E '^0\t' minimizers_results.tsv | head -n 1 | cut -f3 |grep -oE "HPAI|LPAI" | head -n 1 || true)
    fi
    if [[ "\$h_tag" == "H9" ]]; then
        pathotype="LPAI"
    fi
    if [[ -n "\${h_tag}" && -n "\${n_tag}" ]]; then
        subtype="\${h_tag}\${n_tag}"
    elif [[ -n "\${h_tag}" ]]; then
        subtype="\${h_tag}Nx"
    elif [[ -n "\${n_tag}" ]]; then
        subtype="Hx\${n_tag}"
    else
        subtype="Incomplete"
    fi
    printf '%s,%s,%s\n' "${sample_id}" "\${subtype}" "\${pathotype}" > inferred_subtypes_${sample_id}.csv
    """
}