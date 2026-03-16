#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SubtypeDetection {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(ha_fasta), path(na_fasta)

    output:
    tuple val(sample_id), path("inferred_subtypes_${sample_id}.csv")
    
    script:
    def logDir = file(params.outDir)
    """
    input_fasta="${sample_id}_HA_NA.fasta"
    cat ${ha_fasta} ${na_fasta} > "\${input_fasta}" # Combine HA and NA fastas for subtyping

    if [[ "${params.protocol}" == "AVIAN" ]]; then
        minimizer_index="${params.protocols.AVIAN.resources}/Avian_minimizers.json" 
    elif [[ "${params.protocol}" == "SWINE" ]]; then
        minimizer_index="${params.protocols.SWINE.resources}/Swine_minimizers.json"
    else
        echo "No valid protocol specified for subtype detection: ${params.protocol}" >> "${logDir}/errors.log"
        printf '%s,%s,%s\n' "${sample_id}" "Incomplete" "" > inferred_subtypes_${sample_id}.csv
        exit 1 ## If no valid protocol, program cannot proceed with subtyping, so genotyping also cannot proceed.
    fi         ## ASK ALEJANDRA if we want to exit with error or just create an empty results file and exit with 0
    
    nextclade sort -m "\${minimizer_index}" -r minimizers_results.tsv "\${input_fasta}"

    # Extract H and N tags from minimizer results, determine pathotype for H5/H7/H9
    h_tag=\$(grep -E '^[0-9]\t${sample_id}[_|]HA' minimizers_results.tsv | head -n 1 | cut -f3 | grep -oE 'H[0-9]+' | head -n 1 || true)
    n_tag=\$(grep -E '^[0-9]\t${sample_id}[_|]NA' minimizers_results.tsv | head -n 1 | cut -f3 | grep -oE 'N[0-9]+' | head -n 1 || true)
    pathotype="" # Default pathotype for non-H5/H7, will be updated later if needed
    if [[ "\$h_tag" == "H5" || "\$h_tag" == "H7" ]]; then
        pathotype=\$(grep -E '^0\t' minimizers_results.tsv | head -n 1 | cut -f3 |grep -oE "HPAI|LPAI" | head -n 1 || true)
    fi
    if [[ "\$h_tag" == "H9" ]]; then
        pathotype="LPAI"
    fi
    if [[ -n "\${h_tag}" && -n "\${n_tag}" ]]; then
        subtype="\${h_tag}\${n_tag}"
    elif [[ -n "\${h_tag}" && -z "\${n_tag}" ]]; then # If only H is detected, assign N as "Nx" to indicate unknown N subtype
        echo "N subtype not detected for sample ${sample_id}, assigning as Nx." >> "${logDir}/errors.log"
        subtype="\${h_tag}Nx"
    elif [[ -n "\${n_tag}" && -z "\${h_tag}" ]]; then # Same for H if only N is detected
        echo "H subtype not detected for sample ${sample_id}, assigning as Hx. Cannot proceed with clade inference." >> "${logDir}/errors.log"
        subtype="Hx\${n_tag}"
    else
        subtype="Incomplete"
    fi
    printf '%s,%s,%s\n' "${sample_id}" "\${subtype}" "\${pathotype}" > inferred_subtypes_${sample_id}.csv
    """
}