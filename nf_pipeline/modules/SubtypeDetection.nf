#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SubtypeDetection {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(ha_fasta), path(na_fasta)

    output:
    tuple val(sample_id), path("inferred_subtypes_${sample_id}.csv"), emit: results
    tuple val(sample_id), path("SDerrors.log"), optional: true, emit: errors
    
    script:
    """
    input_fasta="${sample_id}_HA_NA.fasta"
    
    # Combine existing files and ignore missing ones
    cat ${ha_fasta} ${na_fasta} 2>/dev/null > "\${input_fasta}" || true

    if [[ ! -s "\${input_fasta}" ]]; then
        echo "${sample_id},HxNx," > inferred_subtypes_${sample_id}.csv
        echo "Missing sequences for ${sample_id}" > SDerrors.log
        exit 0
    fi
         
    minimizer_index="${params.protocols[params.protocol].resources}/${params.protocol}_minimizers.json"
    nextclade sort -m "\${minimizer_index}" -r min.tsv "\${input_fasta}"

    h_tag=\$(grep '[_|]HA' min.tsv | cut -f3 | grep -oE 'H[0-9]+' | head -n 1 || true)
    n_tag=\$(grep '[_|]NA' min.tsv | cut -f3 | grep -oE 'N[0-9]+' | head -n 1 || true)
    
    pathotype=""
    if [[ "\$h_tag" =~ ^(H5|H7)\$ ]]; then
        pathotype=\$(grep '^0\t' min.tsv | cut -f3 | grep -oE "HPAI|LPAI" | head -n 1 || true)
    else
        pathotype=""
    fi

    if [[ -n "\$h_tag" && -n "\$n_tag" ]]; then subtype="\$h_tag\$n_tag"
    elif [[ -n "\$h_tag" ]]; then subtype="\${h_tag}Nx"; echo "SubtypeDetection: Missing N for ${sample_id}" > SDerrors.log
    elif [[ -n "\$n_tag" ]]; then subtype="Hx\${n_tag}"; echo "SubtypeDetection: Missing H for ${sample_id}" > SDerrors.log
    else subtype="Incomplete"; fi

    if [[ "${params.protocol.toUpperCase()}" == "HUMAN" ]]; then
        # Check if subtype matches any of the valid human combinations using regex
        if [[ ! "\$subtype" =~ ^(H1N1|H3N2|H1Nx|H3Nx|HxN1|HxN2|Incomplete)\$ ]]; then
            echo "SubtypeDetection: Unrecognized subtype for the human protocol: \${subtype} for ${sample_id}. Reclassified as Unknown." >> SDerrors.log
            subtype="Unknown"
        fi
    fi

    echo "${sample_id},\${subtype},\${pathotype}" > inferred_subtypes_${sample_id}.csv
    """
}