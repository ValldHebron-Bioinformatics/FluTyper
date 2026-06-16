#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SubtypeDetection {
    // This process detects the subtype of influenza samples based on the provided HA and NA FASTA files.
    // To do so it uses Nextclade with the appropriate minimizer index for the specified protocol.
    // The subtype is inferred from the H and N tags extracted from the Nextclade output, and a pathotype is determined for H5 and H7 subtypes.
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
    
    # Each protocol has its own minimizer index for Nextclade, which is used to sort the sequences and extract H and N tags.
    minimizer_index="${params.protocols[params.protocol].resources}/${params.protocol}_minimizers.json"
    nextclade sort -m "\${minimizer_index}" -r min.tsv "\${input_fasta}"

    # Extract H and N tags from the Nextclade output, and determine the pathotype for H5 and H7 subtypes.
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
        if [[ "\$subtype" == "H1N1" ]]; then
            subtype="A(H1N1)pdm09"
        elif [[ "\$subtype" == "H3N2" ]]; then
            subtype="A(H3N2)"
        else
            subtype="A(\${subtype})"
        fi
    fi

    echo "${sample_id},\${subtype},\${pathotype}" > inferred_subtypes_${sample_id}.csv
    """
}
