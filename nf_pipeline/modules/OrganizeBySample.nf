#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process OrganizeBySample {
    errorStrategy 'ignore' // Ignora errors i continua
    
    input:
    tuple val(sample_id), path(input_fasta)

    output:
    tuple val(sample_id), path("samples/${sample_id}")

    script:
    """
    mkdir -p "samples/${sample_id}"

    # seqkit -r regex and -p pattern to extract all records for the sample.
    seqkit grep -r -p ${sample_id} ${input_fasta} > "samples/${sample_id}/${sample_id}.fasta"

    # Convert the params.segments list to a space-separated string with join(' ') for Bash loop iteration
    for seg in ${params.segments.join(' ')}; do
        mkdir -p "samples/${sample_id}/segments/\${seg}"
        seqkit grep -r -p "${sample_id}[_|]\${seg}[_|]" ${input_fasta} > "samples/${sample_id}/segments/\${seg}/${sample_id}_\${seg}.fasta"
    done
    """
    }