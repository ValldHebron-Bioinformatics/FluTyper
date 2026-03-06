#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process OrganizeBySample {
    errorStrategy 'ignore' // Ignora errors i continua
    

    input:
    tuple val(sample_id), path(input_fasta)

    output:
    tuple val(sample_id), path("${sample_id}")

    script:
    """
    mkdir -p ${sample_id}

    seqkit grep -r -p ${sample_id} ${input_fasta} > ${sample_id}/${sample_id}.fasta

    # 2. per cada segment, crear carpeta i fitxer fasta
    for seg in PB2 PB1 PA HA NP NA MP NS; do
        mkdir -p ${sample_id}/segments/\${seg}
        seqkit grep -r -p "${sample_id}_\${seg}_" ${input_fasta} > ${sample_id}/segments/\${seg}/${sample_id}_\${seg}.fasta 
    done
    """
    }