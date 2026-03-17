#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CompileErrors {
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(logs)

    output:
    tuple val(sample_id), path("samples/${sample_id}/${sample_id}_final_errors.log")

    script:
    """
    mkdir -p "samples/${sample_id}"
    cat $logs > "samples/${sample_id}/${sample_id}_final_errors.log"
    """
}