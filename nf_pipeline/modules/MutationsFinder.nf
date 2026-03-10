#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MutationsFinder {
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(dirSample)

    output:
    tuple val(sample), path("mutations_translated.csv")

    script:
    """
    cd ${params.workDir}/..

    python3 ${params.workDir}/${params.programs.mutationsDictionary} \
      --subtype ${params.mutationsSubtype ?: ''} \
      --output ${task.workDir}/mutations_translated.csv
    """
}
