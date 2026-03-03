#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TranslateToProtein {
    errorStrategy 'ignore'

    input:
    path(sequences_dir)

    output:
    path('sequences')

    script:
    """
    python3 ${params.workDir}/${params.programs.translateSegments} \
      --sequences-dir ${sequences_dir} \
      --output-dir ${sequences_dir}
    """
}
