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
    rm -rf ./_seq_input
    cp -rL ${sequences_dir} ./_seq_input
    rm -rf ./sequences
    mv ./_seq_input ./sequences

    python3 ${params.workDir}/${params.programs.translateSegments} \
      --sequences-dir ./sequences \
      --output-dir ./sequences \
      --protocol ${params.protocol}
    """
}
