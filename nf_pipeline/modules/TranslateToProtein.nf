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
    cp -rL ${sequences_dir} sequences_work

    python3 ${params.workDir}/${params.programs.translateSegments} \
      --sequences-dir sequences_work \
      --output-dir sequences_work

    rm -f sequences
    mv sequences_work sequences
    """
}
