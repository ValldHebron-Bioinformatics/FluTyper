#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// No hi ha absolutament res definitiu aquí, només era per poder incloure a main.nf
process MutationsFinder {
    errorStrategy 'ignore' // Ignora errors i continua

    input:
    tuple val(sample), path(dirSample) // Rep nom del fitxer i directori

    output:
    tuple val(sample), path("mutations_${sample}.csv") // Output resultant

    script:
    """
    # Descarrega el dataset de referència
    nextclade dataset get --name 'community/moncla-lab/iav-h5/ha/
    """
}
