#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GetCDS {
    errorStrategy 'ignore'

    input:
    path(sequences_dir)
    path(inferred_subtypes)

    output:
    path("CDS")

    script:
    """
    export REFERENCES="${params.protocols.${params.protocol}.resources}/CDS_references.fasta"  
    
    mkdir -p CDS

    """
    
}