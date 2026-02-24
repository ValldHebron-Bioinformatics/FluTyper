#!/usr/bin/env nextflow

// Activa la sintaxi DSL2
nextflow.enable.dsl=2
params {
    
}
// Procés per analitzar clades amb Nextclade
process GenotypingNextclade {
    errorStrategy 'ignore' // Ignora errors i continua

    input:
    tuple val(params.sample), path(params.dirSample) // Rep nom del fitxer i directori

    output:
    tuple val(params.sample), path("seqid_clade_${params.sample}.csv") // Output resultant

    script:
    """
    # Filtra el FASTA original per quedar-se només amb les capçaleres que contenen |HA| o (HA)
    # El format de l'awk busca les línies que comencen per '>' i contenen el patró
    awk '/^>/ {f=(\$0 ~ /\\|HA\\|/ || \$0 ~ /\\(HA\\)/)} f' ${params.dirSample}/${params.sample} > filtered_HA.fasta

    # Descarrega el dataset de referència
    nextclade dataset get --name 'community/moncla-lab/iav-h5/ha/2.3.4.4' --output-dir nextclade_dataset
    
    # Executa analisi Nextclade
    nextclade run \
        --input-dataset nextclade_dataset \
        --output-json nextclade_results_${params.sample}.json \
        filtered_HA.fasta

        # Això és extra, només em serveix ara per verificar que l'assignació de clades de Nextclade és fiable, s'acabarà eliminant.
        python ${params.workDir}/${params.programs.extractClades} \
        nextclade_results_${params.sample}.json > seqid_clade_${params.sample}.csv
    """
}



