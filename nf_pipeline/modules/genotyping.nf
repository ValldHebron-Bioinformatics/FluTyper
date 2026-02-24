#!/usr/bin/env nextflow

// Activa la sintaxi DSL2
nextflow.enable.dsl=2

// Defineix paràmetres per defecte, es poden sobreescriure desd de la CLI
params {
    sample: String = "*.fasta" // Nom del fitxer d'entrada per defecte, es pot sobreescriure amb --sample
    dirSample: Path = "docs"
}

// Procés per analitzar clades amb Nextclade
process GenotypingNextclade {
    errorStrategy 'ignore' // Ignora errors i continua

    input:
    tuple val(sample), path(dirSample) // Rep nom del fitxer i directori

    output:
    tuple val(sample), path("seqid_clade_${sample}.csv") // Output resultant

    script:
    """
    # Filtra el FASTA original per quedar-se només amb les capçaleres que contenen |HA| o (HA)
    # El format de l'awk busca les línies que comencen per '>' i contenen el patró
    awk '/^>/ {f=(\$0 ~ /\\|HA\\|/ || \$0 ~ /\\(HA\\)/)} f' ${dirSample}/${sample} > filtered_HA.fasta

    # Descarrega el dataset de referència
    nextclade dataset get --name 'community/moncla-lab/iav-h5/ha/2.3.4.4' --output-dir nextclade_dataset
    
    # Executa analisi Nextclade
    nextclade run \
        --input-dataset nextclade_dataset \
        --output-json nextclade_results_${sample}.json \
        filtered_HA.fasta

        # Això és extra, només em serveix ara per verificar que l'assignació de clades de Nextclade és fiable, s'acabarà eliminant.
        python $projectDir/${params.programs.extractClades} \
        nextclade_results_${sample}.json > seqid_clade_${sample}.csv
    """
}



