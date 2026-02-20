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
    # Descarrega el dataset de referència
    nextclade dataset get --name 'community/moncla-lab/iav-h5/ha/2.3.4.4' --output-dir nextclade_dataset
    
    # Executa analisi Nextclade
    nextclade run \\
        --input-dataset nextclade_dataset \\
        --output-json nextclade_results_${sample}.json \\
        ${dirSample}/${sample}

        python /home/vhir/Desktop/FluTyper/nf_pipeline/bin/extract_clades.py \\
        nextclade_results_${sample}.json > seqid_clade_${sample}.csv
    """
}

// Flux de treball principal
workflow {
    main:
    // Crea el canal d'entrada des dels paràmetres
    input_ch = channel.of( [ params.sample, params.dirSample ] )

    // Executa el procés amb el canal creat
    GenotypingNextclade(input_ch)

    publish:
    res = GenotypingNextclade.out
}

// Bloc final de publicació de resultats
output {
    res {
        // Usa el primer element de la tupla (sample) per crear la carpeta
        path { "${params.sample}" }
        mode "copy"
    }
}