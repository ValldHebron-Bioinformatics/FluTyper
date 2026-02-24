#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GenotypingNextclade } from './modules/genotyping'
include { OrganizeBySpecies   } from './modules/FolderCreation'
include { MutationsFinder     } from './modules/mutations'

// Flux de treball principal
workflow {
    main:

    // Crea el canal d'entrada des dels paràmetres
    input_ch = channel.of( [ params.sample, params.dirSample ] )

    // Executa el procés amb el canal creat
    OrganizeBySpecies(input_ch)
    GenotypingNextclade(input_ch)

    publish:
    res = GenotypingNextclade.out
    folder = OrganizeBySpecies.out
}
// Bloc final de publicació de resultats
output {
    res {
        // Usa el primer element dl tuple (sample) per crear la carpeta
        path { "${params.sample}" }
        mode "copy"
    }
    folder {
        path { "${params.sample}" }
        mode "copy"
    }
}