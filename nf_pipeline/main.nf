#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GenotypingNextclade } from './modules/genotyping'
include { OrganizeBySpecies   } from './modules/FolderCreation'
include { MutationsFinder     } from './modules/mutations'
include { TranslateToProtein  } from './modules/Translation'
include { SubtypeDetection    } from './modules/SubtypeDetection'

// Flux de treball principal
workflow {
    main:

    // Crea el canal d'entrada des dels paràmetres
    input_ch = channel.of( [ params.sample, params.dirSample ] )

    // Executa el procés amb el canal creat
    OrganizeBySpecies(input_ch)
    GenotypingNextclade(input_ch)
    TranslateToProtein(OrganizeBySpecies.out)
    SubtypeDetection(input_ch)

    // Mutacions opcional: només si es passa --mutationsSubtype
    def mut_out = channel.empty()
    if (params.mutationsSubtype) {
        MutationsFinder(input_ch)
        mut_out = MutationsFinder.out
    } else {
        log.info "MutationsFinder omès: passa --mutationsSubtype per activar-lo."
    }

    publish:
    res = GenotypingNextclade.out
    folder = TranslateToProtein.out
    subtype_minimizers = SubtypeDetection.out.minimizers
    subtype_inferred = SubtypeDetection.out.inferred
    mut = mut_out
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
    subtype_minimizers {
        path { "${params.sample}" }
        mode "copy"
    }
    subtype_inferred {
        path { "${params.sample}" }
        mode "copy"
    }
    mut {
        path { "${params.sample}" }
        mode "copy"
    }
}