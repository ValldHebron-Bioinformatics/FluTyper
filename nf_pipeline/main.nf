#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GenotypingNextclade } from './modules/genotyping'
include { OrganizeBySpecies   } from './modules/FolderCreation'
include { MutationsFinder     } from './modules/mutations'
include { TranslateToProtein  } from './modules/Translation'

// Flux de treball principal
workflow {
    main:

    // Valida que el subtipus final estigui definit i sigui compatible.
    def validMutationSubtypes = ['H1', 'H3', 'H5', 'H7', 'H9'] as Set
    def resolvedMutationSubtype = params.mutationsSubtype ?: (params.protocol ? params.protocol.toString().trim().toUpperCase() : null)
    if (!resolvedMutationSubtype) {
        error "No s'ha pogut resoldre 'mutationsSubtype'. Revisa --protocol (per defecte H5) o passa --mutationsSubtype."
    }
    if (!(resolvedMutationSubtype in validMutationSubtypes)) {
        error "mutationsSubtype '${resolvedMutationSubtype}' no és vàlid. Valors admesos: ${validMutationSubtypes.join(', ')}"
    }

    // Crea el canal d'entrada des dels paràmetres
    input_ch = channel.of( [ params.sample, params.dirSample ] )
    def resolvedProtocol = (params.protocol ? params.protocol.toString().trim().toUpperCase() : 'H5')

    // Executa el procés amb el canal creat
    OrganizeBySpecies(input_ch)
    GenotypingNextclade(input_ch)
    TranslateToProtein(OrganizeBySpecies.out)

    // Si el protocol és H5, s'omet la traducció de mutacions.
    def mut_out = channel.empty()
    if (resolvedProtocol != 'H5') {
        MutationsFinder(input_ch)
        mut_out = MutationsFinder.out
    }

    publish:
    res = GenotypingNextclade.out
    folder = TranslateToProtein.out
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
    mut {
        path { "${params.sample}" }
        mode "copy"
    }
}