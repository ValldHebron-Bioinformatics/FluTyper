#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GenotypingNextclade } from './modules/genotyping'
include { OrganizeBySample   } from './modules/OrganizeBySample'
include { MutationsFinder     } from './modules/mutations'
include { TranslateToProtein  } from './modules/Translation'
include { SubtypeDetection    } from './modules/SubtypeDetection'
include { GetCDS              } from './modules/GetCDS'

// Flux de treball principal
workflow {
    main:

    // Crea el canal d'entrada des dels paràmetres
    SampleInput_ch = channel
    .fromPath(params.inputFasta, checkIfExists: true)
    .splitFasta(record: [id: true])
    .map { rec -> tuple(rec.id.tokenize('[|_]')[0], file(params.inputFasta)) }
    .unique { record -> record[0] }
    
    // Executa el procés amb el canal creat
    OrganizeBySample(SampleInput_ch)
    //SubtypeDetection(OrganizeBySample.out)
    //GenotypingNextclade(SampleInput_ch, SubtypeDetection.out)
    //GetCDS(OrganizeBySample.out, SubtypeDetection.out)
    //TranslateToProtein(GetCDS.out)

    // Mutacions opcional: només si es passa --mutationsSubtype
    //def mut_out = channel.empty()
    //if (params.mutationsSubtype) {
    //    MutationsFinder(SampleId_ch)
    //    mut_out = MutationsFinder.out
    //} else {
    //    log.info "MutationsFinder omès: passa --mutationsSubtype per activar-lo."
    //}

    publish:
    folder = OrganizeBySample.out
    //subtype = SubtypeDetection.out
    //res = GenotypingNextclade.out
    //prot = TranslateToProtein.out
    //mut = mut_out
}
// Bloc final de publicació de resultats
output {
    //res {
        // Usa el primer element dl tuple (sample) per crear la carpeta
    //    path { "${params.outDir}" }
    //    mode "copy"
    //}
    folder {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    //subtype {
    //    path { "${params.outDir}" }
    //    mode "copy"
    //}
    //mut {
    //    path { "${params.outDir}" }
    //    mode "copy"
    //}
}