#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GenotypingNextclade } from './modules/GenotypingNextclade'
include { OrganizeBySample   } from './modules/OrganizeBySample'
include { MutationsFinder     } from './modules/MutationsFinder'
include { TranslateToProtein  } from './modules/TranslateToProtein'
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
    SubtypeInput_ch = OrganizeBySample.out.map { sample_id, sample_dir ->
        tuple(
            sample_id,
            file("${sample_dir}/segments/HA/${sample_id}_HA.fasta"),
            file("${sample_dir}/segments/NA/${sample_id}_NA.fasta")
        )
    }
    SubtypeDetection(SubtypeInput_ch)

    // Agafa la sortida del procés, elimina el sample_id i conserva només el fitxer TSV
    // de cada mostra per poder-los fusionar en un únic fitxer final
    SubtypeMerged_ch = SubtypeDetection.out
        .map { _sample_id, subtype_file -> subtype_file }
        // Uneix tots els TSV individuals en un únic fitxer de resultats
        .collectFile(
            // Nom del fitxer agregat final
            name: 'inferred_subtypes.tsv',
            // Header inicial que s'escriu abans del contingut recopilat
            seed: 'seqName\tinferred_subtype\n',
            // Directori on es desa el fitxer final
            storeDir: "${launchDir}/${params.outDir}",
            newLine: false
        )
        .first()

    GenotypingNextclade(SampleInput_ch, SubtypeMerged_ch)
    
    GetCDS(OrganizeBySample.out, SubtypeDetection.out)
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
    subtype = SubtypeMerged_ch
    genotyping = GenotypingNextclade.out
    CDS = GetCDS.out
    //prot = TranslateToProtein.out
    //mut = mut_out
}
// Bloc final de publicació de resultats
output {
    genotyping {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    folder {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    //subtype {
    subtype {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    //mut {
    //    path { "${params.outDir}" }
    //    mode "copy"
    //}
}