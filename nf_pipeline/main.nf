#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GenotypingNextclade } from './modules/GenotypingNextclade'
include { OrganizeBySample    } from './modules/OrganizeBySample'
include { MutationsFinder     } from './modules/MutationsFinder'
include { TranslateToProtein  } from './modules/TranslateToProtein'
include { SubtypeDetection    } from './modules/SubtypeDetection'
include { GetDatasets         } from './modules/GetDatasets'
include { GetCDS              } from './modules/GetCDS'
include { GenotypingResults   } from './modules/GenotypingResults'

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
    // Debug: view output of SubtypeDetection
    SubtypeMerged_ch = SubtypeDetection.out
        .map { arr -> arr[1] }
        // Uneix tots els TSV individuals en un únic fitxer de resultats
        .collectFile(
            // Nom del fitxer agregat final
            name: 'inferred_subtypes.tsv',
            // Header inicial que s'escriu abans del contingut recopilat
            seed: 'seqName\tinferred_subtype\n',
            // Directori on es desa el fitxer final
            storeDir: "${launchDir}/${params.outDir}",
        )

    GetDatasets(SubtypeMerged_ch)
    GenotypingHfile_ch = OrganizeBySample.out.map { sample_id, sample_dir -> 
        tuple(sample_id, file("${sample_dir}/segments/HA/${sample_id}_HA.fasta"))
    }

    // Keep the sample_id (arr[0]) so the join works, but extract the H and N tags from the inferred_subtypes.tsv content
    GenotypingTags_ch = SubtypeDetection.out.map { sample_id, tsv_file ->
        // Read the line (format: sample_id \t subtype)
        def line = tsv_file.readLines()[0]
        def full_subtype = line.split('\t')[1] // e.g., "H5N1"

        // Extract H and N parts using regex
        def h_match = (full_subtype =~ /H\d+/)
        def n_match = (full_subtype =~ /N\d+/)

        def h_tag = h_match ? h_match[0] : "Hx"
        def n_tag = n_match ? n_match[0] : "Nx"

        return tuple(sample_id, h_tag, n_tag)
    }

    GenotypingInput_ch = GenotypingHfile_ch
        .join(GenotypingTags_ch) // Joins on sample_id
        .combine(GetDatasets.out.flatten())
        .map { sample_id, input_fasta, h_tag, n_tag, dataset_dir ->
            // Filter: Only proceed if the dataset directory matches the H-tag
            if (dataset_dir.name.contains(h_tag)) {
                return tuple(sample_id, input_fasta, h_tag, n_tag, dataset_dir)
            }
        }

    GenotypingNextclade(GenotypingInput_ch)
    GenotypingMerged_ch = GenotypingNextclade.out
        .map { arr -> arr[1] } // Extract the path to the CSV file from the tuple
        .collectFile(
            name: 'genotyping_results.csv',
            keepHeader: true,
            skip: 1,
            storeDir: "${launchDir}/${params.outDir}"
        )
      
    GenotypingResults(GenotypingNextclade.out, GetDatasets.out)
    
    
    // Unim els canals una sola vegada i creem una tupla neta
    //GetCDS_ch = OrganizeBySample.out.join(SubtypeDetection.out)
    //    .map { sample_id, sample_dir, subtype_file ->
    //        def segments_dir = file("${sample_dir}/segments")
    //        tuple(sample_id, segments_dir, subtype_file)
    //    }
    //
    //// Passem el canal sencer al procés
    //GetCDS(GetCDS_ch)
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
    datasets = GetDatasets.out
    genotyping = GenotypingMerged_ch
    results = GenotypingResults.out
    //CDS = GetCDS.out
    //prot = TranslateToProtein.out
    //mut = mut_out
}
// Bloc final de publicació de resultats
output {
    genotyping {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    datasets {
        path { "${launchDir}/protocols/${params.protocol}/v1/resources" }
        mode "copy"
    }
    folder {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    
    subtype {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    results {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    //CDS {
    //    path { "${launchDir}/${params.outDir}" }
    //    mode "copy"
    //}
    //mut {
    //    path { "${params.outDir}" }
    //    mode "copy"
    //}
}