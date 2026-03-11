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

    // 1. INPUT & INITIAL ORGANIZATION
    SampleInput_ch = channel
        .fromPath(params.inputFasta, checkIfExists: true)
        .splitFasta(record: [id: true])
        .map { rec -> tuple(rec.id.tokenize('[|_]')[0], file(params.inputFasta)) }
        .unique { record -> record[0] }
    
    OrganizeBySample(SampleInput_ch)


    // 2. SUBTYPE DETECTION
    SubtypeInput_ch = OrganizeBySample.out.map { sample_id, sample_dir ->
        tuple(
            sample_id,
            file("${sample_dir}/segments/HA/${sample_id}_HA.fasta"),
            file("${sample_dir}/segments/NA/${sample_id}_NA.fasta")
        )
    }

    SubtypeDetection(SubtypeInput_ch)

    // Parse subtyping results immediately for use in downstream filtering
    GenotypingInfo_ch = SubtypeDetection.out.map { sample_id, tsv_file ->
        def line = tsv_file.readLines()[0]
        def parts = line.split('\t')
        def full_subtype = parts[1] 
        def pathotype = parts.size() > 2 ? parts[2] : ""

        def h_match = (full_subtype =~ /H\d+/)
        def n_match = (full_subtype =~ /N\d+/)
        def h_tag = h_match ? h_match[0] : "Hx"
        def n_tag = n_match ? n_match[0] : "Nx"

        return tuple(sample_id, h_tag, n_tag, pathotype)
    }

    // Merge individual subtype files into a global report
    SubtypeMerged_ch = SubtypeDetection.out
        .map { arr -> arr[1] }
        .collectFile(
            name: 'inferred_subtypes.tsv',
            seed: 'seqName\tinferred_subtype\tpathotype\n',
            storeDir: "${launchDir}/${params.outDir}",
        )


    // 3. DATASET PREPARATION
    // GetDatasets depends on the merged list to know which H-types to download
    GetDatasets(SubtypeMerged_ch)


    // 4. GENOTYPING ANALYSIS (NEXTCLADE)
    GenotypingHfile_ch = OrganizeBySample.out.map { sample_id, sample_dir -> 
        tuple(sample_id, file("${sample_dir}/segments/HA/${sample_id}_HA.fasta"))
    }

    GenotypingNextcladeInput_ch = GenotypingHfile_ch
        .join(GenotypingInfo_ch)
        .combine(GetDatasets.out.flatten())
        .filter { _sample_id, _input_fasta, h_tag, _n_tag, _pathotype, dataset_dir -> 
            dataset_dir.name.contains(h_tag) 
        }

    GenotypingNextclade(GenotypingNextcladeInput_ch)

    GenotypingMerged_ch = GenotypingNextclade.out
        .collectFile(
            name: 'genotyping_results.csv',
            keepHeader: true,
            storeDir: "${launchDir}/${params.outDir}"
        )


    // 5. FINAL REPORTING & CDS EXTRACTION
    // Re-associate Nextclade files with their IDs for the final join
    NextcladeTuple_ch = GenotypingNextclade.out
        .map { csv_file ->
            def id = csv_file.name.replace('nextclade_results_', '').replace('.csv', '')
            return tuple(id, csv_file)
        }

    GenotypingResultsInput_ch = GenotypingInfo_ch
        .join(NextcladeTuple_ch, remainder: true)
        .map { row ->
            def sample_id = row[0]
            def h_tag     = row[1]
            def n_tag     = row[2]
            def pathotype = row[3]
            def csv_path  = (row.size() > 4 && row[4] != null) ? row[4] : []
            return tuple(sample_id, h_tag, n_tag, pathotype, csv_path)
        }

    GenotypingResults(GenotypingResultsInput_ch, GetDatasets.out.collect())
    
    GenotypingFinal_ch = GenotypingResults.out
        .collectFile(
            name: 'final_genotyping_results.csv',
            keepHeader: true,
            storeDir: "${launchDir}/${params.outDir}"
        )

    // Prepare inputs for sequence extraction
    CDSInput_ch = GenotypingInfo_ch
        .join(OrganizeBySample.out)
        .map { sample_id, h_tag, n_tag, pathotype, sample_dir ->
            return tuple(h_tag, n_tag, sample_id, pathotype, sample_dir)
        }

    GetCDS(CDSInput_ch)
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
    results = GenotypingFinal_ch
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
    CDS {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    //mut {
    //    path { "${params.outDir}" }
    //    mode "copy"
    //}
}