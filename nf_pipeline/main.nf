#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OrganizeBySample    } from './modules/OrganizeBySample'
include { SubtypeDetection    } from './modules/SubtypeDetection'
include { GetDatasets         } from './modules/GetDatasets'
include { GenotypingNextclade } from './modules/GenotypingNextclade'
include { GenotypingResults   } from './modules/GenotypingResults'
include { GetCDS              } from './modules/GetCDS'
include { TranslateToProtein  } from './modules/TranslateToProtein'
include { MutationsFinder     } from './modules/MutationsFinder'

// Flux de treball principal
workflow {
    main:
    // INPUT & INITIAL FOLDER ORGANIZATION
    SampleInput_ch = channel
        .fromPath(params.inputFasta, checkIfExists: true)
        .splitFasta(record: [id: true])
        .map { rec -> tuple(rec.id.tokenize('[|_]')[0], file(params.inputFasta)) } // ASK ALEJANDRA: IS THIS THE MOST EFFICIENT WAY?
        .unique { record -> record[0] }
    
    OrganizeBySample(SampleInput_ch)

    // SUBTYPE DETECTION
    SubtypeInput_ch = OrganizeBySample.out.map { sample_id, sample_dir ->
        def ha_fasta = file("${sample_dir}/segments/HA/${sample_id}_HA.fasta")
        def na_fasta = file("${sample_dir}/segments/NA/${sample_id}_NA.fasta")
        tuple(sample_id, ha_fasta, na_fasta)
    }

    SubtypeDetection(SubtypeInput_ch)

    // Merge individual subtype files into a global report
    SubtypeMerged_ch = SubtypeDetection.out
        .map { tup -> tup[1] }
        .collectFile(
            name: 'inferred_subtypes.csv',
            seed: 'seqName,inferred_subtype,pathotype\n', // Add header to the merged CSV
            storeDir: "${launchDir}/${params.outDir}",
        )

    // Parse subtyping results immediately for use in downstream filtering
    GenotypingInfo_ch = SubtypeDetection.out
        .splitCsv()
        .map { sample_id, row ->
            def h_tag = (row[1] =~ /H\d+/) ? (row[1] =~ /H\d+/)[0] : "Hx"
            def n_tag = (row[1] =~ /N\d+/) ? (row[1] =~ /N\d+/)[0] : "Nx"
            def pathotype = row.size() > 2 ? row[2] : ""

            tuple(sample_id, h_tag, n_tag, pathotype)
        }

    // DATASET PREPARATION
    // GetDatasets depends on the merged list to know which H-types to download, 
    // this way it is only run once and not per sample
    GetDatasets(SubtypeMerged_ch)

    // GENOTYPING ANALYSIS (NEXTCLADE)
    GenotypingHfile_ch = SubtypeInput_ch.map { sample_id, ha_fasta, _na_fasta -> 
        tuple(sample_id, ha_fasta) 
        }
    // Join with subtyping info to filter datasets based on H-type
    GenotypingNextcladeInput_ch = GenotypingHfile_ch
        .join(GenotypingInfo_ch)
        .combine(GetDatasets.out.flatten())
        .filter { _sample_id, _input_fasta, h_tag, _n_tag, _pathotype, dataset_dir -> 
            dataset_dir.name.contains(h_tag) // Only keep datasets that match the H-type of the sample
        }

    GenotypingNextclade(GenotypingNextcladeInput_ch)
    // I think this is not needed, ASK ALEJANDRA
    // Merge individual genotyping results into a global report
    //GenotypingMerged_ch = GenotypingNextclade.out 
    //    .collectFile(
    //        name: 'genotyping_results.csv',
    //        keepHeader: true,
    //    )

    // RESULTS REPORTING & CDS EXTRACTION
    // Re-associate Nextclade files with their IDs for the final join
    NextcladeTuple_ch = GenotypingNextclade.out
        .map { csv_file ->
            def id = csv_file.name.replace('nextclade_results_', '').replace('.csv', '') // Extract sample ID from filename
            return tuple(id, csv_file)
        }
    // Join genotyping info with Nextclade results to prepare the final report
    GenotypingResultsInput_ch = GenotypingInfo_ch
        .join(NextcladeTuple_ch, remainder: true) // remainder: true to keep samples without Nextclade results (e.g. due to missing datasets) for reporting as well
        .map { sample_id, h_tag, n_tag, pathotype, csv_file -> 
            tuple(sample_id, h_tag, n_tag, pathotype, csv_file ?: [])
        }
    // GenotypingResults will need to handle cases where csv_path is empty (no Nextclade result) and report accordingly
    GenotypingResults(GenotypingResultsInput_ch, GetDatasets.out.collect()) // .collect() is key to pass the full list of datasets to each sample so no error when trying to access the dataset dir for a sample with no valid H subtype
    // Merge final genotyping results 
    GenotypingFinal_ch = GenotypingResults.out
        .collectFile(
            name: 'final_genotyping_results.csv',
            keepHeader: true,
        )

    // Prepare inputs for sequence extraction
    CDSInput_ch = GenotypingInfo_ch
        .join(OrganizeBySample.out)
        .map { sample_id, h_tag, n_tag, pathotype, sample_dir ->
            tuple(h_tag, n_tag, sample_id, pathotype, sample_dir)
        }

    GetCDS(CDSInput_ch)
    TranslateToProtein_ch = GetCDS.out
        .join(OrganizeBySample.out)
        .map { sample_id, cds_files, sample_dir ->
            tuple(sample_id, cds_files, sample_dir)
        }

    
    TranslateToProtein(TranslateToProtein_ch)

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
    //genotyping = GenotypingMerged_ch
    results = GenotypingFinal_ch
    CDS = GetCDS.out
    prot = TranslateToProtein.out
    //mut = mut_out
}
// Bloc final de publicació de resultats
output {
    //genotyping {
    //    path { "${launchDir}/${params.outDir}" }
    //    mode "copy"
    //}
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
    prot {
        path { "${launchDir}/${params.outDir}" }
        mode "copy"
    }
    //mut {
    //    path { "${params.outDir}" }
    //    mode "copy"
    //}
}