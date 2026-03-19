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
include { MutationsCompiler   } from './modules/MutationsCompiler'
include { CompileErrors       } from './modules/CompileErrors'

workflow {
    main:
    // PRE-RUN PROTOCOL VALIDATION
    if (params.protocol == "SWINE") {
        exit 1, "PROTOCOL ERROR: The SWINE protocol is currently under development and cannot be used."
    } else if (params.protocol != "AVIAN") {
        def available = params.protocols.keySet() // List available protocols from the config for the error message, just in case we add more later.
        exit 1, "PROTOCOL ERROR: Invalid protocol specified ('${params.protocol}'). Available protocols are: ${available}."
    }
    // INPUT & INITIAL FOLDER ORGANIZATION
    SampleInput_ch = channel
        .fromPath(params.inputFasta, checkIfExists: true)
        .splitFasta(record: [id: true])
        .map { rec -> rec.id.tokenize('[|_]')[0] } 
        .unique()
 
    OrganizeBySample(SampleInput_ch)

    // SUBTYPE DETECTION
    SubtypeInput_ch = OrganizeBySample.out.results.map { sample_id, sample_dir ->
        def ha_fasta = file("${sample_dir}/segments/HA/${sample_id}_HA.fasta")
        def na_fasta = file("${sample_dir}/segments/NA/${sample_id}_NA.fasta")
        tuple(sample_id, ha_fasta, na_fasta)
    }

    SubtypeDetection(SubtypeInput_ch)

    // Merge individual subtype files into a global report
    SubtypeMerged_ch = SubtypeDetection.out.results
        .map { tup -> tup[1] }
        .collectFile(
            name: 'inferred_subtypes.csv',
            seed: 'seqName,inferred_subtype,pathotype\n' // Add header to the merged CSV
        )

    // DATASET PREPARATION
    // GetDatasets depends on the merged list to know which H-types to download, 
    // this way it is only run once and not per sample
    GetDatasets(SubtypeMerged_ch)
  
    // GENOTYPING ANALYSIS (NEXTCLADE)
    // Parse subtyping results immediately for use in downstream filtering
    GenotypingInfo_ch = SubtypeDetection.out.results
    .splitCsv()
    .map { sample_id, row ->
        def subtype = row[1]
        def pathotype = row[2]
        def h_tag = subtype.find(/H\d+/) ?: "Hx" // .find() is a groovy method similar to grep -oE
        def n_tag = subtype.find(/N\d+/) ?: "Nx"
        tuple(sample_id, h_tag, n_tag, pathotype)
    }
    GenotypingHfile_ch = SubtypeInput_ch.map { sample_id, ha_fasta, _na_fasta -> 
        tuple(sample_id, ha_fasta) 
        }
    // Join with subtyping info to filter datasets based on H-type
    GenotypingNextcladeInput_ch = GenotypingHfile_ch
        .join(GenotypingInfo_ch)
        .combine(GetDatasets.out.flatMap { datasets -> datasets }) // Spills the single list into individual items for one-on-one filtering. flatMap is best practice according to Nextflow docs
        .filter { _sample_id, _input_fasta, h_tag, _n_tag, _pathotype, dataset_dir -> 
            dataset_dir.name.contains(h_tag) 
        }

    GenotypingNextclade(GenotypingNextcladeInput_ch)
    
    // RESULTS REPORTING & CDS EXTRACTION
    // Join genotyping info with Nextclade results to prepare the final report
    GenotypingResultsInput_ch = GenotypingInfo_ch
        .join(GenotypingNextclade.out.results, remainder: true) // remainder: true to keep samples without Nextclade results (e.g. due to missing datasets) for reporting as well
        .map { sample_id, h_tag, n_tag, pathotype, csv_file -> 
            tuple(sample_id, h_tag, n_tag, pathotype, csv_file ?: [])
        }
    // GenotypingResults will need to handle cases where csv_path is empty (no Nextclade result) and report accordingly
    GenotypingResults(GenotypingResultsInput_ch, GetDatasets.out.collect()) // .collect() is key to pass the full list of datasets to each sample so no error when trying to access the dataset dir for a sample with no valid H subtype
    // Merge final genotyping results 
    GenotypingFinal_ch = GenotypingResults.out.results.map { tup -> tup[1] }
        .collectFile(
            name: 'final_genotyping_results.csv',
            keepHeader: true,
        )

    // CDS EXTRACTION & TRANSLATION
    // Prepare inputs for sequence extraction
    CDSInput_ch = GenotypingInfo_ch
        .join(OrganizeBySample.out.results)
        .map { sample_id, h_tag, n_tag, pathotype, sample_dir ->
            tuple(h_tag, n_tag, sample_id, pathotype, sample_dir)
        }
    
    GetCDS(CDSInput_ch)
    TranslateToProtein(GetCDS.out.results)

    // MUTATION IDENTIFICATION
    // Join translated protein files with genotyping info to prepare for mutation finding
    Mutations_ch = TranslateToProtein.out.results
        .join(GenotypingInfo_ch)
        .map { sample_id, prot_files, h_tag, n_tag, pathotype ->
            tuple(sample_id, prot_files, h_tag, n_tag, pathotype)
        }
    MutationsFinder(Mutations_ch)
    MutationsCompiler_ch = MutationsFinder.out.results
        .map { _sample_id, _mut_files, combined_csv -> combined_csv }
        .collect()
    MutationsCompiler(MutationsCompiler_ch)
        
    
    
    
    // Funnel all optional error channels together, then group by sample_id
    Errors_ch = OrganizeBySample.out.errors
        .mix(
            SubtypeDetection.out.errors,
            GenotypingNextclade.out.errors,
            GenotypingResults.out.errors,
            GetCDS.out.errors,
            TranslateToProtein.out.errors,
            MutationsFinder.out.errors
        )
        .groupTuple()

    // Pass the grouped bundles to a final concatenation process
    CompileErrors(Errors_ch)
    ErrorsMerged_ch = CompileErrors.out
        .map { sample_id, log_file ->
            // Read the raw text from the individual log file
            def content = log_file.text
            
            // Return a formatted string with the sample ID header
            return "========================================\n" +
                   "Errors for Sample: ${sample_id}\n" +
                   "========================================\n" +
                   "${content}\n"
        }
        .collectFile(
            name: 'pipeline_errors.log',
        )   
           
    publish:
    folder = OrganizeBySample.out.results
    subtype = SubtypeMerged_ch
    datasets = GetDatasets.out
    results = GenotypingFinal_ch
    CDS = GetCDS.out.results
    prot = TranslateToProtein.out.results
    mut = MutationsFinder.out.results
    mutations_report = MutationsCompiler.out.results
    errors = CompileErrors.out
    errors_merged = ErrorsMerged_ch
}
// Bloc final de publicació de resultats
output {
    datasets {
        path { "${projectDir}/../protocols/${params.protocol}/v1/resources" }
        mode "copy"
    }
    folder {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    
    subtype {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    results {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    CDS {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    prot {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    mut {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    mutations_report {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    errors {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    errors_merged {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
}