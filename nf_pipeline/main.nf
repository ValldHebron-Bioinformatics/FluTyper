#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OrganizeBySample          } from './modules/OrganizeBySample'
include { SubtypeDetection          } from './modules/SubtypeDetection'
include { GetDatasets               } from './modules/GetDatasets'
include { FluMutDB                  } from './modules/FluMutDB'
include { MarkersFiles              } from './modules/MarkersFiles'
include { GenotypingNextclade       } from './modules/GenotypingNextclade'
include { GenotypingResults         } from './modules/GenotypingResults'
include { GetCDS                    } from './modules/GetCDS'
include { TranslateToProtein        } from './modules/TranslateToProtein'
include { MutationsFinder           } from './modules/MutationsFinder'
include { MutationsCompiler         } from './modules/MutationsCompiler'
include { CompileErrors             } from './modules/CompileErrors'
include { CladeGraphicReport        } from './modules/CladeGraphicReport'
include { MutationsGraphicReport    } from './modules/MutationsGraphicReport'
include { IndividualGraphicReport   } from './modules/IndividualGraphicReport'
include { InteractiveMutationsTable } from './modules/InteractiveMutationsTable'
include { DateGraphicReport         } from './modules/DateGraphicReport'

workflow {
    main:
    // PRE-RUN PROTOCOL VALIDATION
    if (params.protocol == "SWINE") {
        exit 1, "PROTOCOL ERROR: The SWINE protocol is currently under development and cannot be used."
    } else if (params.protocol != "AVIAN" && params.protocol != "HUMAN") {
        def available = params.protocols.keySet() 
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
        def ha_fasta = file("${sample_dir}/segments/${sample_id}_HA.fasta")
        def na_fasta = file("${sample_dir}/segments/${sample_id}_NA.fasta")
        tuple(sample_id, ha_fasta, na_fasta)
    }

    SubtypeDetection(SubtypeInput_ch)

    SubtypeMerged_ch = SubtypeDetection.out.results
        .map { tup -> tup[1] }
        .collectFile(
            name: 'inferred_subtypes.csv',
            seed: 'seqName,inferred_subtype,pathotype\n' 
        )

    // DATASET PREPARATION
    GetDatasets(SubtypeMerged_ch)

    // Initialize empty channels for downstream publish assignments
    ch_database = channel.empty()
    ch_markerfiles = channel.empty()
    ch_cds = channel.empty()
    ch_prot = channel.empty()
    ch_mut = channel.empty()
    ch_mutations_report = channel.empty()
    ch_mutations_graphic_report = channel.empty()
    ch_interactive_mutations_table = channel.empty()
    ch_individual_graphic_report = channel.empty()
    date_report_ch = channel.empty()
  
    // GENOTYPING ANALYSIS (NEXTCLADE)
    GenotypingInfo_ch = SubtypeDetection.out.results
        .splitCsv()
        .map { sample_id, row ->
            def subtype = row[1]
            def pathotype = row[2]
            def h_tag = subtype.find(/H\d+/) ?: "Hx" 
            def n_tag = subtype.find(/N\d+/) ?: "Nx"
            tuple(sample_id, h_tag, n_tag, pathotype)
        }
        
    GenotypingHfile_ch = SubtypeInput_ch.map { sample_id, ha_fasta, _na_fasta -> 
        tuple(sample_id, ha_fasta) 
    }
    
    GenotypingNextcladeInput_ch = GenotypingHfile_ch
        .join(GenotypingInfo_ch)
        .combine(GetDatasets.out.flatMap { datasets -> datasets }) 
        .filter { _sample_id, _input_fasta, h_tag, _n_tag, _pathotype, dataset_dir -> 
            dataset_dir.name.contains(h_tag) 
        }

    GenotypingNextclade(GenotypingNextcladeInput_ch)
    
    GenotypingResultsInput_ch = GenotypingInfo_ch
        .join(GenotypingNextclade.out.results, remainder: true) 
        .map { sample_id, h_tag, n_tag, pathotype, csv_file -> 
            tuple(sample_id, h_tag, n_tag, pathotype, csv_file ?: [])
        }
        
    GenotypingResults(GenotypingResultsInput_ch, GetDatasets.out.collect()) 
    
    GenotypingFinal_ch = GenotypingResults.out.results.map { tup -> tup[1] }
        .collectFile(
            name: 'final_genotyping_results.csv',
            keepHeader: true,
        )

    CladeGraphicReport(GenotypingFinal_ch)

    // MUTATIONS BLOCK CONDITIONAL EXECUTION
    if (params.protocol == "AVIAN") {
        
        FluMutDB(SubtypeMerged_ch)
        ch_database = FluMutDB.out
        
        MarkersFiles(FluMutDB.out) 
        ch_markerfiles = MarkersFiles.out

        CDSInput_ch = GenotypingInfo_ch
            .join(OrganizeBySample.out.results)
            .map { sample_id, h_tag, n_tag, pathotype, sample_dir ->
                tuple(h_tag, n_tag, sample_id, pathotype, sample_dir)
            }
        
        GetCDS(CDSInput_ch)
        ch_cds = GetCDS.out.results.map { _id, path -> path }

        TranslationInput_ch = GetCDS.out.results
            .join(GetCDS.out.aligned)
            .map { sample_id, cds_files, aligned_cds_files -> tuple(sample_id, cds_files, aligned_cds_files) }
            
        TranslateToProtein(TranslationInput_ch) 
        ch_prot = TranslateToProtein.out.results.map { _id, path -> path }

        Mutations_ch = TranslateToProtein.out.aligned
            .join(GenotypingInfo_ch)
            .map { sample_id, prot_files, h_tag, n_tag, pathotype ->
                tuple(sample_id, prot_files, h_tag, n_tag, pathotype)
            }
            
        MutationsFinder(Mutations_ch)
        ch_mut = MutationsFinder.out.results.map { _id, mut_files, combined_csv -> [mut_files, combined_csv] }.flatten()
        
        MutationsCompiler_ch = MutationsFinder.out.results
            .map { _sample_id, _mut_files, combined_csv -> combined_csv }
            .collect()
            
        MutationsCompiler(MutationsCompiler_ch)
        ch_mutations_report = MutationsCompiler.out.results

        MutationsGraphicReport(MutationsCompiler.out.results.map { full, _filtered -> full })
        ch_mutations_graphic_report = MutationsGraphicReport.out.report
        
        InteractiveMutationsTable(MutationsCompiler.out.results.map { full, _filtered -> full })
        ch_interactive_mutations_table = InteractiveMutationsTable.out.table
        
        IndividualMutations_Ch = MutationsFinder.out.results.map { sample_id, _mut_files, combined_csv -> tuple(sample_id, combined_csv) }
        IndividualGraphicReport(IndividualMutations_Ch)
        ch_individual_graphic_report = IndividualGraphicReport.out.report
        
        if (params.metadata) {
            Metadata_ch = channel.fromPath(params.metadata, checkIfExists: true)
            DateGraphicReport(MutationsCompiler.out.results.map { full, _filtered -> full }, Metadata_ch)
            date_report_ch = DateGraphicReport.out.metadata
        }
    }

    // ERROR HANDLING & COMPILATION
    BaseErrors_ch = OrganizeBySample.out.errors
        .mix(
            SubtypeDetection.out.errors,
            GenotypingNextclade.out.errors,
            GenotypingResults.out.errors
        )

    if (params.protocol == "AVIAN") {
        Errors_ch = BaseErrors_ch.mix(
            GetCDS.out.errors,
            TranslateToProtein.out.errors,
            MutationsFinder.out.errors
        ).groupTuple()
    } else {
        Errors_ch = BaseErrors_ch.groupTuple()
    }

    CompileErrors(Errors_ch)
    
    ErrorsMerged_ch = CompileErrors.out
        .map { sample_id, log_file ->
            def content = log_file.text
            return "========================================\n" +
                   "Errors for Sample: ${sample_id}\n" +
                   "========================================\n" +
                   "${content}\n"
        }
        .collectFile(
            name: 'pipeline_errors.log',
        )

    // PUBLISH DECLARATIONS
    publish:
    folder = OrganizeBySample.out.results.map { _id, path -> path }
    subtype = SubtypeMerged_ch
    datasets = GetDatasets.out
    database = ch_database
    markerfiles = ch_markerfiles
    results = GenotypingFinal_ch
    CDS = ch_cds
    prot = ch_prot
    graphic_report = CladeGraphicReport.out.report
    mutations_graphic_report = ch_mutations_graphic_report
    individual_graphic_report = ch_individual_graphic_report
    interactive_mutations_table = ch_interactive_mutations_table
    date_report = date_report_ch
    mut = ch_mut
    mutations_report = ch_mutations_report
    errors = CompileErrors.out.map { _id, log -> log }
    errors_merged = ErrorsMerged_ch
}

output {
    datasets {
        path { "${projectDir}/../protocols/${params.protocol}/v1/resources" }
        mode "copy"
    }
    database {
        path { "${projectDir}/../protocols/${params.protocol}/v1" }
        mode "copy"
    }
    markerfiles {
        path { "${projectDir}/../protocols/${params.protocol}/v1/markers" }
        mode "copy"
    }
    folder {
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
    subtype {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    results {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    mutations_report {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    mut {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    graphic_report {
        path { "${projectDir}/../${params.outDir}/graphic_reports" }
        mode "copy"
    }
    mutations_graphic_report {
        path { "${projectDir}/../${params.outDir}/graphic_reports" }
        mode "copy"
    }
    individual_graphic_report {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    interactive_mutations_table {
        path { "${projectDir}/../${params.outDir}/graphic_reports" }
        mode "copy"
    }
    date_report {
        path { "${projectDir}/../${params.outDir}/graphic_reports" }
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