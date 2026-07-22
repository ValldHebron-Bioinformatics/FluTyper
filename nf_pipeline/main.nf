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
include { MergeHistoricalData       } from './modules/MergeHistoricalData'
include { GeographicReport          } from './modules/GeographicReport'
include { MetadataMerge             } from './modules/MetadataMerge'

// Comprovació invisible per decidir si generem el mapa
def check_location_column(metadata_path) {
    if (!metadata_path) return false
    def f = file(metadata_path)
    if (!f.exists()) return false
    
    def has_loc = false
    f.withReader { reader ->
        def header = reader.readLine()
        if (header && header.toUpperCase().contains("LOCATION")) {
            has_loc = true
        }
    }
    return has_loc
}

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
        .map { rec -> rec.id.tokenize('[|_]')[0] } // Extract sample ID from FASTA header using the first token before '|' or '_'
        .unique()
 
    OrganizeBySample(SampleInput_ch)

    // SUBTYPE DETECTION
    // Prepare channel for subtype detection by mapping sample IDs to their corresponding HA and NA FASTA files
    SubtypeInput_ch = OrganizeBySample.out.results.map { sample_id, sample_dir ->
        def ha_fasta = file("${sample_dir}/segments/${sample_id}_HA.fasta")
        def na_fasta = file("${sample_dir}/segments/${sample_id}_NA.fasta")
        tuple(sample_id, ha_fasta, na_fasta)
    }

    SubtypeDetection(SubtypeInput_ch)

    // Collect inferred subtypes into a single CSV file for downstream processing
    SubtypeMerged_ch = SubtypeDetection.out.results
        .map { tup -> tup[1] }
        .collectFile(
            name: 'inferred_subtypes.csv',
            seed: 'Sample_ID,inferred_subtype,pathotype\n' 
            seed: 'Sample_ID,inferred_subtype,pathotype\n' 
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
    ch_clade_evolution_report = channel.empty()
    date_report_ch = channel.empty()
    ch_geo_report = channel.empty()
    ch_geo_report = channel.empty()
  
    // GENOTYPING ANALYSIS (NEXTCLADE)
    GenotypingInfo_ch = SubtypeDetection.out.results
        .splitCsv()
        .map { sample_id, row ->
            def subtype = row[1]
            def pathotype = row[2]
            def h_tag = subtype.find(/H\d+/) ?: "Hx" // Extract H subtype or default to "Hx" if not found
            def n_tag = subtype.find(/N\d+/) ?: "Nx" // Extract N subtype or default to "Nx" if not found
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
        .join(OrganizeBySample.out.results)
        .map { sample_id, ha_fasta, h_tag, n_tag, pathotype, dataset_dir, sample_dir ->
            tuple(sample_id, ha_fasta, h_tag, n_tag, pathotype, dataset_dir, sample_dir)
        }
    GenotypingNextclade(GenotypingNextcladeInput_ch)
    
    GenotypingResultsInput_ch = GenotypingInfo_ch
    .join(GenotypingNextclade.out.results, remainder: true) 
    .join(GenotypingNextclade.out.genin, remainder: true) 
    .map { sample_id, h_tag, n_tag, pathotype, csv_file, genin_file -> 
        tuple(sample_id, h_tag, n_tag, pathotype, csv_file ?: [], genin_file ?: []) // Ensure that missing files are represented as empty lists
    }
        
    GenotypingResults(GenotypingResultsInput_ch, GetDatasets.out.collect()) 
    // Use .collectFile to gather all genotyping results into a single CSV file
    GenotypingFinal_ch = GenotypingResults.out.results.map { tup -> tup[1] }
        .collectFile(
            name: 'final_genotyping_results.csv',
            keepHeader: true,
        )

    // MARKERS PREPARATION
    // Depending on the protocol, either use FluMutDB for Avian or the predefined markers directory for Human
    if (params.protocol == "AVIAN") {
        FluMutDB(SubtypeMerged_ch)
        ch_database = FluMutDB.out
        MarkersFiles(FluMutDB.out) 
    } else {
        def humanMarkersDir = file("${projectDir}/../protocols/HUMAN/v1/markers")
        MarkersFiles(humanMarkersDir)
    }
    ch_markerfiles = MarkersFiles.out

    // MUTATIONS BLOCK (ALL PROTOCOLS)
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
        
    MutationsFinder(Mutations_ch, ch_markerfiles.collect())
    ch_mut = MutationsFinder.out.results.map { _id, mut_files, combined_csv -> [mut_files, combined_csv] }.flatten()
    
    MutationsCompiler_ch = MutationsFinder.out.results
        .map { _sample_id, _mut_files, combined_csv -> combined_csv }
        .collect()
        
    MutationsCompiler(MutationsCompiler_ch)
    ch_raw_mutations = MutationsCompiler.out.results
    ch_mutations_report = MutationsCompiler.out.results

    // APPEND LOGIC INTERCEPTION
    def meta_str = params.metadata ? file(params.metadata).toAbsolutePath().toString() : "" // Resolve the absolute path for metadata if provided

    if (params.get('append')) {
        append_dir_ch = file(params.append, checkIfExists: true)
        
        MergeHistoricalData(SubtypeMerged_ch, GenotypingFinal_ch, ch_raw_mutations, meta_str, append_dir_ch)
        
        final_subtypes_ch   = MergeHistoricalData.out.subtypes
        final_genotyping_ch = MergeHistoricalData.out.genotyping
        final_mutations_ch  = MergeHistoricalData.out.mutations
        final_metadata_ch   = MergeHistoricalData.out.metadata
    } else {
        final_subtypes_ch   = SubtypeMerged_ch
        final_genotyping_ch = GenotypingFinal_ch
        final_mutations_ch  = ch_raw_mutations
        final_metadata_ch   = params.metadata ? channel.fromPath(params.metadata, checkIfExists: true) : channel.of([]) // Create an empty channel if no metadata is provided
    }

    // METADATA MERGE — runs in parallel, only affects published outputs
    if (params.metadata) {
        MetadataMerge(
            final_subtypes_ch,
            final_genotyping_ch,
            final_mutations_ch,
            final_metadata_ch
        )
        published_subtypes_ch   = MetadataMerge.out.subtypes
        published_genotyping_ch = MetadataMerge.out.genotyping
        published_mutations_ch  = MetadataMerge.out.mutations
    } else {
        published_subtypes_ch   = final_subtypes_ch
        published_genotyping_ch = final_genotyping_ch
        published_mutations_ch  = final_mutations_ch
    }

    // AGGREGATED GRAPHIC REPORTS (Using Merged Data)
    CladeGraphicReport(final_genotyping_ch, final_metadata_ch)
    ch_clade_evolution_report = CladeGraphicReport.out.evolution_report

    MutationsGraphicReport(final_mutations_ch, final_metadata_ch)
    MutationsGraphicReport(final_mutations_ch, final_metadata_ch)
    ch_mutations_graphic_report = MutationsGraphicReport.out.report
    
    InteractiveMutationsTable(final_mutations_ch)
    ch_interactive_mutations_table = InteractiveMutationsTable.out.table
    
    if (params.metadata || params.get('append')) {
        DateGraphicReport(final_mutations_ch, final_metadata_ch)
        date_report_ch = DateGraphicReport.out.metadata
    } else {
        date_report_ch = channel.empty()
    }

    // Avaluació i crida de l'informe geogràfic
    def run_geographic_report = check_location_column(params.metadata)
    if (run_geographic_report && params.coordinates) {
        GeographicReport(final_genotyping_ch, final_metadata_ch, file(params.coordinates, checkIfExists: true))
        ch_geo_report = GeographicReport.out.geo_report
    }

    // CONDITIONALLY RUN INDIVIDUAL GRAPHIC REPORTS
    if (params.get('IndividualReports', false).toString().toLowerCase() == 'true') {
        IndividualMutations_Ch = MutationsFinder.out.results.map { sample_id, _mut_files, combined_csv -> tuple(sample_id, combined_csv) }
        IndividualGraphicReport(IndividualMutations_Ch)
        ch_individual_graphic_report = IndividualGraphicReport.out.report
    }

    // ERROR HANDLING & COMPILATION
    BaseErrors_ch = OrganizeBySample.out.errors
        .mix(
            SubtypeDetection.out.errors,
            GenotypingNextclade.out.errors,
            GenotypingResults.out.errors,
            GetCDS.out.errors,
            TranslateToProtein.out.errors,
            MutationsFinder.out.errors
        )
        
    Errors_ch = BaseErrors_ch.groupTuple()

    CompileErrors(Errors_ch)
    // Merge all individual error logs into a single comprehensive log file
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
    subtype = published_subtypes_ch
    datasets = GetDatasets.out
    database = ch_database
    markerfiles = ch_markerfiles
    results = published_genotyping_ch
    CDS = ch_cds
    prot = ch_prot
    graphic_report = CladeGraphicReport.out.report
    clade_evolution_report = ch_clade_evolution_report
    mutations_graphic_report = ch_mutations_graphic_report
    individual_graphic_report = ch_individual_graphic_report
    interactive_mutations_table = ch_interactive_mutations_table
    date_report = date_report_ch
    geo_report = ch_geo_report
    geo_report = ch_geo_report
    mut = ch_mut
    mutations_report = published_mutations_ch
    merged_metadata = final_metadata_ch
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
    merged_metadata {
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
    clade_evolution_report {
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
    geo_report {
        path { "${projectDir}/../${params.outDir}/graphic_reports" }
        mode "copy"
    }
    geo_report {
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
