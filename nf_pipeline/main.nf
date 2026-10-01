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
include { MergeReports              } from './modules/MergeReports.nf'
include { PhyloHistoricalNextclade; PhyloHistoricalNextcladeReuse } from './modules/phylogeny/PhyloHistoricalNextclade'
include {
    PHYLOGENETICS as PHYLOGENETICS_H3N2; PHYLOGENETICS as PHYLOGENETICS_H1N1PDM09;
    parseAaPositions; parseAaPositionsFile; validatePhyloParams; resolvePhyloConfig; phyloSubtypeTag
} from './subworkflows/Phylogenetics'

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

    // PHYLOGENETICS PARAMETER VALIDATION: fail fast, before any process runs
    def run_phylogenetics = params.get('phylogenetics', false).toString().toLowerCase() == 'true'
    def phylo_subtypes = []
    def phylo_positions_h3n2 = []
    def phylo_positions_h1n1 = []
    if (run_phylogenetics) {
        def phylo_errors = validatePhyloParams(params)
        if (phylo_errors) {
            exit 1, "PHYLOGENETICS ERROR:\n  " + phylo_errors.join("\n  ")
        }
        phylo_subtypes = params.phyloSubtype.toString().split(',').collect { s -> s.trim() }.findAll { s -> s }
        phylo_positions_h3n2 = parseAaPositions(params.aaPositions, params.aaGene)
        phylo_positions_h1n1 = parseAaPositions(params.aaPositionsH1n1, params.aaGene)
        // --aa-positions-file replaces --aa-positions/--aa-positions-h1n1 (and the aaPositions default) for the
        // subtypes it actually has rows for; a subtype it never mentions keeps its --aa-positions(-h1n1) value
        if (params.aaPositionsFile) {
            def file_positions = parseAaPositionsFile(params.aaPositionsFile)
            def default_aa_positions = "HA1:62,HA1:145,HA1:239"
            if ((params.aaPositions != null && params.aaPositions.toString() != default_aa_positions) || params.aaPositionsH1n1 != null) {
                log.warn "PHYLOGENETICS: --aa-positions-file is given; --aa-positions/--aa-positions-h1n1 are ignored for the subtypes it covers."
            }
            if (file_positions.containsKey('H3N2'))      phylo_positions_h3n2 = file_positions['H3N2']
            if (file_positions.containsKey('H1N1pdm09')) phylo_positions_h1n1 = file_positions['H1N1pdm09']
        }
        // The flat phyloReference/phyloReferenceFallback/phyloRootStrains/phyloOutgroup params are H3N2-only
        // overrides (see resolvePhyloConfig()): warn if they were changed from their built-in default but H3N2 isn't requested
        if (!('H3N2' in phylo_subtypes)) {
            def h3n2_defaults = resolvePhyloConfig('H3N2', [:])
            if (params.phyloReference != h3n2_defaults.reference || params.phyloReferenceFallback != h3n2_defaults.referenceFallback ||
                params.phyloRootStrains != h3n2_defaults.rootStrains || params.phyloOutgroup != h3n2_defaults.outgroup) {
                log.warn "PHYLOGENETICS: --phyloReference/--phyloReferenceFallback/--phyloRootStrains/--phyloOutgroup are H3N2 overrides, but H3N2 is not in --phyloSubtype (${params.phyloSubtype}); they have no effect."
            }
        }
    }

    // INPUT & INITIAL FOLDER ORGANIZATION
    SampleInput_ch = channel
        .fromPath(params.inputFasta, checkIfExists: true)
        .splitFasta(record: [id: true])
        .map { rec -> rec.id.tokenize('[|_]')[0] } // Extract sample ID from FASTA header using the first token before '|' or '_'
        .unique()
        .filter { sample_id ->
            // Sample IDs are embedded in shell/Python scripts and file paths: skip IDs with characters that could break them
            def safe = sample_id && !(sample_id =~ /[\s'"`$\\{};&<>*?\/(),]/) && !(sample_id in ['.', '..'])
            if (!safe) log.warn "Skipping sample '${sample_id}': sample IDs cannot be empty or contain whitespace, quotes or any of \$ \\ { } ; & < > * ? / ( ) ,"
            return safe
        }

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
        )

    // DATASET PREPARATION
    // --append + --phylogenetics: a subtype's samples may exist ONLY in the append directory's history (e.g. an
    // H1N1pdm09 tree resumed from a new FASTA that only has H3N2 samples). GetDatasets only scans the CURRENT
    // run's inferred subtypes, so without this its Nextclade dataset would never be fetched and that subtype's
    // PhyloSubclades/PhyloHistoricalNextclade/PhyloReport would silently never run (no failed task to notice).
    // Concatenating the historical inferred_subtypes.csv in here (content only; the real old+new merge for
    // publishing still happens in MergeHistoricalData below) is enough, since GetDatasets just scans for H tags.
    GetDatasetsInput_ch = SubtypeMerged_ch
    if (run_phylogenetics && params.get('append')) {
        def hist_subtypes_file = file("${params.append}/inferred_subtypes.csv")
        if (hist_subtypes_file.exists()) {
            GetDatasetsInput_ch = SubtypeMerged_ch
                .mix(channel.fromPath(hist_subtypes_file))
                .collectFile(name: 'inferred_subtypes_for_datasets.csv', keepHeader: true, sort: { f -> f.name })
        } else {
            def msg = "PHYLOGENETICS: --append directory '${params.append}' has no inferred_subtypes.csv; historical subtypes are not considered when fetching Nextclade datasets."
            println "WARN: ${msg}"
        }
    }
    GetDatasets(GetDatasetsInput_ch)

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
        def humanMarkersDir = file("${projectDir}/../protocols/HUMAN/v2/markers")
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

    // PHYLOGENETICS (OPTIONAL, --phylogenetics): one annotated whole-genome tree per requested --phyloSubtype
    ch_phylogeny_h3n2 = channel.empty()
    ch_phylogeny_h1n1 = channel.empty()
    ch_phylo_errors = channel.empty()
    ch_phylo_report = channel.empty()
    ch_phylo_tree_problems = channel.empty()
    if (run_phylogenetics) {
        // --append: MergeHistoricalData only merges the summary tables, so PHYLOGENETICS would otherwise only see
        // this run's samples. Build phylo-only genotyping/sample-dir/Nextclade channels covering both runs.
        if (params.get('append')) {
            def append_path = params.append.toString()
            def hist_samples_dir = file("${append_path}/samples")
            if (!hist_samples_dir.exists() || !hist_samples_dir.isDirectory()) {
                def msg = "PHYLOGENETICS: --append directory '${append_path}' has no 'samples' folder; only this run's own samples are considered for the tree(s)."
                println "WARN: ${msg}"
            }

            // (a) genotyping info parsed from the merged (old+new) subtypes table
            PhyloGenotypingInfo_ch = final_subtypes_ch
                .splitCsv(header: true)
                .map { row ->
                    def subtype = row['inferred_subtype'] ?: ''
                    def h_tag = subtype.find(/H\d+/) ?: "Hx"
                    def n_tag = subtype.find(/N\d+/) ?: "Nx"
                    tuple(row['Sample_ID'], h_tag, n_tag, row['pathotype'])
                }

            // (b) sample dirs: new OrganizeBySample results, plus historical <appendDir>/samples/<id> dirs that
            // have a segments/ folder and are not superseded by this run (new run wins on a duplicate ID)
            NewSampleIds_ch = OrganizeBySample.out.results.map { sample_id, _dir -> sample_id }.toList()
            AllHistSampleDirs_ch = channel.fromPath("${append_path}/samples/*", type: 'dir')
            AllHistSampleDirs_ch
                .filter { d -> !file("${d}/segments").isDirectory() }
                .map { d -> d.name }
                .toList()
                .subscribe { ids ->
                    if (ids) {
                        def msg = "PHYLOGENETICS: ${ids.size()} historical sample folder(s) under '${append_path}/samples' have no segments/ subfolder and are skipped: ${ids.take(10).join(', ')}${ids.size() > 10 ? ', ...' : ''}."
                        println "WARN: ${msg}"
                    }
                }
            // NewSampleIds_ch.toList() emits a single List item; combine() would otherwise treat that List as
            // several fields to spread rather than one opaque value (Map is never spread this way)
            NewSampleIdsMap_ch = NewSampleIds_ch.map { ids -> [ids: ids] }
            HistSampleDirs_ch = AllHistSampleDirs_ch
                .filter { d -> file("${d}/segments").isDirectory() }
                .map { d -> tuple(d.name, d) }
                .combine(NewSampleIdsMap_ch)
                .filter { sample_id, _dir, new_ids_map -> !(sample_id in new_ids_map.ids) }
                .map { sample_id, dir, _new_ids_map -> tuple(sample_id, dir) }
            PhyloSampleDirs_ch = OrganizeBySample.out.results.mix(HistSampleDirs_ch)

            // (c) Nextclade results: reuse a historical sample's persisted CSV (samples/<id>/nextclade_results.csv,
            // written by GenotypingNextclade) when present; otherwise rerun Nextclade on its HA segment against the
            // current subtype's dataset, so it is never left with subclade NA just because it predates that file.
            HistNextcladeReuse_ch = HistSampleDirs_ch
                .map { sample_id, dir -> tuple(sample_id, dir, file("${dir}/nextclade_results.csv")) }
                .filter { _sample_id, _dir, csv -> csv.exists() }
                .map { sample_id, _dir, csv -> tuple(sample_id, csv) }
            HistNextcladeRerunInput_ch = HistSampleDirs_ch
                .map { sample_id, dir -> tuple(sample_id, dir, file("${dir}/nextclade_results.csv")) }
                .filter { _sample_id, _dir, csv -> !csv.exists() }
                .map { sample_id, dir, _csv -> tuple(sample_id, file("${dir}/segments/${sample_id}_HA.fasta")) }
                .filter { _sample_id, ha_fasta -> ha_fasta.exists() }
                .combine(PhyloGenotypingInfo_ch.map { sample_id, h_tag, _n_tag, _pathotype -> tuple(sample_id, h_tag) }, by: 0)
                .combine(GetDatasets.out.flatMap { dataset_dirs -> dataset_dirs })
                .filter { _sample_id, _ha_fasta, h_tag, dataset_dir -> dataset_dir.name.contains(h_tag) }
                .map { sample_id, ha_fasta, _h_tag, dataset_dir -> tuple(sample_id, ha_fasta, dataset_dir) }
            PhyloHistoricalNextclade(HistNextcladeRerunInput_ch)
            // Reused historical CSVs are on disk as plain "nextclade_results.csv" (no per-sample suffix): rename to
            // nextclade_results_<sample_id>.csv, the name PhyloSubclades' own glob requires to recover the sample_id
            // (see PhyloHistoricalNextcladeReuse's comment).
            PhyloHistoricalNextcladeReuse(HistNextcladeReuse_ch)
            ch_phylo_hist_errors = PhyloHistoricalNextclade.out.errors
            PhyloNextcladeResults_ch = GenotypingNextclade.out.results.mix(PhyloHistoricalNextcladeReuse.out.results, PhyloHistoricalNextclade.out.results)
        } else {
            PhyloGenotypingInfo_ch = GenotypingInfo_ch
            PhyloSampleDirs_ch = OrganizeBySample.out.results
            PhyloNextcladeResults_ch = GenotypingNextclade.out.results
            ch_phylo_hist_errors = channel.empty()
        }

        if ('H3N2' in phylo_subtypes) {
            def cfg_h3n2 = resolvePhyloConfig('H3N2', params)
            PHYLOGENETICS_H3N2(
                PhyloGenotypingInfo_ch,
                PhyloSampleDirs_ch,
                PhyloNextcladeResults_ch,
                GetDatasets.out,
                final_metadata_ch,
                phylo_positions_h3n2.collect { pos -> pos.label },
                'H3N2',
                phyloSubtypeTag('H3N2'),
                cfg_h3n2.reference,
                cfg_h3n2.referenceFallback,
                cfg_h3n2.rootStrains,
                cfg_h3n2.outgroup
            )
            ch_phylogeny_h3n2 = PHYLOGENETICS_H3N2.out.outputs
            ch_phylo_errors = ch_phylo_errors.mix(PHYLOGENETICS_H3N2.out.errors)
            ch_phylo_report = ch_phylo_report.mix(PHYLOGENETICS_H3N2.out.report)
            ch_phylo_tree_problems = ch_phylo_tree_problems.mix(PHYLOGENETICS_H3N2.out.tree_problems)
        }
        if ('H1N1pdm09' in phylo_subtypes) {
            def cfg_h1n1 = resolvePhyloConfig('H1N1pdm09', params)
            PHYLOGENETICS_H1N1PDM09(
                PhyloGenotypingInfo_ch,
                PhyloSampleDirs_ch,
                PhyloNextcladeResults_ch,
                GetDatasets.out,
                final_metadata_ch,
                phylo_positions_h1n1.collect { pos -> pos.label },
                'H1N1pdm09',
                phyloSubtypeTag('H1N1pdm09'),
                cfg_h1n1.reference,
                cfg_h1n1.referenceFallback,
                cfg_h1n1.rootStrains,
                cfg_h1n1.outgroup
            )
            ch_phylogeny_h1n1 = PHYLOGENETICS_H1N1PDM09.out.outputs
            ch_phylo_errors = ch_phylo_errors.mix(PHYLOGENETICS_H1N1PDM09.out.errors)
            ch_phylo_report = ch_phylo_report.mix(PHYLOGENETICS_H1N1PDM09.out.report)
            ch_phylo_tree_problems = ch_phylo_tree_problems.mix(PHYLOGENETICS_H1N1PDM09.out.tree_problems)
        }
        ch_phylo_errors = ch_phylo_errors.mix(ch_phylo_hist_errors)
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
            MutationsFinder.out.errors,
            ch_phylo_errors
        )
        
    Errors_ch = BaseErrors_ch.groupTuple()

    CompileErrors(Errors_ch)
    // Merge all individual error logs into a single comprehensive log file. Also folds in PHYLOGENETICS'
    // tree_problems (trees an 'ignore'd task silently dropped, per subtype: see subworkflows/Phylogenetics.nf),
    // since sample-keyed CompileErrors has no natural place for a tree-level (not sample-level) problem.
    ErrorsMerged_ch = CompileErrors.out
        .map { sample_id, log_file ->
            def content = log_file.text
            return "========================================\n" +
                   "Errors for Sample: ${sample_id}\n" +
                   "========================================\n" +
                   "${content}\n"
        }
        .mix(
            ch_phylo_tree_problems.map { log_file ->
                "========================================\n" +
                "Phylogenetics: trees requested but never produced a report\n" +
                "========================================\n" +
                "${log_file.text}\n"
            }
        )
        .collectFile(
            name: 'pipeline_errors.log',
        )

    all_reports_ch = CladeGraphicReport.out.report
        .mix(ch_clade_evolution_report)
        .mix(ch_mutations_graphic_report)
        .mix(ch_interactive_mutations_table)
        .mix(date_report_ch)
        .mix(ch_geo_report)
        .mix(ch_phylo_report)
        .collect()

    MergeReports(all_reports_ch)
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
    individual_graphic_report = ch_individual_graphic_report
    mut = ch_mut
    mutations_report = published_mutations_ch
    index = MergeReports.out.index
    merged_metadata = final_metadata_ch
    errors = CompileErrors.out.map { _id, log -> log }
    errors_merged = ErrorsMerged_ch
    phylogeny_h3n2 = ch_phylogeny_h3n2
    phylogeny_h1n1 = ch_phylogeny_h1n1
    nextclade_csv = GenotypingNextclade.out.sample_csv

    onComplete:
    // Processes use errorStrategy 'ignore', so a failed task silently drops its outputs: make that visible
    if (workflow.stats.ignoredCount > 0) {
        log.warn "${workflow.stats.ignoredCount} task(s) failed and were ignored, so some samples or reports may be missing. Check pipeline_errors.log and .nextflow.log for details."
    }
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
    individual_graphic_report {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    index {
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
    nextclade_csv {
        path { "${projectDir}/../${params.outDir}" }
        mode "copy"
    }
    phylogeny_h3n2 {
        path { "${projectDir}/../${params.outDir}/phylogeny/H3N2" }
        mode "copy"
    }
    phylogeny_h1n1 {
        path { "${projectDir}/../${params.outDir}/phylogeny/H1N1pdm09" }
        mode "copy"
    }
}
