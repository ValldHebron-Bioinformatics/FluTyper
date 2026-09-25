#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Optional phylogenetics analysis (--phylogenetics): annotated whole-genome tree of the phyloSubtype samples,
// following Koumaty et al. 2026 (Microbial Genomics 12:001810). Each step is one process in modules/phylogeny/.

include { PhyloPreflight  } from '../modules/phylogeny/PhyloPreflight'
include { PhyloCoverageQC } from '../modules/phylogeny/PhyloCoverageQC'
include { PhyloAlignment  } from '../modules/phylogeny/PhyloAlignment'
include { PhyloTree       } from '../modules/phylogeny/PhyloTree'
include { PhyloSubclades  } from '../modules/phylogeny/PhyloSubclades'
include { PhyloPostTree   } from '../modules/phylogeny/PhyloPostTree'
include { PhyloReport     } from '../modules/phylogeny/PhyloReport'

workflow PHYLOGENETICS {
    take:
    genotyping_info     // [sample_id, h_tag, n_tag, pathotype] from SubtypeDetection
    sample_dirs         // [sample_id, sample_dir] from OrganizeBySample
    nextclade_results   // [sample_id, nextclade csv] from GenotypingNextclade
    datasets            // Nextclade dataset folders from GetDatasets
    metadata            // metadata CSV, or [] when --metadata is not given
    aa_positions        // resolved --aa-positions as a list of GENE:POS labels

    main:
    def resources = params.protocols[params.protocol].resources
    def ref_genomes = file("${resources}/${params.protocol}_references.fasta", checkIfExists: true)
    def cds_refs = file("${resources}/CDS_references.fasta", checkIfExists: true)
    def panel_spec = file(params.phyloPanels, checkIfExists: true)
    def lab_map = params.labMap ? file(params.labMap, checkIfExists: true) : []

    // Startup checks: tools present, reference codon map built, --aa-positions validated. Every later step waits for it.
    PhyloPreflight(
        groovy.json.JsonOutput.toJson(params.phyloTools).bytes.encodeBase64().toString(),
        ref_genomes,
        cds_refs,
        params.phyloSubtype,
        params.phyloReference,
        params.phyloReferenceFallback,
        params.phyloRootStrains,
        params.phyloSegmentOrder,
        aa_positions.join(',')
    )

    // Candidate samples: those typed as phyloSubtype
    candidates = genotyping_info
        .filter { _sample_id, h_tag, n_tag, _pathotype -> "${h_tag}${n_tag}" == params.phyloSubtype }
        .map { sample_id, _h_tag, _n_tag, _pathotype -> tuple(sample_id) }
        .join(sample_dirs)
    candidates.count().subscribe { n ->
        if (n == 0) log.warn "PHYLOGENETICS: no ${params.phyloSubtype} samples were found, so no phylogeny is built."
    }

    // Coverage QC gates entry into the whole-genome phylogeny
    PhyloCoverageQC(candidates, PhyloPreflight.out.coordinate_refs, params.phyloMinCoverage, params.phyloSegmentOrder)
    qc_table = PhyloCoverageQC.out.results
        .map { _sample_id, _sample_dir, qc_tsv -> qc_tsv }
        .collectFile(name: 'phylo_coverage_qc.tsv', keepHeader: true, sort: { f -> f.name })
    passing = PhyloCoverageQC.out.results
        .filter { _sample_id, _sample_dir, qc_tsv ->
            qc_tsv.readLines().drop(1).every { line -> line.split('\t')[5] == 'PASS' }
        }
        .map { _sample_id, sample_dir, _qc_tsv -> sample_dir }
        .toSortedList { a, b -> a.name <=> b.name }

    PhyloAlignment(
        passing,
        PhyloPreflight.out.references,
        PhyloPreflight.out.coordinate_refs,
        PhyloPreflight.out.coordinate_table,
        qc_table,
        params.phyloSegmentOrder
    )
    PhyloTree(PhyloAlignment.out.alignment, phyloRefId(params.phyloOutgroup), params.model, params.bootstrap, params.seed)

    // Subclades: per-sample Nextclade results, plus the root strains run against the same dataset
    nextclade_join = candidates
        .map { sample_id, _sample_dir -> tuple(sample_id, true) }
        .join(nextclade_results, remainder: true)
        .filter { _sample_id, is_candidate, _csv -> is_candidate }
    nextclade_join
        .filter { _sample_id, _is_candidate, csv -> csv == null }
        .subscribe { sample_id, _is_candidate, _csv ->
            log.warn "PHYLOGENETICS: no Nextclade result for sample ${sample_id}; its subclade is NA and it is not clustered."
        }
    nextclade_csvs = nextclade_join
        .filter { _sample_id, _is_candidate, csv -> csv != null }
        .map { _sample_id, _is_candidate, csv -> csv }
        .toSortedList { a, b -> a.name <=> b.name }
    def h_tag = params.phyloSubtype.find(/H\d+/)
    dataset = datasets.flatten()
        .filter { dataset_dir -> dataset_dir.name == "nextclade_${h_tag}_dataset" }
        .toList()
        .map { found ->
            if (!found) log.warn "PHYLOGENETICS: no nextclade_${h_tag}_dataset was retrieved by GetDatasets, so subclades, clusters, annotations, the figure and the manifest are not produced."
            found
        }
        .filter { found -> found }
        .map { found -> found[0] }
    PhyloSubclades(nextclade_csvs, PhyloPreflight.out.references, dataset)

    PhyloPostTree(
        PhyloTree.out.tree,
        PhyloAlignment.out.alignment,
        PhyloPreflight.out.gene_map,
        PhyloAlignment.out.segments,
        PhyloSubclades.out.subclades,
        params.clusterThreshold
    )

    // Resolved parameters and input checksums for the run manifest
    def param_names = [
        'protocol', 'inputFasta', 'metadata', 'outDir', 'phyloSubtype', 'phyloReference', 'phyloReferenceFallback',
        'phyloRootStrains', 'phyloOutgroup', 'phyloSegmentOrder', 'phyloMinCoverage', 'model', 'bootstrap', 'seed',
        'clusterThreshold', 'aaPositions', 'aaGene', 'minPanelCoverage', 'labMap', 'periodBin', 'width', 'height',
        'dpi', 'phyloPanels'
    ]
    def run_params = param_names.collectEntries { name -> [(name): params[name]?.toString()] }
    run_params.aaPositionsResolved = aa_positions
    def run_info = [pipeline_version: params.version, started: workflow.start.toString(), parameters: run_params]
    def checksum_inputs = [file(params.inputFasta), ref_genomes, cds_refs, panel_spec]
    if (params.metadata) checksum_inputs << file(params.metadata)
    if (params.labMap) checksum_inputs << file(params.labMap)

    PhyloReport(
        PhyloTree.out.tree,
        PhyloTree.out.support,
        PhyloSubclades.out.subclades,
        PhyloPostTree.out.clusters,
        PhyloPostTree.out.proteins,
        PhyloPreflight.out.positions,
        metadata,
        lab_map,
        panel_spec,
        PhyloPreflight.out.versions,
        qc_table,
        PhyloPreflight.out.coordinate_table,
        PhyloPreflight.out.root_table,
        PhyloAlignment.out.segments,
        checksum_inputs,
        params.periodBin,
        params.minPanelCoverage,
        params.width,
        params.height,
        params.dpi,
        groovy.json.JsonOutput.toJson(run_info).bytes.encodeBase64().toString()
    )

    outputs = qc_table.mix(
        PhyloPreflight.out.versions,
        PhyloPreflight.out.gene_map,
        PhyloPreflight.out.positions,
        PhyloPreflight.out.coordinate_table,
        PhyloPreflight.out.root_table,
        PhyloAlignment.out.alignment,
        PhyloAlignment.out.segments,
        PhyloTree.out.raw_tree,
        PhyloTree.out.tree,
        PhyloTree.out.support,
        PhyloTree.out.report,
        PhyloTree.out.log,
        PhyloSubclades.out.subclades,
        PhyloPostTree.out.nucleotides,
        PhyloPostTree.out.annotated_tree,
        PhyloPostTree.out.proteins,
        PhyloPostTree.out.clusters,
        PhyloPostTree.out.mutations,
        PhyloReport.out.annotations,
        PhyloReport.out.sources,
        PhyloReport.out.aa_states,
        PhyloReport.out.svg,
        PhyloReport.out.png,
        PhyloReport.out.panels,
        PhyloReport.out.manifest,
        PhyloReport.out.html
    )

    emit:
    outputs = outputs                     // every file published to <outDir>/phylogeny
    errors  = PhyloCoverageQC.out.errors  // per-sample QC errors, merged into pipeline_errors.log
    report  = PhyloReport.out.html        // interactive tree, added to the index.html dashboard
}

/*
 * Parse --aa-positions into an ordered list of [gene, pos, label] maps.
 * Accepts a comma-separated string (e.g. "HA1:62,HA1:145,239") or a list of such strings
 * (e.g. from a -params-file YAML list). Bare integers fall back to defaultGene.
 * null, empty or a value-less flag (true) means no amino acid columns.
 * Syntax errors throw IllegalArgumentException; range checks need the reference translation
 * and are done by PhyloPreflight.
 */
def parseAaPositions(value, defaultGene) {
    if (value == null || value == true || value.toString().trim() == '') return []

    def gene_pattern = /^[A-Za-z0-9][A-Za-z0-9-]*$/
    def default_gene = defaultGene == null ? '' : defaultGene.toString().trim().toUpperCase()

    def raw_items = (value instanceof Collection ? value : [value]).collect { v -> v.toString() }
    def tokens = []
    raw_items.each { item -> tokens.addAll(item.split(',', -1).collect { t -> t.trim() }) }

    def positions = []
    def seen = [] as Set
    tokens.each { token ->
        if (!token) {
            throw new IllegalArgumentException("--aa-positions: empty entry in '${raw_items.join(',')}'. Use a comma-separated list such as HA1:62,HA1:145,239.")
        }
        def gene
        def pos_str
        if (token ==~ /^\d+$/) {
            if (!(default_gene ==~ gene_pattern)) {
                throw new IllegalArgumentException("--aa-gene: invalid gene name '${defaultGene}' used for bare position '${token}'.")
            }
            gene = default_gene
            pos_str = token
        } else {
            def m = token =~ /^([A-Za-z0-9][A-Za-z0-9-]*):(\d+)$/
            if (!m.matches()) {
                throw new IllegalArgumentException("--aa-positions: malformed entry '${token}'. Expected GENE:POS (e.g. HA1:62) or a bare integer.")
            }
            gene = m.group(1).toUpperCase()
            pos_str = m.group(2)
        }
        def pos = pos_str.toBigInteger()
        if (pos < 1 || pos > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("--aa-positions: position '${token}' must be between 1 and the length of the gene.")
        }
        def label = "${gene}:${pos}".toString()
        if (seen.contains(label)) {
            throw new IllegalArgumentException("--aa-positions: position ${label} is given more than once.")
        }
        seen << label
        positions << [gene: gene, pos: pos.intValue(), label: label]
    }
    return positions
}

/*
 * Validate the phylogenetics parameters. Returns a list of error messages (empty when valid).
 */
def validatePhyloParams(Map p) {
    def errors = []

    if (p.protocol != 'HUMAN') {
        errors << "--phylogenetics is only available for the HUMAN protocol (got '${p.protocol}')."
    }
    if (p.phyloSubtype != 'H3N2') {
        errors << "--phyloSubtype '${p.phyloSubtype}' is not supported: only H3N2 (reference A/Darwin/6/2021) is available."
    }
    [ 'phyloMinCoverage', 'minPanelCoverage' ].each { k ->
        if (!isPhyloNumber(p[k]) || p[k].toString().toBigDecimal() < 0 || p[k].toString().toBigDecimal() > 1) {
            errors << "--${k} must be a number between 0 and 1 (got '${p[k]}')."
        }
    }
    if (!isPhyloNumber(p.clusterThreshold) || p.clusterThreshold.toString().toBigDecimal() <= 0) {
        errors << "--clusterThreshold must be a positive number of substitutions per site (got '${p.clusterThreshold}')."
    }
    if (!(p.bootstrap.toString() ==~ /^\d+$/) || (p.bootstrap.toString().toBigInteger() != 0 && p.bootstrap.toString().toBigInteger() < 1000)) {
        errors << "--bootstrap must be 0 (no support values) or at least 1000 ultrafast bootstrap replicates (got '${p.bootstrap}')."
    }
    if (!(p.seed.toString() ==~ /^\d+$/)) {
        errors << "--seed must be a non-negative integer (got '${p.seed}')."
    }
    if (!p.model || p.model == true || !(p.model.toString() ==~ /^[A-Za-z0-9+{}.,\/-]+$/)) {
        errors << "--model must be an IQ-TREE model string such as GTR+G4 (got '${p.model}')."
    }
    if (!(p.periodBin in ['month', 'week', 'none'])) {
        errors << "--periodBin must be one of month, week or none (got '${p.periodBin}')."
    }
    [ 'width', 'height', 'dpi' ].each { k ->
        if (!isPhyloNumber(p[k]) || p[k].toString().toBigDecimal() <= 0) {
            errors << "--${k} must be a positive number (got '${p[k]}')."
        }
    }
    def order = p.phyloSegmentOrder ? p.phyloSegmentOrder.toString().split(',').collect { s -> s.trim() } : []
    if (order.sort(false) != (p.segments ?: []).sort(false)) {
        errors << "--phyloSegmentOrder must list each of ${p.segments} exactly once (got '${p.phyloSegmentOrder}')."
    }
    def roots = p.phyloRootStrains ? p.phyloRootStrains.toString().split(',').collect { s -> s.trim() }.findAll { s -> s } : []
    if (!(p.phyloOutgroup in roots)) {
        errors << "--phyloOutgroup '${p.phyloOutgroup}' must be one of --phyloRootStrains."
    }
    if (p.labMap && !file(p.labMap.toString()).exists()) {
        errors << "--labMap file not found: ${p.labMap}"
    }
    try {
        parseAaPositions(p.aaPositions, p.aaGene)
    } catch (IllegalArgumentException e) {
        errors << e.message
    }
    return errors
}

// True when v is a plain number (not null, not a value-less flag)
def isPhyloNumber(v) {
    return v != null && v != true && v.toString().isBigDecimal()
}

// Tip identifier used in the alignment and tree for a reference strain (e.g. A/Darwin/6/2021 -> REF_A_Darwin_6_2021)
def phyloRefId(String strain) {
    return 'REF_' + strain.replaceAll(/[^A-Za-z0-9.-]/, '_')
}
