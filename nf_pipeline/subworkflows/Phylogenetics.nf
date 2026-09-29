#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Optional phylogenetics analysis (--phylogenetics): annotated whole-genome tree of the phyloSubtype samples,
// following Koumaty et al. 2026 (Microbial Genomics 12:001810). Each step is one process in modules/phylogeny/.

include { PhyloPreflight        } from '../modules/phylogeny/PhyloPreflight'
include { PhyloCoverageQC       } from '../modules/phylogeny/PhyloCoverageQC'
include { PhyloAlignment        } from '../modules/phylogeny/PhyloAlignment'
include { PhyloSegmentAlignment } from '../modules/phylogeny/PhyloSegmentAlignment'
include { PhyloTree             } from '../modules/phylogeny/PhyloTree'
include { PhyloSubclades        } from '../modules/phylogeny/PhyloSubclades'
include { PhyloPostTree         } from '../modules/phylogeny/PhyloPostTree'
include { PhyloReport           } from '../modules/phylogeny/PhyloReport'

workflow PHYLOGENETICS {
    take:
    genotyping_info      // [sample_id, h_tag, n_tag, pathotype] from SubtypeDetection
    sample_dirs          // [sample_id, sample_dir] from OrganizeBySample
    nextclade_results    // [sample_id, nextclade csv] from GenotypingNextclade
    datasets             // Nextclade dataset folders from GetDatasets
    metadata             // metadata CSV, or [] when --metadata is not given
    aa_positions         // resolved --aa-positions/--aa-positions-h1n1 as a list of GENE:POS labels
    subtype              // e.g. "H3N2" or "H1N1pdm09"
    subtype_tag          // GenotypingInfo_ch tag for this subtype, e.g. "H3N2" or "H1N1" (see phyloSubtypeTag())
    reference             // resolved alignment coordinate reference for this subtype
    reference_fallback    // resolved coordinate reference for segments the alignment reference lacks
    root_strains          // resolved comma-separated root strains for this subtype
    outgroup               // resolved root strain used to root this subtype's tree

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
        subtype,
        reference,
        reference_fallback,
        root_strains,
        params.phyloSegmentOrder,
        aa_positions.join(',')
    )

    // Candidate samples: those typed as this subtype. A sample_id here with no matching sample_dirs entry (e.g. an
    // --append historical sample dropped for lacking segments/, or an appendDir folder that never existed) would
    // otherwise be silently dropped by a plain inner join(); warn instead of losing it quietly.
    def subtype_matched_ch = genotyping_info
        .filter { _sample_id, h_tag, n_tag, _pathotype -> "${h_tag}${n_tag}" == subtype_tag }
        .map { sample_id, _h_tag, _n_tag, _pathotype -> tuple(sample_id) }
    // Plain inner join, exactly as before this round's changes: this is what actually drives the tree(s), so its
    // behaviour for a normal (non-append) run must stay identical to before.
    candidates = subtype_matched_ch.join(sample_dirs)
    candidates.count().subscribe { n ->
        if (n == 0) {
            def msg = "PHYLOGENETICS(${subtype}): no ${subtype_tag} samples were found, so no phylogeny is built."
            log.warn msg
            println "WARN: ${msg}"
        }
    }

    // Side-channel warning only (does not feed `candidates`): a plain inner join() above silently drops any
    // subtype-matched sample_id with no matching sample_dirs entry (e.g. an --append historical sample dropped for
    // lacking segments/, or an appendDir folder that never existed). Using join(remainder: true) here instead would
    // behave like a FULL outer join and incorrectly reintroduce samples of OTHER subtypes that merely have a dir, so
    // the missing set is computed separately via a list difference instead. Built straight from genotyping_info
    // (not from subtype_matched_ch's tuple(sample_id) items): a single-element tuple's arity is ambiguous under a
    // 1-arg .map{}, so re-deriving the plain sample_id string here with its real 4-arg destructuring avoids that.
    def matched_ids_ch = genotyping_info
        .filter { _sample_id, h_tag, n_tag, _pathotype -> "${h_tag}${n_tag}" == subtype_tag }
        .map { sample_id, _h_tag, _n_tag, _pathotype -> sample_id }
        .toList()
        .map { list -> [ids: list] }
    def dir_ids_ch = sample_dirs.map { sample_id, _dir -> sample_id }.toList()
        .map { list -> [ids: list] }
    matched_ids_ch.combine(dir_ids_ch)
        .map { matched_m, dir_m -> matched_m.ids - dir_m.ids }
        .filter { missing -> missing }
        .subscribe { missing ->
            def msg = "PHYLOGENETICS(${subtype}): ${missing.size()} sample(s) typed as ${subtype_tag} have no sample folder (dropped from --append history, or missing segments/) and are excluded from every tree: ${missing.take(10).join(', ')}${missing.size() > 10 ? ', ...' : ''}."
            log.warn msg
            println "WARN: ${msg}"
        }

    // Coverage QC gates entry into every tree: the whole-genome tree needs every segment to PASS, a segment tree
    // only needs its own segment to PASS (so its candidate pool is usually larger than the whole-genome one)
    PhyloCoverageQC(candidates, PhyloPreflight.out.coordinate_refs, params.phyloMinCoverage, params.phyloSegmentOrder)
    qc_table = PhyloCoverageQC.out.results
        .map { _sample_id, _sample_dir, qc_tsv -> qc_tsv }
        .collectFile(name: 'phylo_coverage_qc.tsv', keepHeader: true, sort: { f -> f.name })

    // --phyloWholeGenome / --phyloSegmentTrees: which tree(s) to build (validatePhyloParams already checked at
    // least one is requested and that --phyloSegmentTrees only names real segments)
    def whole_genome = params.phyloWholeGenome != null && params.phyloWholeGenome != false && params.phyloWholeGenome.toString().toLowerCase() != 'false'
    def all_segments = params.phyloSegmentOrder.toString().split(',').collect { s -> s.trim() }.findAll { s -> s }
    def requested_trees = (params.phyloSegmentTrees ?: 'none').toString().trim()
    def segment_trees = requested_trees.equalsIgnoreCase('none') ? [] :
                         requested_trees.equalsIgnoreCase('all') ? all_segments :
                         requested_trees.split(',').collect { s -> s.trim() }.findAll { s -> s }

    // Whole-genome tree (optional, --phyloWholeGenome): reproduces Koumaty et al. 2026
    def whole_genome_trees_ch = channel.empty()
    if (whole_genome) {
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
        whole_genome_trees_ch = PhyloAlignment.out.alignment
            .combine(PhyloAlignment.out.segments)
            .map { alignment, segments_tsv -> tuple(subtype, "whole_genome", alignment, phyloRefId(outgroup), segments_tsv) }
    }

    // Segment trees (--phyloSegmentTrees): one per requested segment, gated on >= 3 samples passing that segment's
    // own coverage QC. A root strain missing the segment is dropped by PhyloSegmentAlignment, which also picks a
    // fallback outgroup when the configured one lacks the segment (e.g. H3N2 PA: A/Darwin/6/2021 has none).
    def segment_candidates_ch = segment_trees.collect { seg ->
        PhyloCoverageQC.out.results
            .filter { _sample_id, _sample_dir, qc_tsv ->
                def row = qc_tsv.readLines().drop(1).find { line -> line.split('\t')[1] == seg }
                row != null && row.split('\t')[5] == 'PASS'
            }
            .map { _sample_id, sample_dir, _qc_tsv -> sample_dir }
            .toSortedList { a, b -> a.name <=> b.name }
            .map { dirs -> tuple(seg, dirs) }
    }
    def segment_passing_ch = segment_candidates_ch ? segment_candidates_ch.inject(channel.empty()) { acc, ch -> acc.mix(ch) } : channel.empty()
    segment_passing_ch
        .filter { _seg, dirs -> dirs.size() < 3 }
        .subscribe { seg, dirs ->
            def reason = dirs.size() == 0 ? "0 samples passed coverage QC for ${seg}" : "only ${dirs.size()} sample(s) passed coverage QC for ${seg}"
            // log.warn alone lands only in .nextflow.log (nf-test's workflow.stdout/stderr do not capture it), so
            // also print it: this is exactly the kind of silently-dropped-tree visibility this warning exists for.
            def msg = "PHYLOGENETICS(${subtype}): segment tree ${seg} skipped: ${reason} (need at least 3)."
            log.warn msg
            println "WARN: ${msg}"
        }
    def segment_gated_ch = segment_passing_ch.filter { _seg, dirs -> dirs.size() >= 3 }

    PhyloSegmentAlignment(
        segment_gated_ch,
        PhyloPreflight.out.references,
        PhyloPreflight.out.coordinate_refs,
        PhyloPreflight.out.coordinate_table,
        outgroup,
        root_strains
    )
    def segment_outgroup_ch = PhyloSegmentAlignment.out.outgroup.map { seg, f -> tuple(seg, f.text.trim()) }
    def segment_trees_ch = PhyloSegmentAlignment.out.alignment
        .join(segment_outgroup_ch)
        .join(PhyloSegmentAlignment.out.segments)
        .map { seg, alignment, outgroup_tip, segments_tsv -> tuple("${subtype}_${seg}".toString(), seg, alignment, outgroup_tip, segments_tsv) }

    // Every tree of this subtype, tagged (tree_id, tree_folder, alignment, outgroup_tip, segments_tsv):
    // tree_id is globally unique (used for filenames), tree_folder nests this tree's files under phylogeny/<subtype>/
    def all_trees_ch = whole_genome_trees_ch.mix(segment_trees_ch)
    all_trees_ch.count().subscribe { n ->
        if (n == 0) log.warn "PHYLOGENETICS(${subtype}): no tree was built (check --phyloWholeGenome/--phyloSegmentTrees and whether enough samples passed coverage QC)."
    }

    PhyloTree(
        all_trees_ch.map { tree_id, tree_folder, alignment, outgroup_tip, segments_tsv -> tuple(tree_id, tree_folder, alignment, outgroup_tip, segments_tsv) },
        params.model, params.bootstrap, params.seed
    )

    // Subclades: per-sample Nextclade results, plus the root strains run against the same dataset. Shared by every
    // tree of this subtype: subclade assignment (from HA) does not depend on which tree is being built.
    nextclade_join = candidates
        .map { sample_id, _sample_dir -> tuple(sample_id, true) }
        .join(nextclade_results, remainder: true)
        .filter { _sample_id, is_candidate, _csv -> is_candidate }
    nextclade_join
        .filter { _sample_id, _is_candidate, csv -> csv == null }
        .subscribe { sample_id, _is_candidate, _csv ->
            log.warn "PHYLOGENETICS(${subtype}): no Nextclade result for sample ${sample_id}; its subclade is NA and it is not clustered."
        }
    nextclade_csvs = nextclade_join
        .filter { _sample_id, _is_candidate, csv -> csv != null }
        .map { _sample_id, _is_candidate, csv -> csv }
        .toSortedList { a, b -> a.name <=> b.name }
    def h_tag = subtype.find(/H\d+/)
    dataset = datasets.flatten()
        .filter { dataset_dir -> dataset_dir.name == "nextclade_${h_tag}_dataset" }
        .toList()
        .map { found ->
            if (!found) log.warn "PHYLOGENETICS(${subtype}): no nextclade_${h_tag}_dataset was retrieved by GetDatasets, so subclades, clusters, annotations, the figure and the manifest are not produced."
            found
        }
        .filter { found -> found }
        .map { found -> found[0] }
    PhyloSubclades(nextclade_csvs, PhyloPreflight.out.references, dataset)

    // PhyloPostTree runs once per tree, fanned out over PhyloTree.out.for_post_tree: that emit already carries
    // everything PhyloPostTree needs (tree_id, tree_folder, tree, support, alignment, segments) in the exact shape
    // it takes, so no join by tree_id is needed here. (Chained/grouped joins of several independently-scheduled
    // per-tree channels proved unreliable under load: some trees silently never reached PhyloReport.)
    PhyloPostTree(
        PhyloTree.out.for_post_tree,
        PhyloPreflight.out.gene_map,
        PhyloSubclades.out.subclades,
        params.clusterThreshold
    )

    // Resolved parameters and input checksums for the run manifest. phyloSubtype/Reference/ReferenceFallback/
    // RootStrains/Outgroup are overwritten below with this subtype's resolved values (params holds the raw,
    // possibly multi-subtype, CLI values; see resolvePhyloConfig()).
    def param_names = [
        'protocol', 'inputFasta', 'metadata', 'outDir', 'phyloSubtype', 'phyloReference', 'phyloReferenceFallback',
        'phyloRootStrains', 'phyloOutgroup', 'phyloSegmentOrder', 'phyloMinCoverage', 'model', 'bootstrap', 'seed',
        'clusterThreshold', 'minSupportLabel', 'aaPositions', 'aaPositionsH1n1', 'aaPositionsFile', 'aaGene',
        'minPanelCoverage', 'labMap', 'periodBin', 'width', 'height', 'dpi', 'phyloPanels'
    ]
    def run_params = param_names.collectEntries { name -> [(name): params[name]?.toString()] }
    run_params.phyloSubtype = subtype
    run_params.phyloReference = reference
    run_params.phyloReferenceFallback = reference_fallback
    run_params.phyloRootStrains = root_strains
    run_params.phyloOutgroup = outgroup
    run_params.aaPositionsResolved = aa_positions
    def run_info = [pipeline_version: params.version, started: workflow.start.toString(), parameters: run_params]
    def checksum_inputs = [file(params.inputFasta), ref_genomes, cds_refs, panel_spec]
    if (params.metadata) checksum_inputs << file(params.metadata)
    if (params.labMap) checksum_inputs << file(params.labMap)

    // PhyloReport runs once per tree, fanned out over PhyloPostTree.out.for_report: that emit already carries
    // everything PhyloReport needs (tree_id, tree_folder, tree, support, clusters, proteins, segments) in the exact
    // shape it takes, so no join by tree_id is needed here (see the comment on PhyloPostTree's call above).
    PhyloReport(
        PhyloPostTree.out.for_report,
        PhyloSubclades.out.subclades,
        PhyloPreflight.out.positions,
        metadata.first(),
        lab_map,
        panel_spec,
        PhyloPreflight.out.versions,
        qc_table.first(),
        PhyloPreflight.out.coordinate_table,
        PhyloPreflight.out.root_table,
        checksum_inputs,
        params.periodBin,
        params.minPanelCoverage,
        params.width,
        params.height,
        params.dpi,
        params.minSupportLabel,
        groovy.json.JsonOutput.toJson(run_info).bytes.encodeBase64().toString()
    )

    // Every phylo step below PhyloAlignment/PhyloSegmentAlignment uses errorStrategy 'ignore': a failed task just
    // silently drops its tree, with nothing in pipeline_errors.log and no non-zero exit for onComplete to notice
    // (e.g. "no root strain has segment X", a MAFFT/IQ-TREE/TreeTime crash, or PhyloAlignment's own <3-samples
    // exit). Compare what was actually REQUESTED (segments gated in above, plus the whole-genome tree if enabled)
    // against what reached each later stage, and name the earliest stage a tree never got past.
    //
    // Every "list of ids" channel below is built with .toList(), so each emits exactly one List item. Nextflow's
    // .combine() (like .join()) treats a channel item that is itself a List as a multi-field tuple to spread, not
    // as one opaque value - so combining two such channels directly silently (and wrongly) unpacks each List's
    // elements into separate fields instead of keeping the whole List as a single field. Wrapping each list in a
    // one-key Map first (Maps are never spread this way) keeps every .combine() step unambiguous. (A local closure
    // variable can't be called like a function, e.g. `idMap(x)`, inside a workflow body, so this is inlined.)
    def expected_segment_ids_ch = segment_gated_ch.map { seg, _dirs -> "${subtype}_${seg}".toString() }.toList()
        .map { list -> [ids: list] }
    def expected_wg_ids_ch = channel.of(whole_genome ? subtype : null).filter { it != null }.toList()
        .map { list -> [ids: list] }
    def expected_tree_ids_ch = expected_segment_ids_ch.combine(expected_wg_ids_ch)
        .map { seg_m, wg_m -> seg_m.ids + wg_m.ids }
    def alignment_ok_ids_ch = all_trees_ch.map { tree_id, _f, _al, _og, _seg -> tree_id }.toList()
        .map { list -> [ids: list] }
    def tree_ok_ids_ch      = PhyloTree.out.tree.map { tree_id, _f -> tree_id }.toList()
        .map { list -> [ids: list] }
    def posttree_ok_ids_ch  = PhyloPostTree.out.for_report.map { tree_id, _f, _t, _s, _c, _p, _seg -> tree_id }.toList()
        .map { list -> [ids: list] }
    def report_ok_ids_ch    = PhyloReport.out.html.map { tree_id, _f -> tree_id }.toList()
        .map { list -> [ids: list] }
    def tree_problems_ch = expected_tree_ids_ch.map { list -> [ids: list] }
        .combine(alignment_ok_ids_ch).combine(tree_ok_ids_ch)
        .combine(posttree_ok_ids_ch).combine(report_ok_ids_ch)
        .map { expected_m, alignment_m, tree_m, posttree_m, report_m ->
            def expected = expected_m.ids
            def alignment_ok = alignment_m.ids
            def tree_ok = tree_m.ids
            def posttree_ok = posttree_m.ids
            def report_ok = report_m.ids
            def missing = expected - report_ok
            def lines = missing.collect { tree_id ->
                def stage = !alignment_ok.contains(tree_id) ? "PhyloAlignment/PhyloSegmentAlignment (alignment build - check for 'no root strain has segment', a MAFFT failure, or fewer than 3 samples)" :
                            !tree_ok.contains(tree_id)      ? "PhyloTree (IQ-TREE)" :
                            !posttree_ok.contains(tree_id)  ? "PhyloPostTree (TreeTime/TreeCluster)" :
                                                               "PhyloReport (annotation/figure)"
                def msg = "PHYLOGENETICS(${subtype}): tree '${tree_id}' was requested but never produced a report; ${stage} silently failed for it (errorStrategy 'ignore'). Check .nextflow.log for the failed task, or look for a missing phylogeny/${subtype}/<tree>/ folder."
                log.warn msg
                msg
            }
            lines.join('\n')
        }
        .filter { text -> text }
        .collectFile(name: "phylo_tree_problems_${subtype}.log", newLine: true)

    // Every per-tree emit is a (tree_id, file) tuple (one item per tree built); strip the tag for the flat
    // 'outputs' channel published to <outDir>/phylogeny/<subtype>/. Subtype-level (shared) emits are plain files.
    outputs = qc_table.mix(
        PhyloPreflight.out.versions,
        PhyloPreflight.out.gene_map,
        PhyloPreflight.out.positions,
        PhyloPreflight.out.coordinate_table,
        PhyloPreflight.out.root_table,
        all_trees_ch.map { _tree_id, _folder, alignment, _outgroup_tip, _seg -> alignment },
        all_trees_ch.map { _tree_id, _folder, _al, _outgroup_tip, segments_tsv -> segments_tsv },
        PhyloSegmentAlignment.out.outgroup.map { _seg, f -> f },
        PhyloSegmentAlignment.out.root_used.map { _seg, f -> f },
        PhyloTree.out.raw_tree.map { _tree_id, f -> f },
        PhyloTree.out.tree.map { _tree_id, f -> f },
        PhyloTree.out.support.map { _tree_id, f -> f },
        PhyloTree.out.report.map { _tree_id, f -> f },
        PhyloTree.out.log.map { _tree_id, f -> f },
        PhyloSubclades.out.subclades,
        PhyloPostTree.out.nucleotides.map { _tree_id, f -> f },
        PhyloPostTree.out.annotated_tree.map { _tree_id, f -> f },
        PhyloPostTree.out.proteins.map { _tree_id, f -> f },
        PhyloPostTree.out.clusters.map { _tree_id, f -> f },
        PhyloPostTree.out.mutations.map { _tree_id, f -> f },
        PhyloReport.out.annotations.map { _tree_id, f -> f },
        PhyloReport.out.sources.map { _tree_id, f -> f },
        PhyloReport.out.aa_states.map { _tree_id, f -> f },
        PhyloReport.out.svg.map { _tree_id, f -> f },
        PhyloReport.out.png.map { _tree_id, f -> f },
        PhyloReport.out.panels.map { _tree_id, f -> f },
        PhyloReport.out.manifest.map { _tree_id, f -> f },
        PhyloReport.out.html.map { _tree_id, f -> f }
    )

    emit:
    outputs        = outputs                          // every file published to <outDir>/phylogeny/<subtype>
    errors         = PhyloCoverageQC.out.errors        // per-sample QC errors, merged into pipeline_errors.log
    report         = PhyloReport.out.html.map { _tree_id, f -> f }  // interactive tree(s), added to the index.html dashboard
    tree_problems  = tree_problems_ch                  // 0 or 1 file: requested trees that a later 'ignore'd step silently dropped
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
 * Parse --aa-positions-file: a CSV with (case-insensitive) columns 'protein' and 'position', and an optional
 * 'subtype' column (H3N2 / H1N1pdm09). A row with no subtype column, or a blank subtype cell, applies to every
 * requested subtype; a row with a subtype applies only to it. Returns a Map [subtype: [gene, pos, label] list, ...]
 * with only the subtypes actually touched by at least one row (so a subtype the file never mentions is left for
 * --aa-positions/--aa-positions-h1n1 to fill in, see the PHYLOGENETICS invocation in main.nf).
 * null, empty or a value-less flag (true) means no file was given, i.e. no override ([:]).
 * Syntax/lookup errors throw IllegalArgumentException; range checks need the reference translation and are done
 * by PhyloPreflight, same as --aa-positions.
 */
def parseAaPositionsFile(value) {
    if (value == null || value == true || value.toString().trim() == '') return [:]

    def known_subtypes = phyloKnownSubtypes()
    def gene_pattern = /^[A-Za-z0-9][A-Za-z0-9-]*$/
    def f = file(value.toString())
    if (!f.exists()) {
        throw new IllegalArgumentException("--aa-positions-file: file not found: ${value}")
    }
    def lines = f.readLines().findAll { line -> line.trim() }
    if (!lines) {
        throw new IllegalArgumentException("--aa-positions-file: ${value} is empty.")
    }
    def header = lines[0].split(',', -1).collect { h -> h.trim().toLowerCase() }
    def protein_idx = header.indexOf('protein')
    def position_idx = header.indexOf('position')
    def subtype_idx = header.indexOf('subtype')
    if (protein_idx < 0 || position_idx < 0) {
        throw new IllegalArgumentException("--aa-positions-file: ${value} must have 'protein' and 'position' columns (got header: ${lines[0]}).")
    }

    def result = [:]
    def seen = [:]
    lines.drop(1).eachWithIndex { line, i ->
        def row_num = i + 2 // 1-based, the header is row 1
        def cols = line.split(',', -1)
        if (cols.size() <= Math.max(protein_idx, position_idx)) {
            throw new IllegalArgumentException("--aa-positions-file: row ${row_num} '${line}' is missing the protein or position column.")
        }
        def gene_raw = cols[protein_idx].trim()
        def gene = gene_raw.toUpperCase()
        def pos_str = cols[position_idx].trim()
        def row_subtype = (subtype_idx >= 0 && subtype_idx < cols.size()) ? cols[subtype_idx].trim() : ''
        if (!(gene ==~ gene_pattern)) {
            throw new IllegalArgumentException("--aa-positions-file: row ${row_num} has an invalid protein name '${gene_raw}'.")
        }
        if (!(pos_str ==~ /^\d+$/) || pos_str.toBigInteger() < 1) {
            throw new IllegalArgumentException("--aa-positions-file: row ${row_num} position '${pos_str}' must be a positive integer.")
        }
        def pos = pos_str.toBigInteger().intValue()
        def label = "${gene}:${pos}".toString()
        // Blank/absent subtype cell -> every requested subtype; otherwise only the named one
        def target_subtypes
        if (row_subtype) {
            if (!(row_subtype in known_subtypes)) {
                throw new IllegalArgumentException("--aa-positions-file: row ${row_num} has an unknown subtype '${row_subtype}': only ${known_subtypes.join(', ')} are supported.")
            }
            target_subtypes = [row_subtype]
        } else {
            target_subtypes = known_subtypes
        }
        target_subtypes.each { s ->
            if (!result.containsKey(s)) {
                result[s] = []
                seen[s] = [] as Set
            }
            if (seen[s].contains(label)) {
                throw new IllegalArgumentException("--aa-positions-file: position ${label} (subtype ${s}) is given more than once.")
            }
            seen[s] << label
            result[s] << [gene: gene, pos: pos, label: label]
        }
    }
    return result
}

/*
 * Validate the phylogenetics parameters. Returns a list of error messages (empty when valid).
 */
def validatePhyloParams(Map p) {
    def errors = []

    if (p.protocol != 'HUMAN') {
        errors << "--phylogenetics is only available for the HUMAN protocol (got '${p.protocol}')."
    }
    def known_subtypes = phyloKnownSubtypes()
    def subtypes = (p.phyloSubtype ?: '').toString().split(',').collect { s -> s.trim() }.findAll { s -> s }
    def unknown = subtypes.findAll { s -> !(s in known_subtypes) }
    if (!subtypes) {
        errors << "--phyloSubtype must list at least one subtype (${known_subtypes.join(', ')})."
    } else if (unknown) {
        errors << "--phyloSubtype '${p.phyloSubtype}' has unknown subtype(s) ${unknown.join(', ')}: only ${known_subtypes.join(', ')} are supported."
    }
    // Per-subtype reference / fallback / root strains / outgroup bundle (resolvePhyloConfig(): flat H3N2 overrides,
    // phyloSubtypeConfig, then the built-in defaults)
    subtypes.findAll { s -> s in known_subtypes }.each { s ->
        def cfg = resolvePhyloConfig(s, p)
        if (!cfg.reference) {
            errors << "${s}: no reference strain configured (set phyloSubtypeConfig.${s}.reference)."
        }
        def roots = cfg.rootStrains ? cfg.rootStrains.toString().split(',').collect { r -> r.trim() }.findAll { r -> r } : []
        if (!roots) {
            errors << "${s}: no root strains configured (set phyloSubtypeConfig.${s}.rootStrains)."
        }
        if (!(cfg.outgroup in roots)) {
            if (s == 'H3N2') {
                // Kept as the historical H3N2 message so existing --phyloOutgroup/--phyloRootStrains users see it unchanged
                errors << "--phyloOutgroup '${cfg.outgroup}' must be one of --phyloRootStrains."
            } else {
                errors << "${s}: outgroup '${cfg.outgroup}' must be one of its root strains (${cfg.rootStrains})."
            }
        }
    }
    [ 'phyloMinCoverage', 'minPanelCoverage' ].each { k ->
        if (!isPhyloNumber(p[k]) || p[k].toString().toBigDecimal() < 0 || p[k].toString().toBigDecimal() > 1) {
            errors << "--${k} must be a number between 0 and 1 (got '${p[k]}')."
        }
    }
    if (!isPhyloNumber(p.clusterThreshold) || p.clusterThreshold.toString().toBigDecimal() <= 0) {
        errors << "--clusterThreshold must be a positive number of substitutions per site (got '${p.clusterThreshold}')."
    }
    if (!isPhyloNumber(p.minSupportLabel) || p.minSupportLabel.toString().toBigDecimal() < 0 || p.minSupportLabel.toString().toBigDecimal() > 100) {
        errors << "--minSupportLabel must be a number between 0 and 100 (got '${p.minSupportLabel}')."
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
    // --phyloWholeGenome / --phyloSegmentTrees: at least one of the whole-genome tree or a segment tree must be built
    def whole_genome = p.phyloWholeGenome != null && p.phyloWholeGenome != false && p.phyloWholeGenome.toString().toLowerCase() != 'false'
    def requested_trees = (p.phyloSegmentTrees ?: 'none').toString().trim()
    def segment_trees = []
    if (requested_trees.toLowerCase() == 'none') {
        segment_trees = []
    } else if (requested_trees.toLowerCase() == 'all') {
        segment_trees = order
    } else {
        segment_trees = requested_trees.split(',').collect { s -> s.trim() }.findAll { s -> s }
        def bad_segments = segment_trees.findAll { s -> !(s in order) }
        if (!segment_trees) {
            errors << "--phyloSegmentTrees must be 'all', 'none' or a comma-separated list of segments (${order})."
        } else if (bad_segments) {
            errors << "--phyloSegmentTrees has unknown segment(s) ${bad_segments.join(', ')}: must be one of ${order}."
        }
    }
    if (!whole_genome && !segment_trees) {
        errors << "--phyloWholeGenome is false and --phyloSegmentTrees is 'none': nothing would be built."
    }
    if (p.labMap && !file(p.labMap.toString()).exists()) {
        errors << "--labMap file not found: ${p.labMap}"
    }
    try {
        parseAaPositions(p.aaPositions, p.aaGene)
    } catch (IllegalArgumentException e) {
        errors << e.message
    }
    try {
        parseAaPositions(p.aaPositionsH1n1, p.aaGene)
    } catch (IllegalArgumentException e) {
        errors << "--aa-positions-h1n1: ${e.message}"
    }
    try {
        parseAaPositionsFile(p.aaPositionsFile)
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

// Subtypes the phylogenetics module can build a tree for. Keep in sync with MergeReports.nf's
// PHYLO_KNOWN_SUBTYPES (Groovy/Python, so not shareable directly): both lists gate which subtypes the pipeline
// (this file) and the dashboard (MergeReports.nf) recognise.
def phyloKnownSubtypes() {
    return ['H3N2', 'H1N1pdm09']
}

// Built-in reference / fallback / root strains / outgroup per subtype, mirroring nextflow.config's
// params.phyloSubtypeConfig. Used as the fallback when params.phyloSubtypeConfig is not set (e.g. in unit tests
// that call these functions directly without loading nextflow.config). Keep in sync with nextflow.config.
def phyloSubtypeDefaults() {
    return [
        H3N2: [reference: "A/Darwin/6/2021", referenceFallback: "A/Massachusetts/18/2022",
               rootStrains: "A/Darwin/6/2021,A/Massachusetts/18/2022,A/Croatia/10136RV/2023,A/Singapore/GP20238/2024,A/Sydney/1359/2024",
               outgroup: "A/Darwin/6/2021"],
        H1N1pdm09: [reference: "A/Wisconsin/67/2022", referenceFallback: "A/Wisconsin/67/2022",
               rootStrains: "A/Victoria/2570/2019,A/Sydney/5/2021,A/Wisconsin/67/2022",
               outgroup: "A/Victoria/2570/2019"]
    ]
}

// GenotypingInfo_ch tag for `subtype` (e.g. "H1N1pdm09" -> "H1N1"): SubtypeDetection emits e.g. "A(H1N1)pdm09" ->
// h_tag "H1", n_tag "N1", so H1N1pdm09 samples carry tag "H1N1". Throws for an unknown subtype.
def phyloSubtypeTag(String subtype) {
    def tags = [H3N2: 'H3N2', H1N1pdm09: 'H1N1']
    if (!(subtype in tags.keySet())) {
        throw new IllegalArgumentException("Unknown phylogenetics subtype '${subtype}': only ${phyloKnownSubtypes().join(', ')} are supported.")
    }
    return tags[subtype]
}

/*
 * Resolve the reference / referenceFallback / rootStrains / outgroup bundle for `subtype`, in priority order:
 *   1. the flat phyloReference/phyloReferenceFallback/phyloRootStrains/phyloOutgroup params (H3N2 only, kept for
 *      backward compatibility: a value given there always wins for H3N2)
 *   2. p.phyloSubtypeConfig[subtype] (nextflow.config, or a -params-file override)
 *   3. phyloSubtypeDefaults()[subtype]
 */
def resolvePhyloConfig(String subtype, Map p) {
    def base = phyloSubtypeDefaults()[subtype] ?: [:]
    def configured = (p.phyloSubtypeConfig ?: [:])[subtype] ?: [:]
    def resolved = [
        reference:         base.reference,
        referenceFallback: base.referenceFallback,
        rootStrains:       base.rootStrains,
        outgroup:          base.outgroup,
    ]
    configured.each { k, v -> if (v != null) resolved[k] = v }
    if (subtype == 'H3N2') {
        if (p.phyloReference)         resolved.reference = p.phyloReference
        if (p.phyloReferenceFallback) resolved.referenceFallback = p.phyloReferenceFallback
        if (p.phyloRootStrains)       resolved.rootStrains = p.phyloRootStrains
        if (p.phyloOutgroup)          resolved.outgroup = p.phyloOutgroup
    }
    return resolved
}
