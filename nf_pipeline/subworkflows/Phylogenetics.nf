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
include { PhyloLog              } from '../modules/phylogeny/PhyloLog'
include { PhyloAuspice          } from '../modules/phylogeny/PhyloAuspice'

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
    auspice_geography      // --auspice: value channel [geo.tsv, lat_longs.tsv, geography.log] from PhyloGeography (run ONCE per
                           // run in main.nf and shared by every subtype and tree); channel.empty() without --auspice

    main:
    def resources = params.protocols[params.protocol].resources
    def ref_genomes = file("${resources}/${params.protocol}_references.fasta", checkIfExists: true)
    def cds_refs = file("${resources}/CDS_references.fasta", checkIfExists: true)
    def panel_spec = file(params.phyloPanels, checkIfExists: true)
    def lab_map = params.labMap ? file(params.labMap, checkIfExists: true) : []
    def auspice_on = isPhyloFlag(params.auspice)

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
    def whole_genome_log_ch = channel.empty()
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
        whole_genome_log_ch = PhyloAlignment.out.log.map { f -> tuple(subtype, "whole_genome", f) }
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
    // (No skipped-segment warnings when there are no candidate samples at all: "no samples were found" says it.)
    segment_passing_ch
        .filter { _seg, dirs -> dirs.size() < 3 }
        .combine(candidates.count())
        .filter { _seg, _dirs, n_candidates -> n_candidates > 0 }
        .subscribe { seg, dirs, _n_candidates ->
            def reason = dirs.size() == 0 ? "0 samples passed coverage QC for ${seg}" : "only ${dirs.size()} sample(s) passed coverage QC for ${seg}"
            def msg = "PHYLOGENETICS(${subtype}): segment tree ${seg} skipped: ${reason} (need at least 3)."
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
    // (Not when there are no candidate samples at all: the "no samples were found" warning above already says so.)
    all_trees_ch.count().combine(candidates.count()).subscribe { n_trees, n_candidates ->
        if (n_trees == 0 && n_candidates > 0) println "WARN: PHYLOGENETICS(${subtype}): no tree was built (check --phyloWholeGenome/--phyloSegmentTrees and whether enough samples passed coverage QC)."
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
            println "WARN: PHYLOGENETICS(${subtype}): no Nextclade result for sample ${sample_id}; its subclade is NA and it is not clustered."
        }
    nextclade_csvs = nextclade_join
        .filter { _sample_id, _is_candidate, csv -> csv != null }
        .map { _sample_id, _is_candidate, csv -> csv }
        .toSortedList { a, b -> a.name <=> b.name }
    def h_tag = subtype.find(/H\d+/)
    // `dataset` must stay a value channel (toList/filter/map only): PhyloSubclades and everything fed by its outputs
    // (PhyloPostTree, PhyloReport, PhyloLog) run once per tree, and a queue channel here would make them run once in
    // total. Hence the missing-dataset warning is a separate subscribe, not part of this chain.
    dataset = datasets.flatten()
        .filter { dataset_dir -> dataset_dir.name == "nextclade_${h_tag}_dataset" }
        .toList()
        .filter { found -> found }
        .map { found -> found[0] }
    // (Not when there are no candidate samples at all: "no samples were found" already says so.)
    datasets.flatten()
        .filter { dataset_dir -> dataset_dir.name == "nextclade_${h_tag}_dataset" }
        .toList()
        .map { found -> [found: found] }
        .combine(candidates.count())
        .subscribe { found_m, n_candidates ->
            if (!found_m.found && n_candidates > 0) println "WARN: PHYLOGENETICS(${subtype}): no nextclade_${h_tag}_dataset was retrieved by GetDatasets, so subclades, clusters, annotations, the figure and the manifest are not produced."
        }
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
    // everything PhyloReport needs (tree_id, tree_folder, tree, support, clusters, proteins, segments, plus the
    // alignment, ancestral_nt and gene_map rows it only passes on) in the exact shape it takes, so no join by
    // tree_id is needed here (see the comment on PhyloPostTree's call above). Its `for_auspice` emit (one bundle per
    // tree) feeds PhyloAuspice below when --auspice is set; nothing of it is published.
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

    // Auspice export (--auspice, plan Phase 5): one Nextstrain dataset per tree, fanned out over PhyloReport's own
    // for_auspice bundle (no join by tree_id). The run-level inputs are value channels / plain values: PhyloGeography's
    // geo.tsv and lat_longs.tsv (resolved once per run in main.nf), the metadata, PhyloPreflight's root_strains.tsv.
    def auspice_json_ch = channel.empty()
    def auspice_log_parts_ch = channel.empty()
    def auspice_status_ch = channel.value([items: null])   // null = --auspice off (tree_problems/summary skip it)
    if (auspice_on) {
        PhyloAuspice(
            PhyloReport.out.for_auspice,
            auspice_geography.map { geo, _ll, _log -> geo },
            auspice_geography.map { _geo, ll, _log -> ll },
            metadata.first(),
            file(params.auspiceConfig, checkIfExists: true),
            file(params.auspiceClockRates, checkIfExists: true),
            PhyloPreflight.out.root_table,
            params.auspiceRootDates ? file(params.auspiceRootDates, checkIfExists: true) : [],
            params.auspiceMinDatedTips,
            auspiceTraitList(params.auspiceTraits).join(','),
            params.version,
            params.minSupportLabel
        )
        auspice_json_ch = PhyloAuspice.out.json.map { _tree_id, f -> f }
        // phylo.log parts: the run-level geography.log (the same file for every tree whose report exists) and the tree's
        // own auspice.log. geography.log is attached to for_auspice, not to PhyloAuspice's output, so a PhyloAuspice task
        // that never ran still leaves a visible "did not finish" note in that tree's phylo.log.
        auspice_log_parts_ch = PhyloReport.out.for_auspice
            .map { tree_id, tree_folder, _t, _a, _n, _p, _s, _g, _an, _su -> tuple(tree_id, tree_folder) }
            .combine(auspice_geography.map { _geo, _ll, log -> log })
            .mix(PhyloAuspice.out.step_log.map { tree_id, f -> tuple(tree_id, tree_id == subtype ? "whole_genome" : tree_id.substring(subtype.length() + 1), f) })
        // Each status file: tree_id status mode dataset reason (one row). Read here, wrapped in a one-key Map for combine().
        auspice_status_ch = PhyloAuspice.out.status
            .map { tree_id, _folder, f ->
                def row = f.readLines()[1].split('\t', -1)
                [id: tree_id, status: row[1], mode: row[2], dataset: row[3], reason: row[4]]
            }
            .toList()
            .map { list -> [items: list] }
    }

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
    def posttree_ok_ids_ch  = PhyloPostTree.out.for_report.map { tree_id, _f, _t, _s, _c, _p, _seg, _aln, _nt, _gm -> tree_id }.toList()
        .map { list -> [ids: list] }
    def report_ok_ids_ch    = PhyloReport.out.html.map { tree_id, _f -> tree_id }.toList()
        .map { list -> [ids: list] }
    def tree_problems_ch = expected_tree_ids_ch.map { list -> [ids: list] }
        .combine(alignment_ok_ids_ch).combine(tree_ok_ids_ch)
        .combine(posttree_ok_ids_ch).combine(report_ok_ids_ch).combine(auspice_status_ch)
        .map { expected_m, alignment_m, tree_m, posttree_m, report_m, auspice_m ->
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
                def msg = "PHYLOGENETICS(${subtype}): tree '${tree_id}' was requested but never produced a report; ${stage} silently failed for it (errorStrategy 'ignore'). Check phylogeny/${subtype}/<tree>/phylo.log (present when any step of that tree ran) and .nextflow.log for the failed task."
                println "WARN: ${msg}"
                msg
            }
            // --auspice: every tree that has a report must also have a dataset. PhyloAuspice reports its own expected
            // failures in its status (status failed + reason, exit 0); a tree with no status at all means the task died
            // before its Python ran (e.g. the Nextstrain conda env could not be built) or PhyloGeography failed.
            if (auspice_m.items != null) {
                expected.findAll { tree_id -> report_ok.contains(tree_id) }.each { tree_id ->
                    def st = auspice_m.items.find { it.id == tree_id }
                    def folder = tree_id == subtype ? "whole_genome" : tree_id.substring(subtype.length() + 1)
                    def msg = null
                    if (st == null) {
                        msg = "PHYLOGENETICS(${subtype}): tree '${tree_id}' has no Auspice dataset: PhyloAuspice never reported a status (the task failed before it ran, e.g. the Nextstrain conda env could not be built, or PhyloGeography failed). Check .nextflow.log."
                    } else if (st.status != 'ok') {
                        msg = "PHYLOGENETICS(${subtype}): tree '${tree_id}' has no Auspice dataset: PhyloAuspice failed: ${st.reason} (see phylogeny/${subtype}/${folder}/phylo.log)."
                    }
                    if (msg) {
                        println "WARN: ${msg}"
                        lines << msg
                    }
                }
            }
            lines.join('\n')
        }
        .filter { text -> text }
        .collectFile(name: "phylo_tree_problems_${subtype}.log", newLine: true)

    // Per-tree phylo.log + tree.nwk. Every phylogenetics step writes its messages to a small part file instead of
    // the terminal; the parts of one tree (alignment, IQ-TREE log + treefile, TreeTime/clusters, report) are grouped
    // by tree_id and assembled by PhyloLog together with the subtype-level preflight/QC/subclade notes. groupTuple
    // never drops a tree: a group that is still incomplete when the steps are done (a later step failed) is flushed
    // with the parts that exist, so its phylo.log shows where it stopped.
    // tree_folder is recovered from tree_id ("H3N2" -> whole_genome, "H3N2_HA" -> HA) for the steps that only carry the id.
    def log_parts_ch = whole_genome_log_ch
        .mix(PhyloSegmentAlignment.out.log.map { seg, f -> tuple("${subtype}_${seg}".toString(), seg, f) })
        .mix(PhyloTree.out.step_log.map { tree_id, f -> tuple(tree_id, tree_id == subtype ? "whole_genome" : tree_id.substring(subtype.length() + 1), f) })
        .mix(PhyloTree.out.raw_tree.map { tree_id, f -> tuple(tree_id, tree_id == subtype ? "whole_genome" : tree_id.substring(subtype.length() + 1), f) })
        .mix(PhyloPostTree.out.step_log.map { tree_id, f -> tuple(tree_id, tree_id == subtype ? "whole_genome" : tree_id.substring(subtype.length() + 1), f) })
        .mix(PhyloReport.out.step_log.map { tree_id, f -> tuple(tree_id, tree_id == subtype ? "whole_genome" : tree_id.substring(subtype.length() + 1), f) })
        .mix(auspice_log_parts_ch)
        // 5 parts per tree, + geography.log and auspice.log with --auspice
        .groupTuple(by: [0, 1], size: auspice_on ? 7 : 5, remainder: true)
    // qc_table.first(): collectFile emits one item on a queue channel, which would make PhyloLog run only once
    PhyloLog(log_parts_ch, PhyloPreflight.out.log, PhyloSubclades.out.log, qc_table.first(), params.phyloSegmentOrder)

    // One terminal line per subtype once all its trees are done. Built = has a report (like tree_problems above).
    // Lists are wrapped in one-key Maps before combine(), see the note above tree_problems_ch.
    def built_ch = PhyloLog.out.summary
        .map { tree_id, tree_folder, n -> [id: tree_id, folder: tree_folder, n: n.toString().trim()] }
        .toList()
        .map { list -> [items: list] }
    def skipped_ch = segment_passing_ch
        .filter { _seg, dirs -> dirs.size() < 3 }
        .map { seg, _dirs -> seg }
        .toList()
        .map { list -> [ids: list] }
    built_ch.combine(report_ok_ids_ch).combine(skipped_ch)
        .map { built_m, report_m, skipped_m ->
            def order = params.phyloSegmentOrder.toString().split(',').collect { s -> s.trim() }
            def built = built_m.items.findAll { it.id in report_m.ids }
            def wg = built.find { it.folder == 'whole_genome' }
            def segs = built.findAll { it.folder != 'whole_genome' }.sort { a, b -> order.indexOf(a.folder) <=> order.indexOf(b.folder) }
            def parts = []
            if (wg) parts << "whole-genome tree built (${wg.n} samples)"
            if (segs) parts << "${segs.size()} segment tree${segs.size() == 1 ? '' : 's'} built (${segs.collect { "${it.folder} ${it.n}" }.join(', ')})"
            def line = parts ? "PHYLOGENETICS(${subtype}): ${parts.join(' and ')} -> phylogeny/${subtype}/" : ''
            if (line && skipped_m.ids) line += " (skipped: ${skipped_m.ids.sort { order.indexOf(it) }.join(', ')})"
            line
        }
        .filter { line -> line }
        .subscribe { line -> println line }

    // --auspice: one terminal line per subtype, e.g. "AUSPICE(H3N2): 8 datasets -> auspice/ (view: npx auspice@3.0.0 view
    // <outDir>/auspice)", plus the divergence-only and failed counts (failed = status failed, or no status at all)
    if (auspice_on) {
        report_ok_ids_ch.combine(auspice_status_ch)
            .map { report_m, auspice_m ->
                if (!report_m.ids) return ''
                def items = auspice_m.items.findAll { it.id in report_m.ids }
                def ok = items.findAll { it.status == 'ok' }
                def divergence = ok.count { it.mode != 'timetree' }
                def failed = report_m.ids.size() - ok.size()
                def line = "AUSPICE(${subtype}): ${ok.size()} dataset${ok.size() == 1 ? '' : 's'} -> auspice/ (view: npx auspice@3.0.0 view ${params.outDir}/auspice)"
                if (divergence) line += "; ${divergence} divergence-only (no time tree, see phylo.log)"
                if (failed) line += "; ${failed} failed (see pipeline_errors.log)"
                line
            }
            .filter { line -> line }
            .subscribe { line -> println line }
    }

    // Published to <outDir>/phylogeny/<subtype>/<tree_folder>/: only the tree, its log and the HTML report.
    // Everything else stays in work/ and is still used internally exactly as before.
    outputs = PhyloLog.out.tree.map { _tree_id, f -> f }.mix(
        PhyloLog.out.log.map { _tree_id, f -> f },
        PhyloReport.out.html.map { _tree_id, f -> f }
    )

    emit:
    outputs        = outputs                          // tree.nwk, phylo.log and the HTML of every tree, published to <outDir>/phylogeny/<subtype>
    auspice        = auspice_json_ch                  // --auspice: validated flutyper_<subtype>_<segment>.json files, published to <outDir>/auspice
    errors         = PhyloCoverageQC.out.errors        // per-sample QC errors, merged into pipeline_errors.log
    report         = PhyloReport.out.html.map { _tree_id, f -> f }  // interactive tree(s), added to the index.html dashboard
    tree_problems  = tree_problems_ch                  // 0 or 1 file: requested trees that a later 'ignore'd step silently dropped (incl. a missing Auspice dataset)
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
    errors.addAll(validateAuspiceParams(p))
    return errors
}

/*
 * Validate the --auspice parameters (only when --auspice is set; --auspice without --phylogenetics is rejected in
 * main.nf before this runs). Only existence and syntax are checked here: the contents of the clock-rate and root-date
 * tables are checked per tree by PhyloAuspice, which reports a bad row as that tree's failure.
 */
def validateAuspiceParams(Map p) {
    def errors = []
    if (!isPhyloFlag(p.auspice)) return errors
    [ 'auspiceConfig', 'auspiceClockRates' ].each { k ->
        if (!p[k] || p[k] == true || !file(p[k].toString()).isFile()) {
            errors << "--${k} file not found: ${p[k]}"
        }
    }
    [ 'auspiceRootDates', 'auspiceLatLongs' ].each { k ->
        if (p[k] && (p[k] == true || !file(p[k].toString()).isFile())) {
            errors << "--${k} file not found: ${p[k]}"
        }
    }
    if (!(p.auspiceMinDatedTips.toString() ==~ /^\d+$/)) {
        errors << "--auspiceMinDatedTips must be a non-negative integer (got '${p.auspiceMinDatedTips}')."
    }
    def trait_columns = ['region', 'country', 'division', 'location', 'subclade', 'cluster', 'season', 'age_group', 'sex']
    def traits = auspiceTraitList(p.auspiceTraits)
    def bad_traits = traits == null ? [] : traits.findAll { t -> !(t in trait_columns) }
    if (traits == null) {
        errors << "--auspiceTraits must be 'none' or a comma-separated list of ${trait_columns.join(', ')} (got '${p.auspiceTraits}')."
    } else if (bad_traits) {
        errors << "--auspiceTraits has unknown column(s) ${bad_traits.join(', ')}: must be among ${trait_columns.join(', ')} (or 'none')."
    }
    def continents = ['AF', 'AN', 'AS', 'EU', 'NA', 'OC', 'SA']
    def codes = (p.auspiceContinent && p.auspiceContinent != true) ? p.auspiceContinent.toString().split(',').collect { c -> c.trim().toUpperCase() }.findAll { c -> c } : []
    if (!codes || codes.any { c -> !(c in continents) }) {
        errors << "--auspiceContinent must be a comma-separated list of GeoNames continent codes ${continents.join(', ')} (got '${p.auspiceContinent}')."
    }
    return errors
}

// --auspiceTraits as a list of column names: 'none' (or empty) -> [] (no DTA); a value-less flag -> null (invalid)
def auspiceTraitList(v) {
    if (v == true) return null
    if (v == null || v.toString().trim() == '' || v.toString().trim().equalsIgnoreCase('none')) return []
    return v.toString().split(',').collect { t -> t.trim() }.findAll { t -> t }
}

// True for a boolean-like flag that is set (true, 'true'); false for null/false/'false'
def isPhyloFlag(v) {
    return v != null && v.toString().toLowerCase() == 'true'
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
