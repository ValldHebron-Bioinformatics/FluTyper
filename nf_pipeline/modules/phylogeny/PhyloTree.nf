#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloTree {
    // Maximum likelihood tree of the whole-genome alignment with IQ-TREE (--model, default GTR+G4) and ultrafast
    // bootstrap support (--bootstrap, default 1000; 0 disables it). Branch lengths are substitutions per site and
    // the tree is rooted on the outgroup root strain (--phyloOutgroup). Internal nodes are then named NODE_0000001...
    // in preorder so that TreeTime, TreeCluster and the figure refer to the same nodes; support values are kept in
    // node_support.tsv. A fixed --seed and thread count make the run reproducible;
    // the thread count is capped at the cores available, with a warning when that lowers it.
    label 'big_mem'
    errorStrategy 'ignore'
    // Cache on file content: PhyloPreflight reruns when --aa-positions changes, but its reference outputs do not
    cache 'deep'

    input:
    // tree_id e.g. "H3N2" (whole genome) or "H3N2_HA" (segment tree): globally unique, used for filenames downstream.
    // tree_folder is where THIS tree's files live under the subtype's own output folder: "whole_genome" or the bare
    // segment name e.g. "HA" (see PHYLOGENETICS in subworkflows/Phylogenetics.nf). segments_tsv is only threaded
    // through (unused here) so PhyloPostTree can receive it straight from `for_post_tree` below, with no separate
    // join by tree_id needed: chained/grouped joins of several independently-scheduled per-tree channels proved
    // unreliable (some trees silently missing downstream), so every step instead re-emits what it was given.
    tuple val(tree_id), val(tree_folder), path(alignment), val(outgroup_tip), path(segments_tsv)
    val(model)
    val(bootstrap)
    val(seed)

    output:
    tuple val(tree_id), path("${tree_folder}/phylo.treefile"),   emit: raw_tree
    tuple val(tree_id), path("${tree_folder}/phylo_tree.nwk"),   emit: tree
    tuple val(tree_id), path("${tree_folder}/node_support.tsv"), emit: support
    tuple val(tree_id), path("${tree_folder}/phylo.iqtree"),     emit: report
    tuple val(tree_id), path("${tree_folder}/phylo.log"),        emit: log       // IQ-TREE's own log (internal; not the published phylo.log)
    tuple val(tree_id), path("iqtree.log"),                      emit: step_log  // this step's messages, collected into the tree's phylo.log
    // Everything PhyloPostTree needs, in the exact shape it takes: no join required to call it
    tuple val(tree_id), val(tree_folder), path("${tree_folder}/phylo_tree.nwk"), path("${tree_folder}/node_support.tsv"), path(alignment), path(segments_tsv), emit: for_post_tree

    script:
    def bootstrap_opt = bootstrap.toString().toInteger() > 0 ? "-B ${bootstrap}" : ""
    """
    # Step messages go to iqtree.log (part of the tree's phylo.log), not the terminal
    : > iqtree.log
    mkdir -p "${tree_folder}"
    # IQ-TREE aborts when -T exceeds the machine's cores, so cap it; a different thread count can change the tree slightly
    threads=${task.cpus}
    cores=\$(nproc)
    if [ "\$threads" -gt "\$cores" ]; then
        echo "PhyloTree: WARNING ${task.cpus} cpus requested but only \$cores cores available; IQ-TREE runs with \$cores threads." >> iqtree.log
        threads=\$cores
    fi
    iqtree2 -s "${alignment}" -m "${model}" ${bootstrap_opt} -seed ${seed} -T \$threads \\
        -o "${outgroup_tip}" --prefix "${tree_folder}/phylo" -quiet
    # Surface IQ-TREE warnings (e.g. identical sequences) that would otherwise stay in phylo.log
    grep -E "^WARNING" "${tree_folder}/phylo.log" | sed "s/^/PhyloTree: IQ-TREE /" >> iqtree.log || true

    python3 - <<'PYEOF'
from Bio import Phylo

tree = Phylo.read("${tree_folder}/phylo.treefile", "newick")
rows = []
for i, clade in enumerate(tree.find_clades(order="preorder")):
    if clade.is_terminal():
        continue
    clade.name = f"NODE_{i:07d}"
    rows.append((clade.name, "" if clade.confidence is None else f"{clade.confidence:g}"))
    clade.confidence = None
Phylo.write(tree, "${tree_folder}/phylo_tree.nwk", "newick", format_branch_length="%.10g")
with open("${tree_folder}/node_support.tsv", "w") as f:
    f.write("node\\tsupport\\n")
    for name, support in rows:
        f.write(f"{name}\\t{support}\\n")
PYEOF
    """
}
