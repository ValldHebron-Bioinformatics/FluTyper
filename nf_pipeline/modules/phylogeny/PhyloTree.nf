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
    debug true

    input:
    path(alignment)
    val(outgroup_tip)
    val(model)
    val(bootstrap)
    val(seed)

    output:
    path("phylo.treefile"),   emit: raw_tree
    path("phylo_tree.nwk"),   emit: tree
    path("node_support.tsv"), emit: support
    path("phylo.iqtree"),     emit: report
    path("phylo.log"),        emit: log

    script:
    def bootstrap_opt = bootstrap.toString().toInteger() > 0 ? "-B ${bootstrap}" : ""
    """
    # IQ-TREE aborts when -T exceeds the machine's cores, so cap it; a different thread count can change the tree slightly
    threads=${task.cpus}
    cores=\$(nproc)
    if [ "\$threads" -gt "\$cores" ]; then
        echo "PhyloTree: WARNING ${task.cpus} cpus requested but only \$cores cores available; IQ-TREE runs with \$cores threads."
        threads=\$cores
    fi
    iqtree2 -s "${alignment}" -m "${model}" ${bootstrap_opt} -seed ${seed} -T \$threads \\
        -o "${outgroup_tip}" --prefix phylo -quiet
    # Surface IQ-TREE warnings (e.g. identical sequences) that would otherwise stay in phylo.log
    grep -E "^WARNING" phylo.log | sed "s/^/PhyloTree: IQ-TREE /" || true

    python3 - <<'PYEOF'
from Bio import Phylo

tree = Phylo.read("phylo.treefile", "newick")
rows = []
for i, clade in enumerate(tree.find_clades(order="preorder")):
    if clade.is_terminal():
        continue
    clade.name = f"NODE_{i:07d}"
    rows.append((clade.name, "" if clade.confidence is None else f"{clade.confidence:g}"))
    clade.confidence = None
Phylo.write(tree, "phylo_tree.nwk", "newick", format_branch_length="%.10g")
with open("node_support.tsv", "w") as f:
    f.write("node\\tsupport\\n")
    for name, support in rows:
        f.write(f"{name}\\t{support}\\n")
PYEOF
    """
}
