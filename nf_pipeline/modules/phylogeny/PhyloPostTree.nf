#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloPostTree {
    // Analyses that run on the finished ML tree. It takes no figure option, so changing --aa-positions or the
    // panels does not rerun it.
    //  1. Ancestral state reconstruction with TreeTime (parsimony) on the whole-genome alignment. The protein of every
    //     node is obtained by translating its nucleotide sequence through the reference codon map (gene_map.tsv), so
    //     internal nodes carry amino acid states too. Tips keep their observed sequence: ambiguous or missing codons
    //     stay X / '-' and are never replaced by inferred states.
    //  2. Genetic-distance clustering within each subclade: the tree is pruned to the samples of one subclade and
    //     TreeCluster single-linkage groups tips whose patristic distances chain below --clusterThreshold. Tips that
    //     do not join a cluster stay 'unclustered'; samples without a subclade and root strains are not clustered (NA).
    //     Labels are run-specific (<subclade>_C<n>, numbered by size), not a nomenclature. The amino acid changes on
    //     the branch into each cluster's common ancestor are listed with how exclusive they are within the subclade.
    errorStrategy 'ignore'
    // Cache on file content: PhyloPreflight reruns when --aa-positions changes, but its reference outputs do not
    cache 'deep'
    debug true

    input:
    path(tree)          // phylo_tree.nwk with named internal nodes
    path(alignment)     // phylo_alignment.fasta
    path(gene_map)      // gene_map.tsv
    path(segments)      // phylo_segments.tsv
    path(subclades)     // subclades.tsv
    val(threshold)      // --clusterThreshold

    output:
    path("ancestral_nt.fasta"),    emit: nucleotides
    path("annotated_tree.nexus"),  emit: annotated_tree
    path("ancestral_aa"),          emit: proteins
    path("clusters.tsv"),          emit: clusters
    path("cluster_mutations.tsv"), emit: mutations

    script:
    """
    # ---- 1. Ancestral reconstruction ----
    treetime ancestral --aln "${alignment}" --tree "${tree}" --method-anc parsimony --outdir treetime_out > treetime.log 2>&1 \\
        || { cat treetime.log >&2; exit 1; }
    grep -iE "warn" treetime.log | sed "s/^/PhyloPostTree: TreeTime /" || true
    cp treetime_out/ancestral_sequences.fasta ancestral_nt.fasta
    cp treetime_out/annotated_tree.nexus annotated_tree.nexus

    python3 - <<'PYEOF'
import csv, os
from Bio import SeqIO
from Bio.Data import CodonTable

observed = {r.id: str(r.seq).upper() for r in SeqIO.parse("${alignment}", "fasta")}
reconstructed = {r.id: str(r.seq).upper() for r in SeqIO.parse("ancestral_nt.fasta", "fasta")}
# Tips: observed sequences. Internal nodes: TreeTime reconstruction.
nodes = {name: seq for name, seq in reconstructed.items() if name not in observed}
nodes.update(observed)
is_tip = {name: name in observed for name in nodes}

seg_start = {}
with open("${segments}") as f:
    for row in csv.DictReader(f, delimiter="\\t"):
        seg_start[row["segment"]] = int(row["start"])

codons = {}  # gene -> list of 0-based alignment column triplets, in protein order
with open("${gene_map}") as f:
    for row in csv.DictReader(f, delimiter="\\t"):
        if row["ref_aa"] == "?":
            cols = None
        else:
            base = seg_start[row["segment"]] - 2  # 1-based segment coordinate -> 0-based alignment column
            cols = tuple(base + int(row[k]) for k in ("nt1", "nt2", "nt3"))
        codons.setdefault(row["gene"], []).append(cols)

# Standard genetic code as a lookup table: only unambiguous ACGT codons are translated
table = CodonTable.unambiguous_dna_by_id[1]
CODE = dict(table.forward_table)
CODE.update({stop: "*" for stop in table.stop_codons})

def translate(seq, cols):
    if cols is None:
        return "X"
    codon = seq[cols[0]] + seq[cols[1]] + seq[cols[2]]
    if codon == "---":
        return "-"
    return CODE.get(codon, "X")

order = sorted(nodes, key=lambda n: (not is_tip[n], n))
os.makedirs("ancestral_aa", exist_ok=True)
for gene, cols in codons.items():
    with open(f"ancestral_aa/{gene}.fasta", "w") as out:
        for n in order:
            seq = nodes[n]
            out.write(f">{n}\\n{''.join(translate(seq, c) for c in cols)}\\n")
PYEOF

    # ---- 2. Clusters within each subclade ----
    python3 - <<'PYEOF'
import csv, glob, os, subprocess, sys
import treeswift
from Bio import Phylo, SeqIO

threshold = "${threshold}"
full = Phylo.read("${tree}", "newick")
tips = {t.name for t in full.get_terminals()}

subclade_of, is_sample = {}, {}
with open("${subclades}") as f:
    for row in csv.DictReader(f, delimiter="\\t"):
        subclade_of[row["tip"]] = row["subclade"]
        is_sample[row["tip"]] = row["source"] == "sample"

by_subclade = {}
for tip in sorted(tips):
    sc = subclade_of.get(tip, "NA")
    if is_sample.get(tip, False) and sc != "NA":
        by_subclade.setdefault(sc, []).append(tip)

ts_tree = treeswift.read_tree_newick("${tree}")
cluster_of, size_of = {}, {}
for sc in sorted(by_subclade):
    members = by_subclade[sc]
    groups = {}
    if len(members) >= 2:
        sub = ts_tree.extract_tree_with(set(members))
        with open("subtree.nwk", "w") as out:
            out.write(sub.newick() + "\\n")
        res = subprocess.run(["TreeCluster.py", "-i", "subtree.nwk", "-m", "single_linkage", "-t", threshold, "-o", "treecluster.tsv"],
                             capture_output=True, text=True)
        if res.returncode != 0:
            sys.stderr.write(f"PhyloPostTree: TreeCluster failed for subclade {sc}: {res.stderr.strip()[:500]}\\n")
            sys.exit(1)
        with open("treecluster.tsv") as f:
            for row in csv.DictReader(f, delimiter="\\t"):
                if row["ClusterNumber"] != "-1":
                    groups.setdefault(row["ClusterNumber"], []).append(row["SequenceName"])
    # Deterministic labels: largest cluster first, ties broken by the first tip name
    ranked = sorted((sorted(g) for g in groups.values() if len(g) >= 2), key=lambda g: (-len(g), g[0]))
    for k, g in enumerate(ranked, start=1):
        for tip in g:
            cluster_of[tip] = f"{sc}_C{k}"
            size_of[tip] = len(g)
    for tip in members:
        cluster_of.setdefault(tip, "unclustered")

with open("clusters.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["tip", "subclade", "cluster", "cluster_size"])
    for tip in sorted(tips):
        w.writerow([tip, subclade_of.get(tip, "NA"), cluster_of.get(tip, "NA"), size_of.get(tip, "")])

# Amino acid changes on the branch into each cluster's MRCA
proteins = {}
for path in sorted(glob.glob("ancestral_aa/*.fasta")):
    gene = os.path.basename(path)[:-len(".fasta")]
    proteins[gene] = {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}

clusters = {}
for tip, label in cluster_of.items():
    if label not in ("unclustered", "NA"):
        clusters.setdefault(label, []).append(tip)

rows = []
for label in sorted(clusters):
    members = sorted(clusters[label])
    sc = subclade_of[members[0]]
    others = [t for t in by_subclade[sc] if t not in set(members)]
    mrca = full.common_ancestor(members)
    path = full.get_path(mrca)
    if not path:
        continue  # the cluster MRCA is the root: no parent branch
    parent = path[-2] if len(path) >= 2 else full.root
    for gene in sorted(proteins):
        seqs = proteins[gene]
        if mrca.name not in seqs or parent.name not in seqs:
            continue
        for i, (a, b) in enumerate(zip(seqs[parent.name], seqs[mrca.name]), start=1):
            if a == b or a in "X-" or b in "X-":
                continue
            with_c = sum(1 for t in members if seqs[t][i - 1] == b)
            with_o = sum(1 for t in others if seqs[t][i - 1] == b)
            rows.append([label, sc, mrca.name, gene, i, a, b, f"{a}{i}{b}", with_c, len(members), with_o, len(others),
                         "yes" if with_o == 0 else "no"])

with open("cluster_mutations.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["cluster", "subclade", "mrca_node", "gene", "position", "parent_aa", "mrca_aa", "mutation",
                "n_cluster_with", "n_cluster", "n_other_subclade_with", "n_other_subclade", "exclusive_in_subclade"])
    w.writerows(rows)
print(f"PhyloPostTree: {len(clusters)} clusters; {sum(1 for v in cluster_of.values() if v == 'unclustered')} samples unclustered.")
PYEOF
    """
}
