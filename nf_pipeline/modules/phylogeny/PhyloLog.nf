#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloLog {
    // Last step of every tree, and the only one that publishes anything besides the HTML report: it assembles the
    // tree's phylo.log and publishes the ML tree as tree.nwk. The phylogenetics processes do not print their
    // per-task chatter to the terminal; each writes its messages to a small part file, and this process puts them
    // together, in order:
    //   1. the subtype-level preflight notes relevant to this tree (preflight.log; notes tied to another segment are left out)
    //   2. a coverage QC summary for this tree, from the merged QC table: samples in, samples dropped and why
    //   3. the subclade step's notes (subclades.log)
    //   4. the tree's own steps: alignment, IQ-TREE, TreeTime/TreeCluster, report
    // A tree whose later step failed (errorStrategy 'ignore') still gets a phylo.log with the steps that ran, so the
    // log shows where it stopped. The parts arrive already grouped by tree_id (groupTuple with remainder), so no
    // per-tree join can drop a tree.
    errorStrategy 'ignore'

    input:
    tuple val(tree_id), val(tree_folder), path(parts)   // alignment.log, iqtree.log, phylo.treefile, posttree.log, report.log (whichever ran)
    path(preflight_log)                                 // <segments>\t<message> notes of this subtype (PhyloPreflight)
    path(subclades_log)
    path(qc_table)                                      // merged coverage QC table of all candidate samples
    val(segment_order)                                  // comma-separated segment names

    output:
    tuple val(tree_id), path("${tree_folder}/phylo.log"), emit: log
    tuple val(tree_id), path("${tree_folder}/tree.nwk"),  optional: true, emit: tree
    // Number of samples in this tree, printed as the only stdout of the script
    tuple val(tree_id), val(tree_folder), stdout, emit: summary

    script:
    """
#!/usr/bin/env python3
import csv, os, re, shutil
from collections import OrderedDict

tree_id, tree_folder = "${tree_id}", "${tree_folder}"
segments = [s.strip() for s in "${segment_order}".split(",") if s.strip()]
whole_genome = tree_folder == "whole_genome"
os.makedirs(tree_folder, exist_ok=True)

def read_lines(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [line.rstrip("\\n") for line in f if line.strip()]

out = []
kind = "whole-genome tree" if whole_genome else f"{tree_folder} segment tree"
out.append(f"# phylo.log: {tree_id} ({kind})")

# 1. Preflight notes: general ones, and those of this tree's segment (every segment for the whole-genome tree)
out.append("")
out.append("## Preflight")
notes = []
for line in read_lines("${preflight_log}"):
    segs, _, msg = line.partition("\\t")
    segs = [s for s in segs.split(",") if s]
    if not segs or whole_genome or tree_folder in segs:
        notes.append(msg)
out.extend(f"PhyloPreflight: {m}" for m in notes)

# 2. Coverage QC summary
with open("${qc_table}") as f:
    rows = list(csv.DictReader(f, delimiter="\\t"))
by_sample = OrderedDict()
for r in rows:
    by_sample.setdefault(r["sample_id"], {})[r["segment"]] = r
def failure(sample):
    segs = by_sample[sample]
    wanted = segments if whole_genome else [tree_folder]
    reasons = []
    for seg in wanted:
        r = segs.get(seg)
        if r is None:
            reasons.append(f"{seg}: no QC row")
        elif r["status"] != "PASS":
            reason = re.sub(r"coverage [\\d.]+% below ([\\d.]+%)", r"coverage below \\1", r["reason"] or "failed")
            reasons.append(reason if not whole_genome else f"{seg} {reason}")
    return "; ".join(reasons)
passing = [s for s in by_sample if not failure(s)]
dropped = OrderedDict()
for s in by_sample:
    why = failure(s)
    if why:
        dropped.setdefault(why, []).append(s)
out.append("")
out.append("## Coverage QC")
out.append(f"PhyloCoverageQC: {len(passing)} of {len(by_sample)} candidate samples enter this tree; "
           f"{len(by_sample) - len(passing)} dropped.")
for why, ids in dropped.items():
    shown = ", ".join(ids[:10]) + (", ..." if len(ids) > 10 else "")
    out.append(f"PhyloCoverageQC:   {len(ids)} dropped ({why}): {shown}")

# 3. Subclades
sub = read_lines("${subclades_log}")
if sub:
    out.append("")
    out.append("## Subclades")
    out.extend(sub)

# 4. The tree's own steps, in pipeline order
for title, part in [("Alignment", "alignment.log"), ("IQ-TREE", "iqtree.log"),
                    ("TreeTime and clusters", "posttree.log"), ("Report", "report.log")]:
    out.append("")
    out.append(f"## {title}")
    lines = read_lines(part)
    if part == "iqtree.log" and os.path.isfile("phylo.treefile"):
        lines.append("PhyloTree: tree written to tree.nwk (ultrafast bootstrap support as node labels).")
    if not lines and not os.path.isfile(part):
        lines = [f"(step did not finish: no {part} was produced)"] if part != "alignment.log" else []
    out.extend(lines)

with open(f"{tree_folder}/phylo.log", "w") as f:
    f.write("\\n".join(out) + "\\n")
if os.path.isfile("phylo.treefile"):
    shutil.copy("phylo.treefile", f"{tree_folder}/tree.nwk")
print(len(passing), end="")
    """
}
