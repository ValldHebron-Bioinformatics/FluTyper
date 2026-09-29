#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloSegmentAlignment {
    // One segment's alignment (--phyloSegmentTrees): MAFFT --auto of the samples that passed PhyloCoverageQC for
    // this segment (a different, usually larger, set than the whole-genome candidates) plus the root strains that
    // have this segment. A root strain missing the segment is dropped as a tip (not gapped: there is nothing else
    // to place it by, unlike the whole-genome alignment). If the subtype's configured --phyloOutgroup lacks the
    // segment, the oldest root strain (by the trailing year in its name) that has it becomes this tree's outgroup
    // instead; that fallback is recorded in <segment>_root_used.tsv for the run manifest.
    label 'big_mem'
    errorStrategy 'ignore'
    // Cache on file content: PhyloPreflight reruns when --aa-positions changes, but its reference outputs do not
    cache 'deep'
    debug true

    input:
    tuple val(segment), path(sample_dirs, stageAs: "samples/*")  // segment folders of the samples that passed PhyloCoverageQC for this segment
    path(root_refs)          // phylo_references.fasta (>TIP_ID|SEGMENT)
    path(coordinate_refs)    // reference_segments.fasta (>SEGMENT)
    path(coordinate_table)   // reference_segments.tsv
    val(outgroup_strain)     // the subtype's configured --phyloOutgroup (may lack this segment)
    val(root_strains)        // the subtype's comma-separated root strains (candidates for the fallback outgroup)

    output:
    tuple val(segment), path("${segment}/${segment}_alignment.fasta"), emit: alignment
    tuple val(segment), path("${segment}/${segment}_segments.tsv"),    emit: segments
    tuple val(segment), path("${segment}/${segment}_outgroup.txt"),    emit: outgroup
    tuple val(segment), path("${segment}/${segment}_root_used.tsv"),   emit: root_used

    script:
    """
#!/usr/bin/env python3
import csv, io, os, re, subprocess, sys
from Bio import SeqIO

segment = "${segment}"
threads = "${task.cpus}"
os.makedirs(segment, exist_ok=True)

coord = {r.id: str(r.seq).upper() for r in SeqIO.parse("${coordinate_refs}", "fasta")}
coord_strain = {}
with open("${coordinate_table}") as f:
    for row in csv.DictReader(f, delimiter="\\t"):
        coord_strain[row["segment"]] = row["reference_strain"]

# Root strain segment sequences: only strains that HAVE this segment enter a segment tree (no gapping)
root_seqs = {}
for r in SeqIO.parse("${root_refs}", "fasta"):
    tip, seg = r.id.rsplit("|", 1)
    if seg == segment:
        root_seqs[tip] = str(r.seq).upper()

def ref_id(strain):
    return "REF_" + re.sub(r"[^A-Za-z0-9.-]", "_", strain)

def strain_year(strain):
    m = re.search(r"(\\d{4})\$", strain)
    return int(m.group(1)) if m else 9999

configured_outgroup_strain = "${outgroup_strain}"
root_strain_list = [s.strip() for s in "${root_strains}".split(",") if s.strip()]
configured_outgroup_tip = ref_id(configured_outgroup_strain)

fallback_used = False
if configured_outgroup_tip in root_seqs:
    outgroup_strain_used, outgroup_tip = configured_outgroup_strain, configured_outgroup_tip
else:
    candidates = [s for s in root_strain_list if ref_id(s) in root_seqs]
    if not candidates:
        sys.stderr.write(f"PhyloSegmentAlignment: no root strain has segment {segment}; cannot root the {segment} tree.\\n")
        sys.exit(1)
    outgroup_strain_used = min(candidates, key=strain_year)
    outgroup_tip = ref_id(outgroup_strain_used)
    fallback_used = True
    print(f"PhyloSegmentAlignment: {configured_outgroup_strain} has no {segment} segment; rooting the {segment} tree on {outgroup_strain_used} instead.")

with open(f"{segment}/{segment}_outgroup.txt", "w") as f:
    f.write(outgroup_tip + "\\n")
with open(f"{segment}/{segment}_root_used.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["segment", "configured_outgroup", "outgroup_strain_used", "outgroup_tip", "fallback_used"])
    w.writerow([segment, configured_outgroup_strain, outgroup_strain_used, outgroup_tip, "yes" if fallback_used else "no"])

samples = sorted(os.listdir("samples")) if os.path.isdir("samples") else []
root_order = sorted(root_seqs.keys())
clash = set(root_order) & set(samples)
if clash:
    sys.stderr.write(f"PhyloSegmentAlignment: sample IDs clash with root strain tip IDs: {sorted(clash)}\\n")
    sys.exit(1)
tips = root_order + samples

COORD_ID = "__COORDINATE_REFERENCE__"
records = [(COORD_ID, coord[segment])]
records += [(t, root_seqs[t]) for t in root_order]
for s in samples:
    rec = next(SeqIO.parse(f"samples/{s}/segments/{s}_{segment}.fasta", "fasta"))
    records.append((s, str(rec.seq).upper()))
fasta = "".join(f">{name}\\n{seq}\\n" for name, seq in records)
# --threadit 0 keeps the iterative refinement single-threaded so the alignment is reproducible
res = subprocess.run(["mafft", "--auto", "--thread", threads, "--threadit", "0", "--quiet", "-"],
                     input=fasta, capture_output=True, text=True)
if res.returncode != 0:
    sys.stderr.write(f"PhyloSegmentAlignment: MAFFT failed for segment {segment}: {res.stderr.strip()[:500]}\\n")
    sys.exit(1)
aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(io.StringIO(res.stdout), "fasta")}
keep = [i for i, c in enumerate(aln[COORD_ID]) if c != "-"]
removed = len(aln[COORD_ID]) - len(keep)

with open(f"{segment}/{segment}_alignment.fasta", "w") as out:
    for t in tips:
        seq = "".join(aln[t][i] for i in keep) if t in aln else "-" * len(keep)
        out.write(f">{t}\\n{seq}\\n")
with open(f"{segment}/{segment}_segments.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["segment", "reference_strain", "start", "end", "length", "insertion_columns_removed"])
    w.writerow([segment, coord_strain.get(segment, ""), 1, len(keep), len(keep), removed])
print(f"PhyloSegmentAlignment: {segment} alignment of {len(tips)} sequences x {len(keep)} sites "
      f"(outgroup {outgroup_tip}{' [fallback]' if fallback_used else ''}).")
    """
}
