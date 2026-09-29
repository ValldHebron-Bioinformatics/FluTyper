#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloAlignment {
    // Builds the whole-genome alignment: one MAFFT --auto multiple alignment per segment (coordinate reference,
    // root strains and the samples that passed coverage QC), trimmed to the columns of the coordinate reference,
    // then concatenated in --phyloSegmentOrder. Trimming drops insertions relative to the reference so that
    // alignment columns keep reference coordinates; no manual curation is applied. A root strain lacking a segment
    // gets gaps (missing data) for it.
    label 'big_mem'
    errorStrategy 'ignore'
    // Cache on file content: PhyloPreflight reruns when --aa-positions changes, but its reference outputs do not
    cache 'deep'

    input:
    path(sample_dirs, stageAs: "samples/*")  // segment folders of the samples that passed PhyloCoverageQC
    path(root_refs)                          // phylo_references.fasta (>TIP_ID|SEGMENT)
    path(coordinate_refs)                    // reference_segments.fasta (>SEGMENT)
    path(coordinate_table)                   // reference_segments.tsv
    path(qc_table)                           // merged coverage QC table of all candidate samples
    val(segment_order)                       // comma-separated segment names

    output:
    path("whole_genome/phylo_alignment.fasta"), emit: alignment
    path("whole_genome/phylo_segments.tsv"),    emit: segments
    path("alignment.log"),                      emit: log       // step messages, collected into the tree's phylo.log

    script:
    """
#!/usr/bin/env python3
import csv, io, os, subprocess, sys
from Bio import SeqIO

# Step messages go to alignment.log (part of the tree's phylo.log), not the terminal; errors stay on stderr
LOG = open("alignment.log", "w")
def say(msg):
    LOG.write(f"PhyloAlignment: {msg}\\n")
    LOG.flush()

os.makedirs("whole_genome", exist_ok=True)
order = [s.strip() for s in "${segment_order}".split(",") if s.strip()]
threads = "${task.cpus}"
MIN_SAMPLES = 3

qc_samples = set()
with open("${qc_table}") as f:
    for row in csv.DictReader(f, delimiter="\\t"):
        qc_samples.add(row["sample_id"])
samples = sorted(os.listdir("samples")) if os.path.isdir("samples") else []
say(f"{len(samples)} of {len(qc_samples)} candidate samples passed coverage QC and enter the phylogeny "
      f"({len(qc_samples) - len(samples)} dropped, reasons in the coverage QC summary above).")
if len(samples) < MIN_SAMPLES:
    sys.stderr.write(f"PhyloAlignment: only {len(samples)} samples passed coverage QC; at least {MIN_SAMPLES} are needed to build a tree.\\n")
    sys.exit(1)

coord = {r.id: str(r.seq).upper() for r in SeqIO.parse("${coordinate_refs}", "fasta")}
coord_strain = {}
with open("${coordinate_table}") as f:
    for row in csv.DictReader(f, delimiter="\\t"):
        coord_strain[row["segment"]] = row["reference_strain"]

root_order, root_seqs = [], {}
for r in SeqIO.parse("${root_refs}", "fasta"):
    tip, seg = r.id.rsplit("|", 1)
    if tip not in root_order:
        root_order.append(tip)
    root_seqs[(tip, seg)] = str(r.seq).upper()

clash = set(root_order) & set(samples)
if clash:
    sys.stderr.write(f"PhyloAlignment: sample IDs clash with root strain tip IDs: {sorted(clash)}\\n")
    sys.exit(1)

tips = root_order + samples
concat = {t: [] for t in tips}
segment_rows, start = [], 1
COORD_ID = "__COORDINATE_REFERENCE__"

for seg in order:
    records = [(COORD_ID, coord[seg])]
    records += [(t, root_seqs[(t, seg)]) for t in root_order if (t, seg) in root_seqs]
    for s in samples:
        rec = next(SeqIO.parse(f"samples/{s}/segments/{s}_{seg}.fasta", "fasta"))
        records.append((s, str(rec.seq).upper()))
    fasta = "".join(f">{name}\\n{seq}\\n" for name, seq in records)
    # --threadit 0 keeps the iterative refinement single-threaded so the alignment is reproducible
    res = subprocess.run(["mafft", "--auto", "--thread", threads, "--threadit", "0", "--quiet", "-"],
                         input=fasta, capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(f"PhyloAlignment: MAFFT failed for segment {seg}: {res.stderr.strip()[:500]}\\n")
        sys.exit(1)
    aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(io.StringIO(res.stdout), "fasta")}
    keep = [i for i, c in enumerate(aln[COORD_ID]) if c != "-"]
    removed = len(aln[COORD_ID]) - len(keep)
    for t in tips:
        concat[t].append("".join(aln[t][i] for i in keep) if t in aln else "-" * len(keep))
    segment_rows.append([seg, coord_strain.get(seg, ""), start, start + len(keep) - 1, len(keep), removed])
    start += len(keep)

with open("whole_genome/phylo_alignment.fasta", "w") as out:
    for t in tips:
        out.write(f">{t}\\n{''.join(concat[t])}\\n")
with open("whole_genome/phylo_segments.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["segment", "reference_strain", "start", "end", "length", "insertion_columns_removed"])
    w.writerows(segment_rows)
say(f"whole-genome alignment of {len(tips)} sequences x {start - 1} sites.")
    """
}
