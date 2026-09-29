#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloSubclades {
    // Collects the subclade of every tree tip. Samples reuse the per-sample Nextclade results of GenotypingNextclade
    // (column 'subclade'; the row with the best qc.overallScore, as GenotypingResults does). Root strains are not
    // part of the sample set, so their HA segment is run through the same Nextclade dataset here.
    // Tips without a subclade call stay NA.
    errorStrategy 'ignore'
    // Cache on file content: PhyloPreflight reruns when --aa-positions changes, but its reference outputs do not
    cache 'deep'

    input:
    path(nextclade_csvs, stageAs: "nextclade/*")  // nextclade_results_<sample>.csv files
    path(root_refs)                               // phylo_references.fasta (>TIP_ID|SEGMENT)
    path(dataset)                                 // Nextclade dataset folder for the subtype

    output:
    path("subclades.tsv"), emit: subclades
    path("subclades.log"), emit: log   // step messages, collected into every tree's phylo.log

    script:
    """
    : > subclades.log
    python3 - <<'PYEOF'
from Bio import SeqIO
with open("root_ha.fasta", "w") as out:
    for r in SeqIO.parse("${root_refs}", "fasta"):
        tip, seg = r.id.rsplit("|", 1)
        if seg == "HA":
            out.write(f">{tip}\\n{r.seq}\\n")
PYEOF

    if [ -s root_ha.fasta ]; then
        nextclade run --input-dataset "${dataset}" --output-csv root_nextclade.csv root_ha.fasta > nextclade.log 2>&1 \\
            || { cat nextclade.log >&2; exit 1; }
        grep -iE "warn" nextclade.log | sed "s/^/PhyloSubclades: Nextclade /" >> subclades.log || true
    fi

    python3 - <<'PYEOF'
import csv, glob, os, re

def best_row(path):
    best, best_score = None, float("inf")
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter=";"):
            try:
                score = float(row.get("qc.overallScore", ""))
            except ValueError:
                continue
            if score < best_score:
                best, best_score = row, score
    return best

def clean(value):
    value = (value or "").strip()
    return value if value else "NA"

rows, no_call = [], []
for path in sorted(glob.glob("nextclade/nextclade_results_*.csv")):
    sample = re.sub(r"^nextclade_results_(.*)\\.csv\$", r"\\1", os.path.basename(path))
    row = best_row(path) if os.path.getsize(path) > 0 else None
    if row is None:
        no_call.append(sample)
    rows.append([sample, clean(row.get("subclade")) if row else "NA", clean(row.get("clade")) if row else "NA", "sample"])
if no_call:
    with open("subclades.log", "a") as log:
        log.write(f"PhyloSubclades: {len(no_call)} samples have no Nextclade row with a numeric qc.overallScore; subclade set to NA: {', '.join(no_call[:10])}\\n")

if os.path.isfile("root_nextclade.csv"):
    with open("root_nextclade.csv", newline="") as f:
        for row in csv.DictReader(f, delimiter=";"):
            rows.append([row["seqName"], clean(row.get("subclade")), clean(row.get("clade")), "root_strain"])

with open("subclades.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["tip", "subclade", "clade", "source"])
    w.writerows(rows)
PYEOF
    """
}
