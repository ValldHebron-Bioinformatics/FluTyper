#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloCoverageQC {
    // Measures, for one sample, the coverage of each segment against the coordinate reference: the fraction of
    // reference positions covered by a non-N base in a pairwise MAFFT alignment. A sample enters the whole-genome
    // phylogeny only if every segment reaches --phyloMinCoverage. Missing segments and segments with more than one
    // record (possible coinfection) fail the sample instead of being guessed.
    // D24 (user decision 2026-09-30): a sample whose HA Nextclade analysis FAILED is excluded from every phylogeny of its
    // subtype: all its segments FAIL with reason "HA Nextclade analysis failed: <error>" (so it enters neither the
    // whole-genome nor any segment tree, and the reason reaches every phylo.log and the QC table). Failed = the
    // Nextclade CSV's `errors` column is non-empty: Nextclade fills it only when the sequence could not be analysed
    // (seed matching / alignment failed; clade and QC columns are then empty). A merely "bad" qc.overallStatus, or an
    // aligned sequence without a clade, has an empty `errors` column and is kept. No Nextclade CSV (e.g. no HA segment)
    // or an empty one (Nextclade not run for this H type) is not a failure. Real case: VH_24-26 samples 600497186 and
    // 198309324 (non-pdm09 genomes typed A(H1N1)pdm09, Phase 7 finding P7-1).
    errorStrategy 'ignore'
    // Cache on file content: PhyloPreflight reruns when --aa-positions changes, but its reference outputs do not
    cache 'deep'

    input:
    tuple val(sample_id), path(sample_dir), path(nextclade, stageAs: "nextclade/*")   // HA Nextclade CSV (';'), or []
    path(reference_segments)   // one record per segment, header = segment name
    val(min_coverage)
    val(segments)              // comma-separated segment names

    output:
    tuple val(sample_id), path(sample_dir), path("${sample_id}_phylo_coverage.tsv"), emit: results
    tuple val(sample_id), path("PQerrors.log"), optional: true, emit: errors

    script:
    """
#!/usr/bin/env python3
import csv, io, os, re, subprocess
from Bio import SeqIO

sample_id = "${sample_id}"
min_cov = float("${min_coverage}")
segments = [s.strip() for s in "${segments}".split(",") if s.strip()]
refs = {r.id: str(r.seq).upper() for r in SeqIO.parse("${reference_segments}", "fasta")}

def log(msg):
    with open("PQerrors.log", "a") as f:
        f.write(f"PhyloCoverageQC: {msg}\\n")

def nextclade_failure(path):
    # the first non-empty `errors` cell of the sample's HA Nextclade CSV, without the per-sequence prefix; None = not failed
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    with open(path, newline="", encoding="utf-8") as f:
        errs = [(r.get("errors") or "").strip() for r in csv.DictReader(f, delimiter=";")]
    errs = [e for e in errs if e]
    if not errs:
        return None
    return re.sub(r"^When processing sequence #\\d+ '.*?': ", "", errs[0])[:200]

rows = []
failed = nextclade_failure("${nextclade}".strip())
if failed:
    rows = [[sample_id, seg, "", "", "", "FAIL", f"HA Nextclade analysis failed: {failed}"] for seg in segments]
for seg in ([] if failed else segments):
    if seg not in refs:
        log(f"segment {seg} is missing from the coordinate references")
        rows.append([sample_id, seg, "", "", "", "FAIL", "segment missing from the coordinate references"])
        continue
    ref = refs[seg]
    seg_fasta = f"${sample_dir}/segments/{sample_id}_{seg}.fasta"
    records = list(SeqIO.parse(seg_fasta, "fasta")) if os.path.isfile(seg_fasta) else []
    if not records:
        rows.append([sample_id, seg, len(ref), 0, "0.0000", "FAIL", "segment missing"])
        continue
    if len(records) > 1:
        rows.append([sample_id, seg, len(ref), "", "", "FAIL", f"{len(records)} records for this segment (possible coinfection)"])
        continue
    query = str(records[0].seq).upper()
    res = subprocess.run(["mafft", "--auto", "--quiet", "-"], input=f">ref\\n{ref}\\n>query\\n{query}\\n",
                         capture_output=True, text=True)
    aln = list(SeqIO.parse(io.StringIO(res.stdout), "fasta")) if res.returncode == 0 else []
    if len(aln) != 2:
        log(f"MAFFT failed for sample {sample_id} segment {seg}: {res.stderr.strip()[:300]}")
        rows.append([sample_id, seg, len(ref), "", "", "FAIL", "alignment to the reference failed"])
        continue
    ref_row, query_row = str(aln[0].seq).upper(), str(aln[1].seq).upper()
    covered = sum(1 for r, q in zip(ref_row, query_row) if r != "-" and q not in "-N")
    coverage = covered / len(ref)
    status = "PASS" if coverage >= min_cov else "FAIL"
    reason = "" if status == "PASS" else f"coverage {coverage:.2%} below {min_cov:.2%}"
    rows.append([sample_id, seg, len(ref), covered, f"{coverage:.4f}", status, reason])

with open(f"{sample_id}_phylo_coverage.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["sample_id", "segment", "reference_length", "covered_bases", "coverage", "status", "reason"])
    w.writerows(rows)
    """
}
