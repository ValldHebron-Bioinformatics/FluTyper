#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloPreflight {
    // Startup step of the phylogenetics module; every other phylogenetic step waits for it.
    //  1. Checks that every external tool and R package is installed. A missing (or broken) tool stops the run with
    //     the tool name and the expected version; a different version only warns and is recorded in the manifest.
    //  2. Picks the coordinate reference of each segment: the alignment reference (A/Darwin/6/2021 for H3N2) when the
    //     protocol references contain that segment, otherwise the fallback reference (A/Massachusetts/18/2022).
    //  3. Extracts the segments of every root strain (WHO reference / vaccine strains). A segment missing from the
    //     protocol references is left out and becomes missing data (gaps) in the alignment, never imputed.
    //  4. Maps every codon of every protein onto the coordinate references by aligning the protocol CDS references
    //     against them, so amino acid numbering follows the protocol CDS (HA1 = mature HA1 numbering).
    //  5. Validates --aa-positions against that reference translation and fails fast when a position is out of range.
    errorStrategy 'terminate'

    input:
    val(expected_tools_b64)   // base64 JSON {tool: expected_version}, from params.phyloTools
    path(ref_genomes)         // protocol whole-genome references (HUMAN_references.fasta)
    path(cds_refs)            // protocol CDS references (CDS_references.fasta)
    val(subtype)              // e.g. H3N2
    val(reference)            // alignment reference strain, e.g. A/Darwin/6/2021
    val(reference_fallback)   // coordinate reference for segments the alignment reference lacks
    val(root_strains)         // comma-separated strains added to the tree for rooting
    val(segments)             // comma-separated segment names
    val(aa_positions)         // normalised comma-separated GENE:POS list, may be empty

    output:
    path("phylo_tool_versions.json"), emit: versions
    path("phylo_references.fasta"),   emit: references
    path("reference_segments.fasta"), emit: coordinate_refs
    path("reference_segments.tsv"),   emit: coordinate_table
    path("root_strains.tsv"),         emit: root_table
    path("gene_map.tsv"),             emit: gene_map
    path("aa_positions.tsv"),         emit: positions
    path("preflight.log"),            emit: log        // notes for the per-tree phylo.log, see say() below

    script:
    """
#!/usr/bin/env python3
import base64, csv, io, json, re, shutil, subprocess, sys
from Bio import SeqIO
from Bio.Seq import Seq

def fail(msg):
    sys.stderr.write(f"PhyloPreflight: {msg}\\n")
    sys.exit(1)

# Notes are not printed to the terminal: they are collected in preflight.log for each tree's phylo.log. One line per
# note, "<segments>\\t<message>": segments (comma-separated) restricts the note to the trees of those segments (the
# whole-genome tree includes every segment); empty means it concerns every tree of the subtype.
open("preflight.log", "w").close()
def say(msg, segments=""):
    with open("preflight.log", "a") as f:
        f.write(f"{segments}\\t{msg}\\n")

def write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
        w.writerow(header)
        w.writerows(rows)

# ---- 1. Tools ------------------------------------------------------------------------------------------------------
expected = json.loads(base64.b64decode("${expected_tools_b64}").decode("utf-8"))

# How each tool reports its version. Anything not listed is called as '<tool> --version'.
VERSION_CMDS = {
    "mafft":          ["mafft", "--version"],
    "iqtree2":        ["iqtree2", "--version"],
    "treetime":       ["treetime", "--version"],
    "TreeCluster.py": ["TreeCluster.py", "--version"],
    "nextclade":      ["nextclade", "--version"],
    "Rscript":        ["Rscript", "--version"],
}
R_PACKAGES = {"ggtree", "treeio", "ape", "ggplot2", "ggnewscale", "svglite"}

def first_version(text):
    m = re.search(r"(\\d+(?:\\.\\d+)+)", text or "")
    return m.group(1) if m else None

def r_package_version(pkg):
    if shutil.which("Rscript") is None:
        return None
    res = subprocess.run(["Rscript", "-e", f"cat(as.character(packageVersion('{pkg}')))"], capture_output=True, text=True)
    return first_version(res.stdout) if res.returncode == 0 else None

report, missing = {}, []
for tool, want in expected.items():
    if tool in R_PACKAGES:
        found, kind = r_package_version(tool), "R package"
    else:
        cmd, kind = VERSION_CMDS.get(tool, [tool, "--version"]), "binary"
        if shutil.which(cmd[0]) is None:
            found = None
        else:
            res = subprocess.run(cmd, capture_output=True, text=True)
            # A tool that is on PATH but cannot run (broken install) counts as missing
            found = (first_version(res.stdout + res.stderr) or "unknown") if res.returncode == 0 else None
    if found is None:
        missing.append(f"{kind} '{tool}' (expected version {want})")
        status = "missing"
    elif found != str(want):
        say(f"WARNING {kind} '{tool}' is version {found}, expected {want}. Results may differ from the pinned environment.")
        status = "version_mismatch"
    else:
        status = "ok"
    report[tool] = {"expected": str(want), "found": found, "status": status}

with open("phylo_tool_versions.json", "w") as f:
    json.dump(report, f, indent=2, sort_keys=True)
if missing:
    fail("the phylogenetics module needs tools that are not installed:\\n  - " + "\\n  - ".join(missing)
         + "\\nInstall them with: conda env update -f FluTyper_env.yaml (or run without --phylogenetics).")
say("all phylogenetics tools found.")

# ---- 2. Coordinate reference of each segment -----------------------------------------------------------------------
subtype = "${subtype}"
reference = "${reference}"
fallback = "${reference_fallback}"
root_strains = [s.strip() for s in "${root_strains}".split(",") if s.strip()]
segments = [s.strip() for s in "${segments}".split(",") if s.strip()]
requested = [p.strip() for p in "${aa_positions}".split(",") if p.strip()]
ref_subtype = "H1N1" if subtype == "H1N1pdm09" else subtype

PROT_BY_SEGMENT = {
    "HA": ["HA1", "HA2"], "NA": ["NA"], "PB2": ["PB2"], "PB1": ["PB1", "PB1-F2"],
    "PA": ["PA", "PA-X"], "NP": ["NP"], "MP": ["M1", "M2"], "NS": ["NS1", "NS2"],
}

def ref_id(strain):
    return "REF_" + re.sub(r"[^A-Za-z0-9.-]", "_", strain)

# Protocol headers are {SUBTYPE}_{SEGMENT}_{STRAIN}[_EPI_ISL_...]
genomes = list(SeqIO.parse("${ref_genomes}", "fasta"))

def find_segment(strain, seg):
    prefix = f"{ref_subtype}_{seg}_{strain.replace(' ', '_')}"
    hits = [r for r in genomes if r.description == prefix or r.description.startswith(prefix + "_")]
    if len(hits) > 1:
        say(f"WARNING {len(hits)} records match {prefix}; using the first one.", seg)
    return str(hits[0].seq).upper().replace("U", "T") if hits else None

coord = {}
for seg in segments:
    for strain in (reference, fallback):
        seq = find_segment(strain, seg) if strain else None
        if seq:
            coord[seg] = (strain, seq)
            break
    else:
        fail(f"segment {seg} is missing for both the reference {reference} and the fallback {fallback} in ${ref_genomes}.")
    if coord[seg][0] != reference:
        say(f"{reference} has no {seg} segment in the protocol references; {seg} coordinates use {coord[seg][0]}.", seg)

with open("reference_segments.fasta", "w") as out:
    for seg in segments:
        out.write(f">{seg}\\n{coord[seg][1]}\\n")
write_tsv("reference_segments.tsv", ["segment", "reference_strain", "length"],
          [[seg, coord[seg][0], len(coord[seg][1])] for seg in segments])

# ---- 3. Root strain segments (tips of the tree) --------------------------------------------------------------------
root_rows = []
with open("phylo_references.fasta", "w") as out:
    for strain in root_strains:
        present, absent = [], []
        for seg in segments:
            seq = find_segment(strain, seg)
            if seq:
                out.write(f">{ref_id(strain)}|{seg}\\n{seq}\\n")
                present.append(seg)
            else:
                absent.append(seg)
        if not present:
            fail(f"root strain {strain} ({subtype}) is not in ${ref_genomes}.")
        if absent:
            say(f"root strain {strain} has no {','.join(absent)} segment(s); coded as missing data in the alignment.", ",".join(absent))
        root_rows.append([strain, ref_id(strain), ",".join(present), ",".join(absent)])
write_tsv("root_strains.tsv", ["strain", "tip_id", "segments_present", "segments_missing"], root_rows)

# ---- 4. Codon map of every protein onto the coordinate references --------------------------------------------------
cds = {}
for r in SeqIO.parse("${cds_refs}", "fasta"):
    parts = r.id.split("_")
    if len(parts) > 2 and parts[0] in {subtype, ref_subtype}:
        cds.setdefault(parts[1], str(r.seq).upper().replace("U", "T"))

def align_pair(seg_seq, cds_seq):
    # E-INS-i handles the long gaps of spliced proteins (M2, NS2, PA-X) against the full segment
    res = subprocess.run(["mafft", "--genafpair", "--maxiterate", "1000", "--quiet", "-"],
                         input=f">segment\\n{seg_seq}\\n>cds\\n{cds_seq}\\n", capture_output=True, text=True, check=True)
    recs = list(SeqIO.parse(io.StringIO(res.stdout), "fasta"))
    return str(recs[0].seq).upper(), str(recs[1].seq).upper()

rows, gene_len = [], {}
for seg in segments:
    seg_strain, seg_seq = coord[seg]
    for gene in PROT_BY_SEGMENT.get(seg, []):
        if gene not in cds:
            continue
        aln_seg, aln_cds = align_pair(seg_seq, cds[gene])
        cds_to_seg, seg_pos, cds_pos = {}, 0, 0
        for s, c in zip(aln_seg, aln_cds):
            if s != "-":
                seg_pos += 1
            if c != "-":
                cds_pos += 1
                if s != "-":
                    cds_to_seg[cds_pos] = seg_pos
        n_codons = len(cds[gene]) // 3
        gene_rows, matches = [], 0
        for k in range(n_codons):
            nts = [cds_to_seg.get(3 * k + i + 1) for i in range(3)]
            if None in nts:
                codon, aa = "NNN", "?"
            else:
                codon = "".join(seg_seq[n - 1] for n in nts)
                aa = str(Seq(codon).translate())
            gene_rows.append([gene, k + 1, seg] + [n if n else "" for n in nts] + [codon, aa])
            if aa == str(Seq(cds[gene][3 * k:3 * k + 3]).translate()):
                matches += 1
        # Terminal stop codon is not part of the protein numbering
        if gene_rows and gene_rows[-1][-1] == "*":
            gene_rows.pop()
        if n_codons and matches / n_codons < 0.9:
            say(f"WARNING {gene} of {seg_strain} matches only {matches / n_codons:.1%} of the protocol CDS reference; check the reference sequences.", seg)
        unmapped = sum(1 for r in gene_rows if r[-1] == "?")
        if unmapped:
            say(f"WARNING {unmapped} codons of {gene} could not be placed on {seg_strain} segment {seg}.", seg)
        rows.extend(gene_rows)
        gene_len[gene] = len(gene_rows)

write_tsv("gene_map.tsv", ["gene", "aa_position", "segment", "nt1", "nt2", "nt3", "ref_codon", "ref_aa"], rows)

# ---- 5. Validate requested positions against the reference translation --------------------------------------------
ref_aa = {(r[0], r[1]): r[-1] for r in rows}
errors, out_rows = [], []
for i, label in enumerate(requested, start=1):
    gene, pos = label.split(":")
    pos = int(pos)
    if gene not in gene_len:
        errors.append(f"{label}: gene {gene} is not in the {subtype} reference. Available genes: {', '.join(sorted(gene_len))}.")
    elif not 1 <= pos <= gene_len[gene]:
        errors.append(f"{label} is out of range: valid positions for {gene} are 1-{gene_len[gene]}.")
    elif ref_aa[(gene, pos)] == "?":
        errors.append(f"{label}: this codon could not be placed on the reference.")
    else:
        out_rows.append([i, label, gene, pos, ref_aa[(gene, pos)]])
if errors:
    fail("invalid --aa-positions:\\n  " + "\\n  ".join(errors))

write_tsv("aa_positions.tsv", ["order", "label", "gene", "position", "ref_aa"], out_rows)
    """
}
