#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloReport {
    // Final step of the phylogenetics module, the only one that reruns when figure options change.
    //  1. annotations.tsv, one tidy row per tree tip: id, tip type, subclade, cluster, source group, period and the
    //     amino acid state at every requested position (aa_states.tsv also lists internal nodes). The source group
    //     comes from ORIGINATING_LAB, collapsed through --labMap when given (unmapped labs -> Other); the period comes
    //     from DATE binned by --periodBin; season (from ISO week 40, as in the other reports), AGE GROUP and SEX are
    //     added as well. Metadata is never imputed: missing or unparseable values stay NA, a
    //     year-month date gets a month but no week, and root strains have no metadata.
    //  2. The annotated tree (Koumaty et al. 2026, Fig. 2) with R/ggtree: rectangular ML tree, ultrafast bootstrap
    //     support, scale bar in substitutions per site, cluster labels and one heatmap column per panel. Panels come
    //     from the declarative spec (RESOURCES/phylo_panels.tsv) filtered by availability: a non-required panel whose
    //     metadata field is absent, or non-missing for fewer than --minPanelCoverage of the samples, is omitted with an
    //     INFO line and the next panels move left. Residue colours are fixed per amino acid letter.
    //  3. PhylogeneticTreeReport_<subtype_label>.html: interactive version of the tree for the index.html dashboard,
    //     with dropdowns for season (All Time, or from one season to another), age group and sex. A selection prunes
    //     the tree to the matching samples (reference strains always stay), keeping the branch lengths of the full
    //     ML tree.
    //  4. run_manifest.json: tool versions, resolved parameters, input checksums, method choices, coverage QC,
    //     clusters, and the panels rendered or skipped with the reason.
    errorStrategy 'ignore'
    debug true

    input:
    // subtype_label names and tags every tree of the module: "H3N2"/"H1N1pdm09" for the (optional) whole-genome
    // tree, "H3N2_HA" etc. for a segment tree; used for the report filename/title. tree_folder is where this
    // tree's files live under the subtype's own output folder: "whole_genome" or the bare segment name e.g. "HA".
    // tree/support/clusters/proteins/segments are joined by subtype_label before this call, since they come from
    // separate per-tree PhyloTree/PhyloPostTree/PhyloAlignment(-Segment) invocations fanned out over one channel
    // (see PHYLOGENETICS in subworkflows/Phylogenetics.nf).
    tuple val(subtype_label), val(tree_folder), path(tree), path(support), path(clusters), path(proteins), path(segments)
    path(subclades)                             // subclades.tsv (shared: subclade assignment does not depend on the tree)
    path(positions)                             // aa_positions.tsv (shared; positions whose gene has no ancestral_aa/<gene>.fasta in `proteins` are dropped, so a segment tree only shows its own genes)
    path(metadata)                              // metadata CSV, or empty when --metadata is not given
    path(lab_map)                               // --labMap TSV, or empty
    path(panel_spec)                            // phylo_panels.tsv
    path(versions)                              // phylo_tool_versions.json
    path(qc_table)                              // phylo_coverage_qc.tsv
    path(coordinate_table)                      // reference_segments.tsv
    path(root_table)                            // root_strains.tsv
    path(checksum_inputs, stageAs: "inputs/*")  // input files to checksum
    val(period_bin)                             // month, week or none
    val(min_panel_coverage)
    val(width)
    val(height)
    val(dpi)
    val(min_support_label)                      // minimum ultrafast bootstrap support (%) drawn as a node label; 0 shows all
    val(run_info_b64)                           // base64 JSON: resolved parameters and run information

    output:
    tuple val(subtype_label), path("${tree_folder}/annotations.tsv"),        emit: annotations
    tuple val(subtype_label), path("${tree_folder}/annotation_sources.tsv"), emit: sources
    tuple val(subtype_label), path("${tree_folder}/aa_states.tsv"),          emit: aa_states
    tuple val(subtype_label), path("${tree_folder}/phylo_tree.svg"),         emit: svg
    tuple val(subtype_label), path("${tree_folder}/phylo_tree.png"),         emit: png
    tuple val(subtype_label), path("${tree_folder}/panels_status.tsv"),      emit: panels
    tuple val(subtype_label), path("${tree_folder}/run_manifest.json"),      emit: manifest
    tuple val(subtype_label), path("${tree_folder}/PhylogeneticTreeReport_${subtype_label}.html"), emit: html

    script:
    """
    # ---- 1. Annotations ----
    python3 - <<'PYEOF'
import csv, os, re, sys
import pandas as pd
from Bio import Phylo, SeqIO

meta_path = "${metadata}".strip()
lab_map_path = "${lab_map}".strip()
period_bin = "${period_bin}"
LAB_FIELD, DATE_FIELD = "ORIGINATING_LAB", "DATE"

def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\\t"))

tips = sorted(t.name for t in Phylo.read("${tree}", "newick").get_terminals())
sub = {r["tip"]: r for r in read_tsv("${subclades}")}
clu = {r["tip"]: r["cluster"] for r in read_tsv("${clusters}")}
# aa_positions.tsv is resolved once per subtype and shared by every tree (whole genome + every segment tree); a
# segment tree's ancestral_aa/ only has the genes of its own segment(s) (PhyloPostTree), so positions on any other
# gene are silently dropped here -- a segment tree ends up with only the aa columns that belong to it.
positions = [p for p in read_tsv("${positions}") if os.path.isfile(f"${proteins}/{p['gene']}.fasta")]
proteins = {}
for gene in sorted({p["gene"] for p in positions}):
    proteins[gene] = {r.id: str(r.seq) for r in SeqIO.parse(f"${proteins}/{gene}.fasta", "fasta")}
tip_set = set(tips)
aa = {}
with open("aa_states.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["node", "node_type", "label", "gene", "position", "aa"])
    nodes = sorted(next(iter(proteins.values())), key=lambda n: (n not in tip_set, n)) if proteins else []
    for n in nodes:
        for p in positions:
            state = proteins[p["gene"]][n][int(p["position"]) - 1]
            w.writerow([n, "tip" if n in tip_set else "internal", p["label"], p["gene"], p["position"], state])
            if n in tip_set:
                # X (ambiguous codon) and '-' (gap: missing data or deletion) cannot be told apart from missing: NA
                aa[(n, p["label"])] = "NA" if state in ("X", "-", "") else state

labels = [p["label"] for p in positions]
sources = []

# --- Metadata -------------------------------------------------------------------------------
meta, meta_columns = {}, []
if meta_path:
    df = pd.read_csv(meta_path, dtype=str, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    meta_columns = list(df.columns)
    id_col = next((c for c in df.columns if c.lower().replace("_", "").replace(" ", "") in ("id", "sampleid")), None)
    if id_col is None:
        sys.stderr.write(f"PhyloReport: no sample ID column ('ID', 'Sample_ID' or 'SampleID') in the metadata. Columns: {meta_columns}\\n")
        sys.exit(1)
    blank_ids = int(df[id_col].isna().sum())
    if blank_ids:
        print(f"PhyloReport: {blank_ids} metadata rows have an empty {id_col} and are ignored.")
    df = df[df[id_col].notna()].copy()
    df[id_col] = df[id_col].astype(str).str.strip()
    # One metadata row per sample, keeping the first as MetadataMerge does
    df = df.drop_duplicates(subset=[id_col], keep="first")
    for _, row in df.iterrows():
        meta[row[id_col]] = row
lab_present, date_present = LAB_FIELD in meta_columns, DATE_FIELD in meta_columns

def meta_value(tip, field):
    row = meta.get(tip)
    if row is None or field not in row or pd.isna(row[field]):
        return None
    value = str(row[field]).strip()
    return value or None

# --- Source group ---------------------------------------------------------------------------
lab_map = {}
if lab_map_path:
    with open(lab_map_path) as f:
        for n, line in enumerate(f, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.rstrip("\\n").split("\\t")]
            if len(parts) < 2 or not parts[0] or not parts[1]:
                sys.stderr.write(f"PhyloReport: --labMap line {n} must be ORIGINATING_LAB<TAB>display_group: {line.strip()}\\n")
                sys.exit(1)
            key = parts[0].casefold()
            if n == 1 and key == LAB_FIELD.casefold():
                continue  # header
            if key in lab_map and lab_map[key] != parts[1]:
                sys.stderr.write(f"PhyloReport: --labMap maps '{parts[0]}' to both '{lab_map[key]}' and '{parts[1]}'.\\n")
                sys.exit(1)
            lab_map[key] = parts[1]

def source_group(tip):
    lab = meta_value(tip, LAB_FIELD)
    if lab is None:
        return "NA"
    if lab_map_path:
        return lab_map.get(lab.casefold(), "Other")
    return lab

if not meta_path:
    sources.append(["source_group", LAB_FIELD, "no", "no --metadata provided"])
elif not lab_present:
    note = f"metadata has no {LAB_FIELD} column"
    near = [c for c in meta_columns if c.upper().replace(" ", "_") == LAB_FIELD]
    if near:
        note += f" (found '{near[0]}'; the pipeline reads {LAB_FIELD})"
    sources.append(["source_group", LAB_FIELD, "no", note])
else:
    if lab_map_path:
        labs = [meta_value(t, LAB_FIELD) for t in tips if t in meta]
        labs = [l for l in labs if l is not None]
        matched = sum(1 for l in labs if l.casefold() in lab_map)
        print(f"PhyloReport: --labMap matched {matched} of {len(labs)} originating labs; the rest are grouped as Other.")
        sources.append(["source_group", LAB_FIELD, "yes", f"grouped with --labMap ({matched}/{len(labs)} labs matched, rest Other)"])
    else:
        sources.append(["source_group", LAB_FIELD, "yes", "raw values"])

# --- Period ---------------------------------------------------------------------------------
unparsed = []

def period(tip):
    if period_bin == "none":
        return "NA"
    raw = meta_value(tip, DATE_FIELD)
    if raw is None:
        return "NA"
    if re.fullmatch(r"\\d{4}-\\d{2}", raw):  # year-month only: no day to place it in a week
        return raw if period_bin == "month" else "NA"
    if re.fullmatch(r"\\d{4}", raw):
        return "NA"
    date = pd.to_datetime(raw, errors="coerce")
    if pd.isna(date):
        unparsed.append(f"{tip}={raw}")
        return "NA"
    if period_bin == "month":
        return f"{date.year:04d}-{date.month:02d}"
    iso = date.isocalendar()
    return f"{iso[0]:04d}-W{iso[1]:02d}"

if period_bin == "none":
    sources.append(["period", DATE_FIELD, "no", "--periodBin none"])
elif not meta_path:
    sources.append(["period", DATE_FIELD, "no", "no --metadata provided"])
elif not date_present:
    sources.append(["period", DATE_FIELD, "no", f"metadata has no {DATE_FIELD} column"])
else:
    sources.append(["period", DATE_FIELD, "yes", f"binned by {period_bin}"])

# --- Season, age group, sex -----------------------------------------------------------------
# Season as in the other FluTyper reports: it starts in ISO week 40 ("2024-2025" = 2024-W40 to 2025-W39).
# A year-month date gets a season only when the whole month falls in a single season.
AGE_FIELD, SEX_FIELD = "AGE GROUP", "SEX"

def season_of(date):
    iso = date.isocalendar()
    start = iso[0] if iso[1] >= 40 else iso[0] - 1
    return f"{start}-{start + 1}"

def season(tip):
    raw = meta_value(tip, DATE_FIELD)
    if raw is None or re.fullmatch(r"\\d{4}", raw):
        return "NA"
    if re.fullmatch(r"\\d{4}-\\d{2}", raw):
        first = pd.Timestamp(f"{raw}-01")
        seasons = {season_of(first), season_of(first + pd.offsets.MonthEnd(0))}
        return seasons.pop() if len(seasons) == 1 else "NA"
    date = pd.to_datetime(raw, errors="coerce")
    return "NA" if pd.isna(date) else season_of(date)

def plain(tip, field):
    value = meta_value(tip, field)
    return "NA" if value is None else value

for column, field in (("season", DATE_FIELD), ("age_group", AGE_FIELD), ("sex", SEX_FIELD)):
    if not meta_path:
        sources.append([column, field, "no", "no --metadata provided"])
    elif field not in meta_columns:
        sources.append([column, field, "no", f"metadata has no {field} column"])
    else:
        sources.append([column, field, "yes", "season starts in ISO week 40" if column == "season" else "raw values"])

# --- Rows -----------------------------------------------------------------------------------
BASE = ["id", "tip_type", "subclade", "cluster", "source_group", "period", "season", "age_group", "sex"]
rows = []
for tip in tips:
    tip_type = sub.get(tip, {}).get("source", "sample")
    is_sample = tip_type == "sample"
    rows.append([tip, tip_type, sub.get(tip, {}).get("subclade", "NA"), clu.get(tip, "NA"),
                 source_group(tip) if is_sample else "NA",
                 period(tip) if is_sample else "NA",
                 season(tip) if is_sample else "NA",
                 plain(tip, AGE_FIELD) if is_sample else "NA",
                 plain(tip, SEX_FIELD) if is_sample else "NA"]
                + [aa.get((tip, l), "NA") for l in labels])

if unparsed:
    print(f"PhyloReport: {len(unparsed)} DATE values could not be parsed and are left NA: {', '.join(unparsed[:10])}")
missing_meta = [r[0] for r in rows if r[1] == "sample" and meta_path and r[0] not in meta]
if missing_meta:
    print(f"PhyloReport: {len(missing_meta)} samples have no metadata row: {', '.join(missing_meta[:10])}")

for i, label in enumerate(labels):
    col = len(BASE) + i
    states = {r[col] for r in rows if r[1] == "sample" and r[col] != "NA"}
    if len(states) <= 1:
        state = next(iter(states)) if states else "NA"
        print(f"PhyloReport: WARNING position {label} is invariant across the samples ({state}); rendering it anyway.")
        sources.append([label, "", "yes", f"invariant across samples ({state})"])
    else:
        sources.append([label, "", "yes", ""])

with open("annotations.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(BASE + labels)
    w.writerows(rows)
with open("annotation_sources.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
    w.writerow(["column", "metadata_field", "available", "note"])
    w.writerows(sources)
PYEOF

    # ---- 2. Figure ----
    Rscript - <<'REOF'
suppressPackageStartupMessages({
    library(ggtree)
    library(ape)
    library(ggplot2)
    library(ggnewscale)
})

min_cov <- as.numeric("${min_panel_coverage}")
fig_width <- as.numeric("${width}")
fig_height <- as.numeric("${height}")
fig_dpi <- as.numeric("${dpi}")
min_support_label <- as.numeric("${min_support_label}")

NA_COLOUR <- "#D9D9D9"
BASE_COLS <- c("id", "tip_type", "subclade", "cluster", "source_group", "period", "season", "age_group", "sex")
# Fixed colour per amino acid letter (Kelly's contrast colours), identical for every position and run
RESIDUE_COLOURS <- c(A = "#F3C300", C = "#875692", D = "#F38400", E = "#A1CAF1", F = "#BE0032", G = "#C2B280",
    H = "#008856", I = "#E68FAC", K = "#0067A5", L = "#F99379", M = "#604E97", N = "#F6A600", P = "#B3446C",
    Q = "#DCD300", R = "#882D17", S = "#8DB600", T = "#654522", V = "#E25822", W = "#2B3D26", Y = "#00B7C7",
    "*" = "#222222")
CATEGORICAL <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

read_tsv <- function(path) {
    read.delim(path, colClasses = "character", check.names = FALSE, na.strings = character(0), comment.char = "")
}
tree <- read.tree("${tree}")
support <- read_tsv("${support}")
ann <- read_tsv("annotations.tsv")
sources <- read_tsv("annotation_sources.tsv")
spec <- read_tsv("${panel_spec}")
aa_cols <- setdiff(names(ann), BASE_COLS)

# Expand the declarative spec: '@aa_positions' becomes one panel per requested position, in the order given
panels <- list()
for (i in seq_len(nrow(spec))) {
    if (spec[["column"]][i] == "@aa_positions") {
        for (lab in aa_cols) {
            panels[[length(panels) + 1]] <- list(name = lab, column = lab, scale = spec[["scale"]][i], required = spec[["required"]][i])
        }
    } else {
        panels[[length(panels) + 1]] <- list(name = spec[["panel"]][i], column = spec[["column"]][i],
                                             scale = spec[["scale"]][i], required = spec[["required"]][i])
    }
}

# Keep the panels whose data is available
samples <- ann[ann[["tip_type"]] == "sample", , drop = FALSE]
status <- list()
render <- list()
for (p in panels) {
    src <- sources[sources[["column"]] == p[["column"]], , drop = FALSE]
    field <- if (nrow(src) && nzchar(src[["metadata_field"]][1])) src[["metadata_field"]][1] else p[["column"]]
    has_col <- p[["column"]] %in% names(ann)
    coverage <- if (has_col && nrow(samples)) mean(samples[[p[["column"]]]] != "NA") else 0
    reason <- ""
    if (!has_col) {
        reason <- sprintf("column %s is not in annotations.tsv", p[["column"]])
    } else if (nrow(src) && src[["available"]][1] == "no") {
        reason <- src[["note"]][1]
    } else if (p[["required"]] != "yes" && coverage < min_cov) {
        reason <- sprintf("%s is non-missing for %.1f%% of samples, below --minPanelCoverage %.1f%%", field, 100 * coverage, 100 * min_cov)
    }
    if (nzchar(reason)) {
        cat(sprintf("PhyloReport: INFO panel '%s' omitted: %s (field %s, coverage %.1f%%).\\n", p[["name"]], reason, field, 100 * coverage))
        status[[length(status) + 1]] <- c(p[["name"]], p[["column"]], "skipped", sprintf("%.4f", coverage), reason)
    } else {
        note <- if (nrow(src)) src[["note"]][1] else ""
        status[[length(status) + 1]] <- c(p[["name"]], p[["column"]], "rendered", sprintf("%.4f", coverage), note)
        render[[length(render) + 1]] <- p
    }
}

level_order <- function(values) {
    lv <- sort(unique(values[values != "NA"]))
    for (last in c("Other", "unclustered")) {
        if (last %in% lv) lv <- c(setdiff(lv, last), last)
    }
    lv
}
palette_for <- function(scale, levels) {
    n <- length(levels)
    if (n == 0) return(character(0))
    if (scale == "residue") {
        cols <- unname(RESIDUE_COLOURS[levels])
        if (any(is.na(cols))) {
            cat(sprintf("PhyloReport: WARNING residue values without a fixed colour are drawn white: %s\\n", paste(levels[is.na(cols)], collapse = ", ")))
            cols[is.na(cols)] <- "#FFFFFF"
        }
    } else if (scale == "ordinal") {
        cols <- hcl.colors(n, "viridis")
    } else {
        cols <- if (n <= length(CATEGORICAL)) CATEGORICAL[seq_len(n)] else hcl.colors(n, "Dark 3")
    }
    setNames(cols, levels)
}

# Tree, support values and scale bar
n_tips <- length(tree[["tip.label"]])
show_tips <- n_tips <= 60
p <- ggtree(tree, ladderize = TRUE, linewidth = 0.3)
d <- p[["data"]]
tree_w <- max(d[["x"]])
if (show_tips) p <- p + geom_tiplab(size = 2.2, align = TRUE, linesize = 0.15)
sup <- setNames(support[["support"]], support[["node"]])
int <- d[!d[["isTip"]], c("x", "y", "label")]
int[["support"]] <- unname(sup[int[["label"]]])
int <- int[!is.na(int[["support"]]) & nzchar(int[["support"]]), , drop = FALSE]
# node_support.tsv (emitted above) keeps every value; only the drawn labels are filtered by --minSupportLabel
int <- int[as.numeric(int[["support"]]) >= min_support_label, , drop = FALSE]
if (nrow(int)) {
    p <- p + geom_text(data = int, aes(x = x, y = y, label = support), size = 1.8, hjust = 1.15, vjust = -0.35, colour = "grey30")
}
p <- p + geom_treescale(fontsize = 2.5, linesize = 0.4) +
    labs(caption = sprintf("Branch lengths: substitutions per site; node labels: ultrafast bootstrap support ≥ %g", min_support_label))

# Cluster labels alongside the tree
label_space <- if (show_tips) 0.3 * tree_w else 0.02 * tree_w
clusters <- sort(unique(ann[["cluster"]][!(ann[["cluster"]] %in% c("NA", "unclustered"))]))
for (cl in clusters) {
    members <- ann[["id"]][ann[["cluster"]] == cl]
    node <- if (length(members) > 1) getMRCA(tree, members) else match(members, tree[["tip.label"]])
    p <- p + geom_cladelab(node = node, label = cl, offset = label_space, align = TRUE, fontsize = 2.3,
                           barsize = 0.5, angle = 90, hjust = 0.5, offset.text = 0.01 * tree_w)
}

# One heatmap column per rendered panel, packed left to right with no gaps for omitted panels
col_w <- 0.06 * tree_w
gap <- 0.012 * tree_w
offset <- label_space + (if (length(clusters)) 0.08 * tree_w else 0)
colour_rows <- list()
for (k in seq_along(render)) {
    pn <- render[[k]]
    values <- ann[[pn[["column"]]]]
    lv <- level_order(values)
    pal <- c(palette_for(pn[["scale"]], lv), "NA" = NA_COLOUR)
    colour_rows[[length(colour_rows) + 1]] <- data.frame(panel = pn[["name"]], column = pn[["column"]],
                                                         level = names(pal), colour = unname(pal))
    df <- data.frame(value = factor(values, levels = c(lv, "NA")), row.names = ann[["id"]], check.names = FALSE)
    names(df) <- pn[["name"]]
    if (k > 1) p <- p + new_scale_fill()
    p <- suppressMessages(
        gheatmap(p, df, offset = offset, width = col_w / tree_w, colnames = TRUE, colnames_position = "top",
                 colnames_angle = 90, colnames_offset_y = 0.3, font.size = 2.5, hjust = 0, color = "white") +
        scale_fill_manual(name = pn[["name"]], values = pal, breaks = c(lv, "NA"), drop = FALSE,
                          guide = guide_legend(order = k))
    )
    offset <- offset + col_w + gap
}
p <- p + expand_limits(x = tree_w + offset + 0.02 * tree_w, y = n_tips + 4) +
    theme(legend.position = "right", legend.text = element_text(size = 7), legend.title = element_text(size = 8))

ggsave("phylo_tree.svg", p, width = fig_width, height = fig_height, units = "in")
ggsave("phylo_tree.png", p, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi)

st <- as.data.frame(do.call(rbind, status), stringsAsFactors = FALSE)
colnames(st) <- c("panel", "column", "status", "coverage", "reason")
write.table(st, "panels_status.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)
# Colours of every rendered panel, reused by the interactive tree so both figures match
colours <- if (length(colour_rows)) do.call(rbind, colour_rows) else data.frame(panel = character(), column = character(), level = character(), colour = character())
write.table(colours, "panel_colours.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)
cat(sprintf("PhyloReport: rendered panels: %s\\n", paste(vapply(render, function(x) x[["name"]], ""), collapse = " | ")))
REOF

    # ---- 3. Interactive tree (dropdown filters, added to the index.html dashboard) ----
    python3 - <<'PYEOF'
import csv, html, json
from Bio import Phylo
from plotly.offline import get_plotlyjs_version

SUBTYPE_LABEL = "${subtype_label}"  # tree_id, e.g. "H3N2" (whole genome) or "H3N2_HA" (segment tree): names the file
TREE_FOLDER_LABEL = "${tree_folder}"
# "whole_genome" -> "whole-genome phylogeny"; any segment folder (e.g. "HA") -> "HA segment phylogeny"
TREE_KIND_LABEL = "whole-genome phylogeny" if TREE_FOLDER_LABEL == "whole_genome" else f"{TREE_FOLDER_LABEL} segment phylogeny"
# Bare subtype for display (e.g. "H3N2_HA" -> "H3N2"), so the title/heading don't repeat the segment twice
BARE_SUBTYPE_LABEL = SUBTYPE_LABEL if TREE_FOLDER_LABEL == "whole_genome" else SUBTYPE_LABEL[: -(len(TREE_FOLDER_LABEL) + 1)]

def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\\t"))

tree = Phylo.read("${tree}", "newick")
tree.ladderize()
support = {r["node"]: r["support"] for r in read_tsv("${support}")}
ann = {r["id"]: r for r in read_tsv("annotations.tsv")}
sources = {r["column"]: r for r in read_tsv("annotation_sources.tsv")}
rendered = [r for r in read_tsv("panels_status.tsv") if r["status"] == "rendered"]
colour_of = {}
for r in read_tsv("panel_colours.tsv"):
    colour_of.setdefault(r["column"], {})[r["level"]] = r["colour"]
NA_COLOUR = "#D9D9D9"

# Tree as flat preorder arrays (iterative walk: deep trees do not hit the recursion limit)
parent, length, name, is_tip, node_support = [], [], [], [], []
stack = [(tree.root, -1)]
while stack:
    clade, p = stack.pop()
    i = len(parent)
    parent.append(p)
    length.append(float(clade.branch_length or 0.0))
    name.append(clade.name or "")
    is_tip.append(clade.is_terminal())
    node_support.append("" if clade.is_terminal() else support.get(clade.name, ""))
    for child in reversed(clade.clades):
        stack.append((child, i))

BASE = ["id", "tip_type", "subclade", "cluster", "source_group", "period", "season", "age_group", "sex"]
LABELS = {"subclade": "Subclade", "cluster": "Cluster", "source_group": "Source", "period": "Period",
          "season": "Season", "age_group": "Age group", "sex": "Sex"}
aa_cols = [c for c in (next(iter(ann.values())).keys() if ann else []) if c not in BASE]

hover, filt = {}, {}
for tip, r in ann.items():
    lines = [f"<b>{html.escape(tip)}</b>", "Reference strain" if r["tip_type"] == "root_strain" else "Sample"]
    lines += [f"{LABELS[c]}: {html.escape(r[c])}" for c in LABELS if r[c] != "NA"]
    lines += [f"{c}: {html.escape(r[c])}" for c in aa_cols]
    hover[tip] = "<br>".join(lines)
    filt[tip] = {"sample": r["tip_type"] == "sample", "season": r["season"], "age": r["age_group"], "sex": r["sex"]}

panels = []
for p in rendered:
    col = p["column"]
    values = {tip: r.get(col, "NA") for tip, r in ann.items()}
    cmap = colour_of.get(col, {})
    panels.append({"name": p["panel"], "values": values,
                   "colours": {tip: cmap.get(v, NA_COLOUR) for tip, v in values.items()}})

samples = [r for r in ann.values() if r["tip_type"] == "sample"]

def options(column):
    if sources.get(column, {}).get("available") != "yes":
        return None
    present = sorted({r[column] for r in samples if r[column] != "NA"})
    if not present:
        return None
    return present + (["NA"] if any(r[column] == "NA" for r in samples) else [])

data = {
    "parent": parent, "length": length, "name": name, "tip": is_tip, "support": node_support,
    "hover": hover, "filter": filt, "panels": panels,
    "seasons": [s for s in (options("season") or []) if s != "NA"],
    "ages": options("age_group"), "sexes": options("sex"),
    "subcladeColours": colour_of.get("subclade", {}),
    "subclade": {tip: r["subclade"] for tip, r in ann.items()},
    "naColour": NA_COLOUR,
    "subtypeLabel": BARE_SUBTYPE_LABEL,
    "treeKindLabel": TREE_KIND_LABEL,
    "minSupportLabel": float("${min_support_label}"),
}

def note(column, field):
    src = sources.get(column, {})
    return src.get("note", "") if src.get("available") != "yes" else f"no {field} values"

def select(sel_id, label, opts, missing_note, all_label="All"):
    if not opts:
        return (f'<div class="ctl"><label>{label}</label><br><select id="{sel_id}" disabled title="{html.escape(missing_note)}">'
                f'<option value="ALL">{html.escape(missing_note or "not available")}</option></select></div>')
    items = [f'<option value="ALL">{all_label}</option>'] + [f'<option value="{html.escape(o)}">{html.escape(o)}</option>' for o in opts]
    return f'<div class="ctl"><label>{label}</label><br><select id="{sel_id}">{"".join(items)}</select></div>'

season_opts = data["seasons"]
controls = (select("selFrom", "SEASON FROM", season_opts, note("season", "DATE"), "All Time")
            + (f'<div class="ctl"><label>SEASON TO</label><br><select id="selTo">'
               + "".join(f'<option value="{html.escape(s)}">{html.escape(s)}</option>' for s in season_opts)
               + '</select></div>' if season_opts else "")
            + select("selAge", "AGE GROUP", data["ages"], note("age_group", "AGE GROUP"))
            + select("selSex", "SEX", data["sexes"], note("sex", "SEX")))

legend = "".join(
    f'<div class="leg"><b>{html.escape(p["name"])}</b><br>'
    + "".join(f'<span class="chip"><i style="background:{c}"></i>{html.escape(lvl)}</span>'
              for lvl, c in colour_of.get(next(r["column"] for r in rendered if r["panel"] == p["name"]), {}).items())
    + "</div>" for p in panels)

payload = json.dumps(data).replace("<", "\\\\u003c")

page = r'''<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Phylogenetic Tree - SUBTYPELABEL - TREEKINDLABEL</title>
<script src="https://cdn.plot.ly/plotly-PLOTLYVERSION.min.js"></script>
<style>
body { font-family: Arial, Helvetica, sans-serif; margin: 16px; color: #222; }
.bar { display: flex; flex-wrap: wrap; gap: 12px; align-items: flex-end; margin-bottom: 8px; }
.ctl label { font-size: 10px; font-weight: bold; color: #666; }
.ctl select { padding: 6px; border-radius: 4px; min-width: 150px; border: 1px solid #ccc; background: white; }
#info { font-size: 12px; color: #555; margin: 4px 0 8px 0; }
.legends { display: flex; flex-wrap: wrap; gap: 18px; font-size: 12px; margin-top: 8px; }
.chip { display: inline-block; margin: 2px 8px 2px 0; }
.chip i { display: inline-block; width: 11px; height: 11px; margin-right: 4px; vertical-align: middle; border: 1px solid #999; }
</style></head><body>
<h2 id="title">SUBTYPELABEL TREEKINDLABEL</h2>
<div class="bar">CONTROLS</div>
<div id="info"></div>
<div id="tree"></div>
<div class="legends">LEGENDS</div>
<p style="font-size:11px;color:#777">Maximum likelihood tree inferred once from all samples that passed coverage QC; filters prune it to the
selected samples, keeping the branch lengths of the full tree. Node labels: ultrafast bootstrap support &ge; MINSUPPORTLABEL. Reference strains are always shown.</p>
<script>
const D = PAYLOAD;
const n = D.parent.length;
const kids = Array.from({length: n}, () => []);
for (let i = 1; i < n; i++) kids[D.parent[i]].push(i);
const byId = id => document.getElementById(id);

function wanted(tip) {
  const f = D.filter[tip];
  if (!f || !f.sample) return true;
  const from = byId("selFrom") ? byId("selFrom").value : "ALL";
  if (from !== "ALL") {
    const lo = D.seasons.indexOf(from), hi = D.seasons.indexOf(byId("selTo").value);
    const k = D.seasons.indexOf(f.season);
    if (k < 0 || k < Math.min(lo, hi) || k > Math.max(lo, hi)) return false;
  }
  const age = byId("selAge") ? byId("selAge").value : "ALL";
  if (age !== "ALL" && f.age !== age) return false;
  const sex = byId("selSex") ? byId("selSex").value : "ALL";
  if (sex !== "ALL" && f.sex !== sex) return false;
  return true;
}

function layout() {
  // keep[i]: the node has at least one wanted tip below it (children always follow parents in preorder)
  const keep = new Array(n).fill(false);
  for (let i = n - 1; i >= 0; i--) {
    if (D.tip[i]) keep[i] = wanted(D.name[i]);
    else keep[i] = kids[i].some(c => keep[c]);
  }
  const kept = i => kids[i].filter(c => keep[c]);
  const drawn = [];
  let order = 0;
  function place(i, xParent, parentIdx, extra) {
    let total = extra, j = i;
    while (!D.tip[j] && kept(j).length === 1) { j = kept(j)[0]; total += D.length[j]; }
    const node = {j: j, x: xParent + total, parent: parentIdx, kids: []};
    const idx = drawn.length;
    drawn.push(node);
    if (D.tip[j]) { node.y = order++; }
    else {
      for (const c of kept(j)) node.kids.push(place(c, node.x, idx, D.length[c]));
      node.y = (drawn[node.kids[0]].y + drawn[node.kids[node.kids.length - 1]].y) / 2;
    }
    return idx;
  }
  if (keep[0]) place(0, 0, -1, 0);
  return drawn;
}

function niceScale(maxX) {
  const raw = maxX / 5, p = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 5, 10]) if (raw <= m * p) return m * p;
  return raw;
}

function draw() {
  const drawn = layout();
  const tips = drawn.filter(d => D.tip[d.j]);
  const nTips = tips.length, nSamples = tips.filter(d => D.filter[D.name[d.j]] && D.filter[D.name[d.j]].sample).length;
  const totalSamples = Object.values(D.filter).filter(f => f.sample).length;
  const from = byId("selFrom") ? byId("selFrom").value : "ALL";
  const to = byId("selTo") ? byId("selTo").value : "";
  const label = from === "ALL" ? "All Time" : (from === to ? "Season " + from : "Seasons " + from + " to " + to);
  byId("title").textContent = D.subtypeLabel + " " + D.treeKindLabel + " - " + label;
  byId("info").textContent = "Showing " + nSamples + " of " + totalSamples + " samples, plus " + (nTips - nSamples) + " reference strains.";

  const bx = [], by = [];
  for (const d of drawn) {
    if (d.parent >= 0) { bx.push(drawn[d.parent].x, d.x, null); by.push(d.y, d.y, null); }
    if (d.kids.length) { bx.push(d.x, d.x, null); by.push(drawn[d.kids[0]].y, drawn[d.kids[d.kids.length - 1]].y, null); }
  }
  const maxX = Math.max(1e-9, ...drawn.map(d => d.x));
  const showLabels = nTips <= 80;
  const traces = [{x: bx, y: by, mode: "lines", line: {color: "#333", width: 1}, hoverinfo: "skip", xaxis: "x", yaxis: "y"}];
  traces.push({x: tips.map(d => d.x), y: tips.map(d => d.y), mode: "markers", xaxis: "x", yaxis: "y",
    marker: {size: 6, color: tips.map(d => D.subcladeColours[D.subclade[D.name[d.j]]] || D.naColour), line: {color: "#333", width: 0.5}},
    text: tips.map(d => D.hover[D.name[d.j]]), hovertemplate: "%{text}<extra></extra>"});
  // D.support keeps every node's raw value; only the drawn labels are filtered by --minSupportLabel
  const sup = drawn.filter(d => !D.tip[d.j] && D.support[d.j] !== "" && d.parent >= 0 && Number(D.support[d.j]) >= D.minSupportLabel);
  traces.push({x: sup.map(d => d.x), y: sup.map(d => d.y), mode: "text", text: sup.map(d => D.support[d.j]),
    textposition: "top left", textfont: {size: 9, color: "#666"}, hoverinfo: "skip", xaxis: "x", yaxis: "y"});
  if (showLabels) traces.push({x: tips.map(() => maxX * 1.02), y: tips.map(d => d.y), mode: "text", text: tips.map(d => D.name[d.j]),
    textposition: "middle right", textfont: {size: 10}, hoverinfo: "skip", xaxis: "x", yaxis: "y"});

  const K = D.panels.length;
  const height = Math.max(450, nTips * 16 + 200);
  // Heatmap cells are bars, one row tall (base y-0.5, height 1), so a column is a solid band whatever the number of
  // tips; colWidth (x2 units, column pitch = 1) leaves only a thin gap between neighbouring columns
  const colWidth = 0.9;
  D.panels.forEach((p, k) => {
    traces.push({type: "bar", x: tips.map(() => k), y: tips.map(() => 1), base: tips.map(d => d.y - 0.5),
      width: colWidth, xaxis: "x2", yaxis: "y",
      marker: {color: tips.map(d => p.colours[D.name[d.j]] || D.naColour), line: {width: 0}},
      hovertext: tips.map(d => "<b>" + D.name[d.j] + "</b><br>" + p.name + ": " + (p.values[D.name[d.j]] || "NA")),
      hovertemplate: "%{hovertext}<extra></extra>"});
  });
  // Width of the heatmap area: each column gets a fixed share of the plot, capped so wide panel sets still leave room
  const heatFrac = Math.min(0.5, 0.05 * K + 0.02);
  const scale = niceScale(maxX);
  Plotly.react("tree", traces, {
    height: height, showlegend: false, hovermode: "closest", barmode: "overlay", margin: {t: 110, l: 20, r: 20, b: 60},
    xaxis: {domain: [0, 1 - heatFrac - 0.02], range: [-maxX * 0.02, maxX * (showLabels ? 1.45 : 1.05)],
            showgrid: false, zeroline: false, showticklabels: false},
    xaxis2: {domain: [1 - heatFrac, 1], range: [-0.5, Math.max(K, 1) - 0.5], side: "top", tickangle: -90,
             tickvals: D.panels.map((p, k) => k), ticktext: D.panels.map(p => p.name), showgrid: false, zeroline: false},
    yaxis: {autorange: "reversed", showgrid: false, zeroline: false, showticklabels: false, range: [nTips + 1, -1]},
    shapes: [{type: "line", xref: "x", yref: "y", x0: 0, x1: scale, y0: nTips + 0.5, y1: nTips + 0.5, line: {width: 2}}],
    annotations: [{xref: "x", yref: "y", x: scale / 2, y: nTips + 0.5, yshift: -14, showarrow: false,
                   text: scale + " substitutions/site", font: {size: 11}}]
  }, {responsive: true});
}

function onFrom() {
  const from = byId("selFrom"), to = byId("selTo");
  if (to) {
    to.disabled = from.value === "ALL";
    if (from.value !== "ALL" && D.seasons.indexOf(to.value) < D.seasons.indexOf(from.value)) to.value = from.value;
  }
  draw();
}
if (byId("selFrom")) byId("selFrom").addEventListener("change", onFrom);
for (const id of ["selTo", "selAge", "selSex"]) if (byId(id)) byId(id).addEventListener("change", draw);
if (byId("selTo")) byId("selTo").disabled = true;
draw();
</script></body></html>
'''
page = (page.replace("PLOTLYVERSION", get_plotlyjs_version()).replace("CONTROLS", controls)
            .replace("LEGENDS", legend).replace("PAYLOAD", payload).replace("SUBTYPELABEL", BARE_SUBTYPE_LABEL)
            .replace("TREEKINDLABEL", TREE_KIND_LABEL).replace("MINSUPPORTLABEL", "${min_support_label}"))
with open(f"PhylogeneticTreeReport_{SUBTYPE_LABEL}.html", "w", encoding="utf-8") as f:
    f.write(page)
print(f"PhyloReport: interactive tree written with {len(panels)} panels and filters: "
      + ", ".join(k for k, v in (("season", season_opts), ("age group", data["ages"]), ("sex", data["sexes"])) if v) + ".")
PYEOF

    # ---- 4. Run manifest ----
    python3 - <<'PYEOF'
import base64, csv, hashlib, json, os

run = json.loads(base64.b64decode("${run_info_b64}").decode("utf-8"))

def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\\t"))

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

inputs = []
for name in sorted(os.listdir("inputs")) if os.path.isdir("inputs") else []:
    path = os.path.join("inputs", name)
    if os.path.isfile(path):
        inputs.append({"file": name, "sha256": sha256(path)})

qc = {}
for row in read_tsv("${qc_table}"):
    entry = qc.setdefault(row["sample_id"], [])
    if row["status"] != "PASS":
        entry.append(f"{row['segment']}: {row['reason']}")
dropped = [{"sample": s, "reasons": r} for s, r in sorted(qc.items()) if r]

cluster_rows = read_tsv("${clusters}")
cluster_sizes = {}
for row in cluster_rows:
    if row["cluster"] not in ("NA", "unclustered"):
        cluster_sizes[row["cluster"]] = cluster_sizes.get(row["cluster"], 0) + 1

panel_rows = read_tsv("panels_status.tsv")
params = run["parameters"]
with open("${versions}") as f:
    tools = json.load(f)

manifest = {
    "module": "FluTyper phylogenetics",
    "pipeline_version": run.get("pipeline_version"),
    "started": run.get("started"),
    "parameters": params,
    "tools": tools,
    "inputs": inputs,
    "methods": {
        "reference_spec": "Koumaty et al. 2026, Microbial Genomics 12:001810 (doi:10.1099/mgen.0.001810)",
        "coverage_rule": f"every segment covered (non-N, against the coordinate reference) for at least {params['phyloMinCoverage']} of its length",
        "coordinate_references": read_tsv("${coordinate_table}"),
        "root_strains": read_tsv("${root_table}"),
        "alignment": "MAFFT --auto per segment, concatenated in phyloSegmentOrder",
        "alignment_curation": "none (no manual curation); columns inserted relative to the coordinate reference are removed",
        "alignment_segments": read_tsv("${segments}"),
        "tree": f"IQ-TREE model {params['model']}, {params['bootstrap']} ultrafast bootstrap replicates, seed {params['seed']}, rooted on {params['phyloOutgroup']}",
        "subclades": "Nextclade column 'subclade' (samples: GenotypingNextclade results; root strains: same dataset)",
        "aa_state_source": "TreeTime parsimony ancestral reconstruction (internal nodes) translated through the reference codon map; tips use their observed sequence, ambiguous codons are NA",
        "clustering": f"TreeCluster single_linkage within each subclade, threshold {params['clusterThreshold']} substitutions/site",
    },
    "qc": {"candidates": len(qc), "passed": len(qc) - len(dropped), "dropped": dropped},
    "clusters": {"n_clusters": len(cluster_sizes), "sizes": dict(sorted(cluster_sizes.items())),
                 "unclustered": sum(1 for r in cluster_rows if r["cluster"] == "unclustered")},
    "aa_positions": read_tsv("${positions}"),
    "panels": {
        "rendered": [r["panel"] for r in panel_rows if r["status"] == "rendered"],
        "skipped": [{"panel": r["panel"], "reason": r["reason"], "coverage": r["coverage"]} for r in panel_rows if r["status"] == "skipped"],
    },
}
with open("run_manifest.json", "w") as f:
    json.dump(manifest, f, indent=2)
PYEOF

    # ---- Move this tree's files into its own folder, so publishing several trees of the same subtype never collides ----
    mkdir -p "${tree_folder}"
    mv annotations.tsv annotation_sources.tsv aa_states.tsv phylo_tree.svg phylo_tree.png panels_status.tsv \\
       run_manifest.json "PhylogeneticTreeReport_${subtype_label}.html" "${tree_folder}/"
    """
}
