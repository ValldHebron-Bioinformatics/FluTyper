#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloAuspice {
    // Nextstrain/Auspice dataset (augur export v2 JSON) for one tree of the phylogenetics module (plan
    // .claude/plans/auspice-integration.md, Phase 4), built from the tree's own for_auspice bundle (no join by tree_id):
    //  1. meta.tsv, one row per tip (plan D9, D12): subclade, cluster, season, age_group and sex from PhyloReport's
    //     annotations.tsv (one source of truth); DATE and LOCATION (+ optional COUNTRY) from --metadata; country,
    //     division, region and the Auspice location label from PhyloGeography's geo.tsv, joined on (LOCATION, COUNTRY
    //     or NA). Metadata is never imputed: a missing value is written '?' (augur's missing value; '' would be a
    //     state) and every unresolved geography value, and every value without a map point, is a WARNING. Dates stay
    //     at their precision (YYYY-MM -> YYYY-MM-XX: TreeTime infers within the month, the date is never completed).
    //     Root strains have no metadata: their date is YYYY-XX-XX when the last '/' field of the strain name
    //     (root_strains.tsv) is a 4-digit year, else undated with a WARNING; a root_strain_dates.tsv row (sourced
    //     collection date) wins.
    //  2. augur refine --timetree with a fixed clock (D8): the subtype x segment rate of --auspiceClockRates
    //     (nextstrain/seasonal-flu core.smk@a1e704c7; MP/NS generic 0.001 -> WARNING), std-dev = rate/5, no clock
    //     filter; the whole-genome rate is the length-weighted mean of the 8 segment rates (logged). With no rate for
    //     the tree, fewer tips dated to at least the month than --auspiceMinDatedTips, or a failed timetree, the
    //     export falls back to a divergence-only tree with a WARNING (never a guessed rate). TreeTime 0.12.1 always
    //     collapses zero-length branches (Phase 0; augur has no flag); D6 accepts that only as a pure coarsening:
    //     identical tips, every surviving internal node is an input node with the identical tip set, and every
    //     removed branch has 0 nucleotide changes in FluTyper's reconstruction. Anything else fails the tree loudly.
    //  3. nt_muts.json / aa_muts.json (D7): branch mutations = parent->child differences on the REFINED tree from
    //     FluTyper's own TreeTime parsimony reconstruction (ancestral_nt.fasta; ancestral_aa/ with FluTyper gene names;
    //     tips = observed alignment; non-ACGT and '-'/'X' skipped), so Auspice shows the static report's mutations.
    //     The annotations block (entropy panel) comes from this tree's gene_map.tsv rows: alignment column =
    //     segment start - 1 + segment nt; a CDS in several parts (spliced M2/NS2, frameshifted PA-X) is written as
    //     augur "segments" (schema-annotations.json). generated_by carries the FluTyper version (PEP 440).
    //  4. augur traits --confidence on --auspiceTraits (default country,region; D10, D13): discrete trait analysis
    //     (DTA), exploratory and sampling-biased. augur imputes tip states for '?' tips; those tip entries are
    //     removed afterwards, internal-node states are kept (they draw the transmission lines).
    //  5. augur export v2 -> flutyper_<subtype>_<segment>.json (lower case, 'genome' for the whole-genome tree, no
    //     '_' inside a part: Auspice maps '_' to URL levels), then augur validate export-v2. The config is adapted
    //     per dataset: title, divergence-only mode (distance 'div', no num_date colouring), and explicit legend/scale
    //     for age_group (ordered by the lower age bound: the ordinal type alone sorts "0-2", "15-65", ">65" as text)
    //     and season.
    //  6. Bootstrap support (plan D18): the IQ-TREE ultrafast bootstrap of every internal node that survives refine is a
    //     hover attribute ('support', exported through the config's metadata_columns from an export-only metadata table,
    //     so it is not a colouring), and a branch label (Auspice "Branch labels" dropdown) when >= --minSupportLabel
    //     (0 labels every node that has a value). Low support never removes or collapses a branch; the only collapses
    //     are D6's zero-length ones, whose support values are then not shown (logged). No value -> no attribute.
    // Every message goes to auspice.log (step_log -> the tree's phylo.log), prefixed PhyloAuspice, INFO/WARNING/ERROR,
    // including the augur and treetime versions actually used (transitive dependencies are not pinned).
    // Failures are reported, never silent: an expected failure (D6 check, refine/export/validate error, bad input)
    // still exits 0 so that auspice.log and auspice_status.tsv are emitted (errorStrategy 'ignore' would discard
    // them), with status 'failed', an ERROR line and no JSON. Only a crash before Python runs (e.g. the conda env)
    // emits nothing: the caller must also report trees with no status (expected ids minus emitted ids).
    conda "${moduleDir}/../../../envs/nextstrain.yaml"
    errorStrategy 'ignore'

    input:
    // PhyloReport's for_auspice bundle, as-is: tree_id "H3N2" (whole genome) or "H3N2_HA"; tree_folder
    // "whole_genome" or the bare segment
    tuple val(tree_id), val(tree_folder), path(tree), path(alignment), path(ancestral_nt), path(ancestral_aa),
          path(segments), path(gene_map), path(annotations), path(support)
    // Shared, once-per-run inputs, each staged in its own folder so that no user file can collide with a bundle file
    // or an intermediate (meta.tsv, refined.nwk, ...)
    path(geo, stageAs: "in_geo/*")                  // PhyloGeography geo.tsv
    path(lat_longs, stageAs: "in_lat_longs/*")      // PhyloGeography lat_longs.tsv (augur --lat-longs)
    path(metadata, stageAs: "in_metadata/*")        // --metadata CSV (ID, DATE, LOCATION, optional COUNTRY), or []
    path(auspice_config, stageAs: "in_config/*")    // --auspiceConfig (RESOURCES/auspice_config.json)
    path(clock_rates, stageAs: "in_clock/*")        // --auspiceClockRates (RESOURCES/auspice_clock_rates.tsv)
    path(root_strains, stageAs: "in_roots/*")       // PhyloPreflight root_strains.tsv (strain <-> tip_id)
    path(root_dates, stageAs: "in_root_dates/*")    // RESOURCES/root_strain_dates.tsv (override), or []
    val(min_dated_tips)                             // --auspiceMinDatedTips: tips dated to >= the month for a time tree
    val(traits)                                     // --auspiceTraits: comma-separated DTA columns (default country,region)
    val(flutyper_version)                           // params.version (generated_by of the node-data JSONs)
    val(min_support_label)                          // --minSupportLabel (default 70; 0 = label every node with a value)

    output:
    tuple val(tree_id), path("flutyper_*.json"), optional: true, emit: json    // only when exported AND validated
    tuple val(tree_id), val(tree_folder), path("auspice_status.tsv"), emit: status  // tree_id status(ok|failed) mode dataset reason
    tuple val(tree_id), path("auspice.log"), emit: step_log   // this step's messages, collected into the tree's phylo.log

    script:
    """
    python3 - <<'PYEOF'
TREE_ID, TREE_FOLDER = "${tree_id}", "${tree_folder}"
TREE, ALN, ANC_NT, ANC_AA = "${tree}", "${alignment}", "${ancestral_nt}", "${ancestral_aa}"
SEGMENTS, GENE_MAP, ANNOT = "${segments}", "${gene_map}", "${annotations}"
GEO, LATLONGS, META = "${geo}".strip(), "${lat_longs}".strip(), "${metadata}".strip()
CONFIG, CLOCK = "${auspice_config}", "${clock_rates}"
ROOTS, ROOT_DATES = "${root_strains}".strip(), "${root_dates}".strip()
MIN_DATED = int("${min_dated_tips}")
TRAITS, FLUTYPER_VERSION = "${traits}", "${flutyper_version}"
SUPPORT, MIN_SUPPORT = "${support}", float("${min_support_label}")
import csv, glob, json, os, re, subprocess, sys, traceback
from datetime import date as Date
import pandas as pd
from Bio import Phylo, SeqIO

P = "PhyloAuspice"
MISSING = "?"        # augur: '?' = missing (traits does not treat it as a state, export drops it); '' would be a state
NT_OK = set("ACGT")
AA_SKIP = set("-X")
ORDINAL_COLOURS = ["#440154", "#482878", "#3E4A89", "#31688E", "#26828E", "#1F9E89", "#35B779", "#6DCD59",
                   "#B4DE2C", "#FDE725"]   # viridis, 10 steps: ordinal legends (age_group, season)

log_lines = []
def say(level, msg):
    log_lines.append(f"{P}: {level} {msg}")

class Fail(Exception):
    pass

status = {"status": "failed", "mode": "NA", "dataset": "NA", "reason": ""}

def tail(path, n=6):
    try:
        with open(path, errors="replace") as f:
            lines = [l.strip() for l in f if l.strip()]
    except OSError:
        return ""
    return " | ".join(lines[-n:])

def run(cmd, logfile):
    with open(logfile, "w") as f:
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode

def tool_version(tool):
    try:
        r = subprocess.run([tool, "--version"], capture_output=True, text=True)
        out = (r.stdout + r.stderr).strip().splitlines()
        return out[-1].strip() if out else f"{tool} (version unknown)"
    except OSError as e:
        return f"{tool} (not runnable: {e})"

def read_table(path):
    # TSV with a header; '#' lines (comments) and blank lines skipped
    with open(path, encoding="utf-8") as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\\t"))

def na(v):
    # the NA spellings of the upstream tables -> None (never a value)
    if v is None:
        return None
    v = str(v).strip()
    return None if v in ("", "NA", "?") else v

def part(s):
    # one part of the dataset name: lower case, and no '_' (Auspice turns '_' into a URL level)
    return re.sub(r"[^a-z0-9.-]+", "-", s.lower()).strip("-") or "unknown"


def main():
    # ---- Tree identity and dataset name ------------------------------------------------------------------------------
    whole = TREE_FOLDER == "whole_genome"
    if whole:
        subtype, segment = TREE_ID, None
    elif TREE_ID.endswith("_" + TREE_FOLDER):
        subtype, segment = TREE_ID[: -len(TREE_FOLDER) - 1], TREE_FOLDER
    else:
        raise Fail(f"tree_id '{TREE_ID}' does not end with its tree_folder '_{TREE_FOLDER}'")
    label = f"{subtype} {segment or 'whole genome'}"
    dataset = f"flutyper_{part(subtype)}_{part(segment or 'genome')}.json"
    status["dataset"] = dataset
    say("INFO", f"{tool_version('augur')}, {tool_version('treetime')} (Nextstrain env; transitive dependencies are not "
                "pinned, so the versions actually used are logged).")
    say("INFO", f"Auspice export of the {label} tree -> {dataset}.")

    # ---- Tree, tip roles ---------------------------------------------------------------------------------------------
    tree = Phylo.read(TREE, "newick")
    tips = [t.name for t in tree.get_terminals()]
    internal = [c.name for c in tree.get_nonterminals()]
    if any(n is None for n in tips + internal) or len(set(tips + internal)) != len(tips + internal):
        raise Fail("the input tree has unnamed or duplicated node names; augur and the D6 check need unique names")

    ann = {}
    ann_cols = []
    if ANNOT:
        rows = read_table(ANNOT)
        ann = {r["id"]: r for r in rows}
        ann_cols = list(rows[0].keys()) if rows else []
    roots = {}       # tip_id -> strain name
    if ROOTS:
        for r in read_table(ROOTS):
            roots[r["tip_id"].strip()] = r["strain"].strip()
    is_root = {t: t in roots or (ann.get(t, {}).get("tip_type") == "root_strain") for t in tips}
    samples = [t for t in tips if not is_root[t]]

    # ---- Metadata (DATE, LOCATION, COUNTRY) and geography (PhyloGeography geo.tsv) -----------------------------------
    meta = {}
    cols = {}
    if META:
        df = pd.read_csv(META, dtype=str, skipinitialspace=True)
        df.columns = [str(c).strip() for c in df.columns]
        id_col = next((c for c in df.columns if c.lower().replace("_", "").replace(" ", "") in ("id", "sampleid")), None)
        if id_col is None:
            raise Fail(f"no sample ID column ('ID', 'Sample_ID' or 'SampleID') in the metadata; columns: {list(df.columns)}")
        df = df[df[id_col].notna()].copy()
        df[id_col] = df[id_col].astype(str).str.strip()
        df = df.drop_duplicates(subset=[id_col], keep="first")    # first row wins, as MetadataMerge / PhyloReport
        for _, row in df.iterrows():
            meta[row[id_col]] = row
        cols = {c.upper(): c for c in df.columns}
    else:
        say("WARNING", "no --metadata: samples have no date, location or country.")

    def mval(tip, field):
        row = meta.get(tip)
        col = cols.get(field)
        if row is None or col is None or pd.isna(row[col]):
            return None
        return str(row[col]).strip() or None

    geo = {}
    if GEO:
        for r in read_table(GEO):
            geo[(r["location"], r["country_input"])] = r

    no_meta_row = [t for t in samples if META and t not in meta]
    if no_meta_row:
        say("WARNING", f"{len(no_meta_row)} sample tip(s) have no --metadata row (all their metadata stays NA): "
                       f"{', '.join(no_meta_row[:10])}{' ...' if len(no_meta_row) > 10 else ''}")

    # ---- Dates -------------------------------------------------------------------------------------------------------
    def augur_date(raw):
        # -> (augur date string or None, precision). Partial dates stay partial, written as augur's explicit
        # ambiguity (YYYY-MM-XX / YYYY-XX-XX: TreeTime infers within the range); never completed to a day.
        m = re.fullmatch(r"(\\d{4})-(\\d{2}|XX)-(\\d{2}|XX)", raw)
        if m:
            y, mo, d = m.groups()
            try:
                if mo != "XX" and d != "XX":
                    Date(int(y), int(mo), int(d))
                    return raw, "day"
                if mo != "XX" and d == "XX" and 1 <= int(mo) <= 12:
                    return raw, "month"
                if mo == "XX" and d == "XX":
                    return raw, "year"
            except ValueError:
                pass
            return None, None
        m = re.fullmatch(r"(\\d{4})-(\\d{2})", raw)
        if m:
            return (f"{raw}-XX", "month") if 1 <= int(m.group(2)) <= 12 else (None, None)
        if re.fullmatch(r"\\d{4}", raw):
            return f"{raw}-XX-XX", "year"
        d = pd.to_datetime(raw, errors="coerce")      # same parser as PhyloReport's season
        if pd.isna(d):
            return None, None
        return f"{d.year:04d}-{d.month:02d}-{d.day:02d}", "day"

    override = {}
    if ROOT_DATES:
        for r in read_table(ROOT_DATES):
            strain, value = (r.get("strain") or "").strip(), (r.get("date") or "").strip()
            if not strain:
                continue
            parsed, prec = augur_date(value) if value else (None, None)
            if parsed is None or not re.fullmatch(r"\\d{4}-(\\d{2}|XX)-(\\d{2}|XX)", value):
                say("WARNING", f"root_strain_dates row '{strain}': date '{value}' is not YYYY-MM-DD / YYYY-MM-XX / "
                               "YYYY-XX-XX: row ignored.")
                continue
            override[strain] = (parsed, prec, (r.get("source") or "").strip() or "no source given")

    dates, precision = {}, {}
    unparsed = []
    for t in tips:
        if is_root[t]:
            strain = roots.get(t)
            if strain is None:
                say("WARNING", f"root strain {t} is not in root_strains.tsv (no strain name): left undated.")
                continue
            if strain in override:
                dates[t], precision[t], src = override[strain]
                say("INFO", f"root strain {t} ({strain}): date {dates[t]} from {os.path.basename(ROOT_DATES)} (source: {src}).")
                continue
            last = strain.split("/")[-1].strip()
            if re.fullmatch(r"\\d{4}", last):
                dates[t], precision[t] = f"{last}-XX-XX", "year"
                say("INFO", f"root strain {t} ({strain}): date {dates[t]} from the strain name (ambiguous within the year).")
            else:
                say("WARNING", f"root strain {t} ({strain}): the last '/' field is not a 4-digit year: left undated "
                               "(TreeTime infers it); add a sourced row to root_strain_dates.tsv to date it.")
            continue
        raw = mval(t, "DATE")
        if raw is None:
            continue
        d, p = augur_date(raw)
        if d is None:
            unparsed.append(f"{t}={raw}")
            continue
        dates[t], precision[t] = d, p
    if unparsed:
        say("WARNING", f"{len(unparsed)} DATE value(s) could not be parsed and are left undated: {', '.join(unparsed[:10])}")
    n_prec = {k: sum(1 for t in tips if precision.get(t) == k) for k in ("day", "month", "year")}
    n_dated = n_prec["day"] + n_prec["month"]
    say("INFO", f"{len(tips)} tips: {n_prec['day']} dated to the day, {n_prec['month']} to the month, {n_prec['year']} "
                f"to the year only, {len(tips) - len(dates)} undated.")

    # ---- meta.tsv ----------------------------------------------------------------------------------------------------
    ANN_FIELDS = ["subclade", "cluster", "season", "age_group", "sex"]
    GEO_FIELDS = ["region", "country", "division", "location"]
    for c in ANN_FIELDS:
        if ann and c not in ann_cols:
            say("INFO", f"annotations.tsv has no '{c}' column: '{c}' is NA for every tip.")
    values = {t: {} for t in tips}
    geo_missing = {c: {} for c in GEO_FIELDS}    # field -> reason -> n sample tips
    for t in tips:
        a = ann.get(t, {})
        for c in ANN_FIELDS:
            values[t][c] = na(a.get(c))
        if is_root[t]:
            for c in GEO_FIELDS:
                values[t][c] = None
            continue
        city = mval(t, "LOCATION")
        given = mval(t, "COUNTRY")
        g = geo.get((city or "NA", given or "NA"))
        if g is None and (city is not None or given is not None):
            say("WARNING", f"{t}: (LOCATION '{city}', COUNTRY '{given}') is not in geo.tsv: location kept as written, "
                           "country/division/region NA.")
            values[t].update(location=city, country=None, division=None, region=None)
        elif g is None:
            values[t].update(location=None, country=None, division=None, region=None)
        else:
            values[t].update(location=na(g["auspice_location"]), country=na(g["country"]),
                             division=na(g["division"]), region=na(g["region"]))
        for c in GEO_FIELDS:
            if values[t][c] is None:
                reason = f"LOCATION '{city}'" if city else "no LOCATION"
                geo_missing[c][reason] = geo_missing[c].get(reason, 0) + 1
    for c in GEO_FIELDS:
        if geo_missing[c]:
            n = sum(geo_missing[c].values())
            say("WARNING", f"{c}: {n} sample tip(s) have no value ({', '.join(f'{k}: {v}' for k, v in sorted(geo_missing[c].items()))}): "
                           "left NA, never imputed (the DTA gives these tips no state).")
    n_roots = sum(1 for t in tips if is_root[t])
    if n_roots:
        say("INFO", f"{n_roots} root strain tip(s) have no geography, season, age group or sex (reference strains, not imputed).")

    columns = ["strain", "date"] + ANN_FIELDS + GEO_FIELDS
    with open("meta.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
        w.writerow(columns)
        for t in tips:
            w.writerow([t, dates.get(t, "")] + [values[t][c] if values[t][c] is not None else MISSING
                                                for c in ANN_FIELDS + GEO_FIELDS])

    # ---- Map points: every value shown at a map level needs a lat_longs row -----------------------------------------
    cfg = json.load(open(CONFIG, encoding="utf-8"))
    ll = set()
    if LATLONGS:
        with open(LATLONGS, encoding="utf-8") as f:
            for line in f:
                p = line.rstrip("\\n").split("\\t")
                if len(p) >= 4 and not line.startswith("#"):
                    ll.add((p[0].strip().lower(), p[1].strip().lower()))
    else:
        say("WARNING", "no lat_longs.tsv: the map has no points.")
    no_coords = set()     # (level, value) already WARNed here: augur export's own messages about them are not repeated
    for level in [g["key"] for g in cfg.get("geo_resolutions", [])]:
        used = {}
        for t in samples:
            v = values[t].get(level)
            if v is not None:
                used[v] = used.get(v, 0) + 1
        for v, n in sorted(used.items()):
            if (level, v.lower()) not in ll:
                no_coords.add((level, v.lower()))
                say("WARNING", f"{level} '{v}' ({n} tip{'s' if n > 1 else ''}) has no coordinates in lat_longs.tsv: "
                               f"no map point at the {level} level.")

    # ---- Clock rate (D8) ---------------------------------------------------------------------------------------------
    rates = {}
    for r in read_table(CLOCK):
        try:
            rate = float(r["rate"])
        except (TypeError, ValueError):
            raise Fail(f"clock-rate table {os.path.basename(CLOCK)}: rate '{r.get('rate')}' of {r.get('subtype')} {r.get('segment')} is not a number")
        if not rate > 0:
            raise Fail(f"clock-rate table {os.path.basename(CLOCK)}: rate of {r['subtype']} {r['segment']} must be > 0")
        rates[(r["subtype"].strip(), r["segment"].strip())] = (rate, (r.get("generic") or "").strip().lower() == "yes",
                                                               (r.get("source") or "").strip())
    clock, why_no_time = None, None
    if whole:
        segs = [(r["segment"], int(r["length"])) for r in read_table(SEGMENTS)]
        lacking = [s for s, _ in segs if (subtype, s) not in rates]
        if lacking:
            why_no_time = (f"no clock rate for {subtype} {', '.join(lacking)} in {os.path.basename(CLOCK)}: divergence-only "
                           "tree (the whole-genome rate is the length-weighted mean of every segment's rate; no rate is guessed)")
        else:
            total = sum(n for _, n in segs)
            clock = sum(rates[(subtype, s)][0] * n for s, n in segs) / total
            say("INFO", f"clock rate {clock:.6g} (std-dev {clock / 5:.6g}) substitutions/site/year for {label}: "
                        f"length-weighted mean of {', '.join(f'{s} {rates[(subtype, s)][0]:g} x {n} nt' for s, n in segs)}.")
            generic = [s for s, _ in segs if rates[(subtype, s)][1]]
            if generic:
                say("WARNING", f"clock rate of {', '.join(generic)}: generic rate, not segment-specific "
                               f"({'; '.join(sorted({rates[(subtype, s)][2] for s in generic}))}); it enters the whole-genome mean.")
    else:
        if (subtype, segment) not in rates:
            why_no_time = (f"no clock rate for {subtype} {segment} in {os.path.basename(CLOCK)}: divergence-only tree "
                           "(no rate is guessed; add a row to --auspiceClockRates)")
        else:
            clock, generic, source = rates[(subtype, segment)]
            say("INFO", f"clock rate {clock:.6g} (std-dev {clock / 5:.6g}) substitutions/site/year for {label} ({source}).")
            if generic:
                say("WARNING", f"clock rate of {segment}: generic rate, not segment-specific ({source}).")
    if why_no_time is None and n_dated < MIN_DATED:
        why_no_time = (f"only {n_dated} tip(s) dated to at least the month (minimum {MIN_DATED}): divergence-only tree")

    # ---- augur refine (D2, D6, D8) -----------------------------------------------------------------------------------
    mode = "divergence"
    if why_no_time is None:
        rc = run(["augur", "refine", "--tree", TREE, "--alignment", ALN, "--metadata", "meta.tsv",
                  "--timetree", "--keep-root", "--keep-polytomies", "--date-inference", "marginal", "--date-confidence",
                  "--coalescent", "opt", "--clock-rate", repr(clock), "--clock-std-dev", repr(clock / 5),
                  "--divergence-units", "mutations-per-site", "--seed", "1",
                  "--output-tree", "refined.nwk", "--output-node-data", "branch_lengths.json"], "refine.log")
        if rc == 0:
            mode = "timetree"
        else:
            why_no_time = f"augur refine --timetree failed (exit {rc}: {tail('refine.log', 3)}): divergence-only tree"
    if mode != "timetree":
        say("WARNING", why_no_time)
        rc = run(["augur", "refine", "--tree", TREE, "--metadata", "meta.tsv", "--keep-root", "--keep-polytomies",
                  "--divergence-units", "mutations-per-site",
                  "--output-tree", "refined.nwk", "--output-node-data", "branch_lengths.json"], "refine.log")
        if rc != 0:
            raise Fail(f"augur refine (divergence only) failed, exit {rc}: {tail('refine.log')}")
    status["mode"] = mode

    # D6: accept only a pure coarsening of the input tree (TreeTime collapses zero-length branches)
    refined = Phylo.read("refined.nwk", "newick")
    anc = {r.id: str(r.seq).upper() for r in SeqIO.parse(ANC_NT, "fasta")}
    observed = {r.id: str(r.seq).upper() for r in SeqIO.parse(ALN, "fasta")}
    def clade_map(t):
        return {c.name: frozenset(x.name for x in c.get_terminals()) for c in t.get_nonterminals()}
    problems = []
    r_tips = [t.name for t in refined.get_terminals()]
    if sorted(r_tips) != sorted(tips):
        problems.append(f"tips differ (lost: {sorted(set(tips) - set(r_tips))[:5]}, new: {sorted(set(r_tips) - set(tips))[:5]})")
    in_clades, out_clades = clade_map(tree), clade_map(refined)
    for name, members in out_clades.items():
        if name not in in_clades:
            problems.append(f"internal node {name} is not in the input tree")
        elif in_clades[name] != members:
            problems.append(f"internal node {name} changed its tip set")
    parent = {c.name: p.name for p in tree.find_clades() for c in p.clades}
    removed = [n for n in internal if n not in out_clades]
    changed = []
    for n in removed:
        if n not in anc or parent.get(n) not in anc:
            problems.append(f"collapsed node {n} or its parent has no reconstructed sequence")
            continue
        k = sum(1 for a, b in zip(anc[parent[n]], anc[n]) if a != b and a in NT_OK and b in NT_OK)
        if k:
            changed.append(f"{n}: {k} reconstructed nucleotide change(s) on its branch from {parent[n]}")
    if changed:
        problems.append("collapsed branch(es) carry changes, " + "; ".join(changed[:10]))
    if problems:
        raise Fail("D6 coarsening check failed: augur refine did not return a pure coarsening of the input tree: "
                   + " | ".join(problems) + ". Not exported: branch mutations would be misplaced.")
    if removed:
        shown = ", ".join(removed[:10]) + (f" and {len(removed) - 10} more" if len(removed) > 10 else "")
        say("INFO", f"augur refine collapsed {len(removed)} zero-length internal branch(es) ({shown}), all with 0 "
                    "reconstructed nucleotide changes (D6: pure coarsening; shown as polytomies in Auspice).")
    else:
        say("INFO", "augur refine kept every internal node of the input tree.")
    bl = json.load(open("branch_lengths.json"))["nodes"]
    if mode == "timetree":
        root = refined.root.name
        def when(n):
            d = bl.get(n, {})
            ci = d.get("num_date_confidence") or [None, None]
            return f"{d.get('date', '?')} [{ci[0]:.2f}, {ci[1]:.2f}]" if ci[0] is not None else d.get("date", "?")
        msg = f"time tree: root {root} dated {when(root)}"
        dated_samples = [t for t in samples if t in r_tips]
        if len(dated_samples) >= 2:
            mrca = refined.common_ancestor(dated_samples).name
            msg += f"; MRCA of the {len(dated_samples)} samples {mrca} dated {when(mrca)}"
        say("INFO", msg + " (date, confidence interval in decimal years).")

    # ---- Branch mutations from FluTyper's reconstruction (D7) --------------------------------------------------------
    nt = dict(anc)
    nt.update(observed)      # tips: observed alignment, as PhyloPostTree
    edges = [(p.name, c.name) for p in refined.find_clades() for c in p.clades]
    lacking = sorted({n for e in edges for n in e if n not in nt})
    if lacking:
        raise Fail(f"nodes of the refined tree without a nucleotide sequence: {', '.join(lacking[:10])}")
    L = len(next(iter(observed.values())))
    def diffs(a, b, keep):
        return [f"{x}{i}{y}" for i, (x, y) in enumerate(zip(a, b), start=1) if x != y and keep(x) and keep(y)]
    nodes_nt = {refined.root.name: {"muts": []}}
    for p, c in edges:
        nodes_nt[c] = {"muts": diffs(nt[p], nt[c], lambda s: s in NT_OK)}

    seg_start = {r["segment"]: int(r["start"]) for r in read_table(SEGMENTS)}
    gene_cols, unmapped = {}, {}
    for r in read_table(GENE_MAP):
        if r["segment"] not in seg_start:
            continue
        if r["ref_aa"] == "?":
            unmapped[r["gene"]] = unmapped.get(r["gene"], 0) + 1
            continue
        base = seg_start[r["segment"]] - 1          # 1-based alignment column = segment start - 1 + segment nt
        gene_cols.setdefault(r["gene"], []).extend(base + int(r[k]) for k in ("nt1", "nt2", "nt3"))
    annotations = {"nuc": {"start": 1, "end": L, "strand": "+", "type": "source"}}
    for gene, pos in gene_cols.items():
        runs, s = [], pos[0]
        for a, b in zip(pos, pos[1:]):
            if b != a + 1:
                runs.append((s, a))
                s = b
        runs.append((s, pos[-1]))
        if len(runs) == 1:
            annotations[gene] = {"start": runs[0][0], "end": runs[0][1], "strand": "+", "type": "CDS"}
        else:   # spliced (M2, NS2) or frameshifted (PA-X) CDS: augur schema-annotations.json "segments"
            annotations[gene] = {"segments": [{"start": a, "end": b} for a, b in runs], "strand": "+", "type": "CDS"}
    spliced = [g for g, a in annotations.items() if "segments" in a]
    if spliced:
        parts = ", ".join(f"{g} ({len(annotations[g]['segments'])} parts)" for g in spliced)
        say("INFO", f"multi-part CDS written as annotation segments: {parts}.")
    aa_nodes = {n: {"aa_muts": {}} for n in nodes_nt}
    aa_ref, genes_done = {}, []
    for fa in sorted(glob.glob(os.path.join(ANC_AA, "*.fasta"))):
        gene = os.path.basename(fa)[: -len(".fasta")]
        if gene not in annotations:
            say("WARNING", f"ancestral_aa/{gene}.fasta has no coordinates in gene_map.tsv: its mutations are not exported.")
            continue
        aa = {r.id: str(r.seq).upper() for r in SeqIO.parse(fa, "fasta")}
        lacking = [n for n in nodes_nt if n not in aa]
        if lacking:
            raise Fail(f"ancestral_aa/{gene}.fasta lacks node(s) {', '.join(lacking[:10])}")
        ncod = len(gene_cols[gene]) // 3 + unmapped.get(gene, 0)
        if len(aa[refined.root.name]) != ncod:
            say("WARNING", f"{gene}: protein length {len(aa[refined.root.name])} != {ncod} codons in gene_map.tsv.")
        aa_nodes[refined.root.name]["aa_muts"][gene] = []
        for p, c in edges:
            aa_nodes[c]["aa_muts"][gene] = diffs(aa[p], aa[c], lambda s: s not in AA_SKIP)
        aa_ref[gene] = aa[refined.root.name]
        genes_done.append(gene)
    generated_by = {}
    try:
        from packaging.version import Version
        Version(FLUTYPER_VERSION)
        generated_by = {"generated_by": {"program": "FluTyper", "version": FLUTYPER_VERSION}}
    except Exception:
        say("WARNING", f"FluTyper version '{FLUTYPER_VERSION}' is not a PEP 440 version (augur export rejects it): "
                       "node-data written without generated_by.")
    json.dump({**generated_by, "nodes": nodes_nt, "reference": {"nuc": nt[refined.root.name]}},
              open("nt_muts.json", "w"), indent=1)
    json.dump({**generated_by, "annotations": annotations, "nodes": aa_nodes, "reference": aa_ref},
              open("aa_muts.json", "w"), indent=1)
    n_nt = sum(len(v["muts"]) for v in nodes_nt.values())
    n_aa = sum(len(m) for v in aa_nodes.values() for m in v["aa_muts"].values())
    say("INFO", f"branch mutations from FluTyper's parsimony reconstruction on the refined tree: {n_nt} nt, {n_aa} aa "
                f"over {len(edges)} branches; genes: {', '.join(genes_done) or 'none'}.")

    # ---- DTA: augur traits (D10, D13) --------------------------------------------------------------------------------
    node_data = ["branch_lengths.json", "nt_muts.json", "aa_muts.json"]
    wanted = [c.strip() for c in TRAITS.split(",") if c.strip()]
    usable = []
    for c in wanted:
        if c not in ANN_FIELDS + GEO_FIELDS:
            say("WARNING", f"trait '{c}' is not a metadata column of the dataset ({', '.join(ANN_FIELDS + GEO_FIELDS)}): skipped.")
        elif not any(values[t][c] is not None for t in tips):
            say("WARNING", f"trait '{c}': no tip has a value: no ancestral-state reconstruction.")
        else:
            usable.append(c)
    if usable:
        rc = run(["augur", "traits", "--tree", "refined.nwk", "--metadata", "meta.tsv", "--columns", *usable,
                  "--confidence", "--output-node-data", "traits_raw.json"], "traits.log")   # + mugration model files
        if rc != 0:
            say("WARNING", f"augur traits failed (exit {rc}: {tail('traits.log', 3)}): exported without inferred "
                           "ancestral states (no transmission lines on the map); tips keep their metadata.")
        else:
            d = json.load(open("traits_raw.json"))
            n_strip = 0
            for c in usable:
                for t in tips:
                    if values[t][c] is None:
                        node = d["nodes"].get(t, {})
                        for k in [k for k in node if k == c or k.startswith(c + "_")]:
                            del node[k]
                            n_strip += 1
            json.dump(d, open("traits.json", "w"), indent=1)
            node_data.append("traits.json")
            say("INFO", f"ancestral {', '.join(usable)} by DTA (augur traits --confidence; exploratory: sampling-biased); "
                        f"{n_strip} state(s) augur imputed to tips with a missing value were removed (never impute).")

    # ---- Bootstrap support (D18): hover attribute on every surviving internal node, label when >= threshold -------
    support = {}
    for r in read_table(SUPPORT):
        v = na(r.get("support"))
        if v is not None:
            support[r["node"].strip()] = v
    kept = {n: v for n, v in support.items() if n in out_clades}
    dropped = [n for n in support if n in removed]
    labels, odd = {}, []
    for n, v in kept.items():
        try:
            if float(v) >= MIN_SUPPORT:
                labels[n] = {"labels": {"support": v}}
        except ValueError:
            odd.append(f"{n}={v}")
    if odd:
        say("WARNING", f"{len(odd)} support value(s) are not numbers (shown on hover, never labelled): {', '.join(odd[:10])}")
    say("INFO", f"bootstrap support: {len(kept)} internal node(s) carry a value; {len(labels)} labelled (>= {MIN_SUPPORT:g}"
                f"{', --minSupportLabel 0 labels every value' if MIN_SUPPORT <= 0 else ''}); no branch is removed for low support.")
    if dropped:
        say("INFO", f"{len(dropped)} support value(s) belonged to zero-length branches collapsed by augur refine (D6): not shown.")
    if labels:     # no labels (e.g. --bootstrap 0): no node-data file, which augur would only warn about as empty
        json.dump({**generated_by, "nodes": {}, "branches": labels}, open("support.json", "w"), indent=1)
        node_data.append("support.json")
    # export-only metadata: meta.tsv + a 'support' column, with rows for the internal nodes (augur export matches
    # metadata rows by node name); refine and traits never see these rows
    with open("meta_export.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
        w.writerow(columns + ["support"])
        for t in tips:
            w.writerow([t, dates.get(t, "")] + [values[t][c] if values[t][c] is not None else MISSING
                                                for c in ANN_FIELDS + GEO_FIELDS] + [MISSING])
        for n, v in kept.items():
            w.writerow([n, ""] + [MISSING] * len(ANN_FIELDS + GEO_FIELDS) + [v])

    # ---- Dataset config: mode and ordinal legends --------------------------------------------------------------------
    cfg["title"] = f"{cfg.get('title', 'FluTyper')}: {label}"
    if kept and "support" not in cfg.get("metadata_columns", []):
        cfg["metadata_columns"] = cfg.get("metadata_columns", []) + ["support"]
    if mode != "timetree":
        cfg["colorings"] = [c for c in cfg.get("colorings", []) if c.get("key") != "num_date"]
        dd = cfg.setdefault("display_defaults", {})
        dd["distance_measure"] = "div"
        if dd.get("color_by") == "num_date":
            dd.pop("color_by")
        cfg["filters"] = [f for f in cfg.get("filters", []) if f != "num_date"]
    def age_key(v):
        m = re.fullmatch(r"\\s*(\\d+(?:\\.\\d+)?)\\s*(?:-\\s*\\d+(?:\\.\\d+)?|\\+)\\s*", v)
        if m:
            return (0, float(m.group(1)), v)
        m = re.fullmatch(r"\\s*(>=?|<=?)\\s*(\\d+(?:\\.\\d+)?)\\s*", v)
        if m:
            x = float(m.group(2))
            return (0, {">": x + 0.5, ">=": x, "<": -1.0, "<=": -1.0}[m.group(1)], v)
        return (1, 0.0, v)          # not an age range: after the ranges, alphabetically
    for c in cfg.get("colorings", []):
        if c.get("key") not in ("age_group", "season") or "scale" in c or "legend" in c:
            continue
        present = sorted({values[t][c["key"]] for t in tips if values[t].get(c["key"]) is not None},
                         key=age_key if c["key"] == "age_group" else (lambda v: v))
        if not present:
            continue
        if c["key"] == "age_group":
            odd = [v for v in present if age_key(v)[0] == 1]
            if odd:
                say("INFO", f"age_group value(s) {', '.join(odd)} are not age ranges: listed last in the legend.")
        n = len(present)
        idx = [round(i * (len(ORDINAL_COLOURS) - 1) / (n - 1)) if n > 1 else 0 for i in range(n)]
        c["scale"] = [[v, ORDINAL_COLOURS[i]] for v, i in zip(present, idx)]
        c["legend"] = [{"value": v} for v in present]
    json.dump(cfg, open("dataset_config.json", "w"), indent=1)

    # ---- augur export v2 + validate ----------------------------------------------------------------------------------
    cmd = ["augur", "export", "v2", "--tree", "refined.nwk", "--metadata", "meta_export.tsv", "--node-data", *node_data,
           "--auspice-config", "dataset_config.json", "--output", "unvalidated.json"]
    if LATLONGS:
        cmd[-2:-2] = ["--lat-longs", LATLONGS]
    rc = run(cmd, "export.log")
    if rc != 0:
        raise Fail(f"augur export v2 failed, exit {rc}: {tail('export.log')}")
    js = json.load(open("unvalidated.json"))
    exported = {c["key"] for c in js["meta"].get("colorings", [])}
    absent = {c["key"] for c in cfg.get("colorings", []) if c["key"] not in exported}   # WARNed once, below
    for line in open("export.log", errors="replace"):
        m = re.match(r"\\s*(WARNING|ERROR):?\\s*(.*)", line)
        if not m:
            continue
        # augur's two messages about a value without coordinates duplicate our own WARNING above: skip exactly those
        dup = re.search(r"(\\w+)->(.+?) did not have an associated lat/long", m.group(2))
        if dup and (dup.group(1).lower(), dup.group(2).strip().lower()) in no_coords:
            continue
        dup = re.search(r'"(.+?)", a value of the geographic resolution "(\\w+)", appears in the tree but not in the metadata', m.group(2))
        if dup and (dup.group(2).lower(), dup.group(1).strip().lower()) in no_coords:
            continue
        # augur's messages about a colouring/filter with no value duplicate our "colouring ... is not in the dataset"
        dup = (re.search(r"color-by field '(\\w+)' does not exist", m.group(2))
               or re.search(r'filter "(\\w+)" does not appear as a property', m.group(2)))
        if dup and dup.group(1) in absent:
            continue
        say("WARNING", f"augur export: {m.group(2).strip()}")
    rc = run(["augur", "validate", "export-v2", "unvalidated.json"], "validate.log")
    if rc != 0:
        raise Fail(f"augur validate export-v2 rejected the dataset: {tail('validate.log')}")
    out_tips = []
    def walk(n):
        if n.get("children"):
            for ch in n["children"]:
                walk(ch)
        else:
            out_tips.append(n["name"])
    walk(js["tree"])
    if sorted(out_tips) != sorted(tips):
        raise Fail(f"the exported tree's tips differ from the input tree's ({len(out_tips)} vs {len(tips)})")
    for c in cfg.get("colorings", []):
        if c["key"] in absent:
            say("WARNING", f"colouring '{c['key']}' is not in the dataset (no node has a value).")
    os.replace("unvalidated.json", dataset)
    say("INFO", f"augur validate export-v2: {dataset} is valid ({len(out_tips)} tips, "
                f"{'time tree' if mode == 'timetree' else 'divergence-only tree'}).")
    status.update(status="ok", reason="")


try:
    main()
except Fail as e:
    say("ERROR", f"{e} No Auspice dataset for {TREE_ID}.")
    status["reason"] = str(e)
except Exception as e:
    say("ERROR", f"unexpected failure: {type(e).__name__}: {e}. No Auspice dataset for {TREE_ID}.")
    status["reason"] = f"unexpected {type(e).__name__}: {e}"
    traceback.print_exc()
finally:
    if status["status"] != "ok":
        for n in glob.glob("flutyper_*.json"):
            os.remove(n)       # never emit a dataset for a failed tree
        sys.stderr.write("\\n".join(l for l in log_lines if " ERROR " in l) + "\\n")
    with open("auspice.log", "w", encoding="utf-8") as f:
        f.write("\\n".join(log_lines) + "\\n")
    with open("auspice_status.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\\t", lineterminator="\\n")
        w.writerow(["tree_id", "status", "mode", "dataset", "reason"])
        w.writerow([TREE_ID, status["status"], status["mode"], status["dataset"],
                    re.sub(r"\\s+", " ", status["reason"]) or "NA"])
PYEOF
    """
}
