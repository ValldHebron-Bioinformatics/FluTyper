#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GetCDS {
    errorStrategy 'ignore'

    input:
    path(sequences_dir)
    path(inferred_subtypes)

    output:
    path("sequences")

    script:
    """
    cp -rL ${sequences_dir} sequences_work && rm -rf sequences && mv sequences_work sequences
    export REFERENCES_FASTA="${params.protocols[params.protocol].resources}/CDS_references.fasta"
    export INFERRED_SUBTYPES="${inferred_subtypes}"
    export SEQUENCE_DIR="sequences"

    python3 - <<'PY'
import csv, os, re, subprocess, tempfile
from pathlib import Path

def parse_fasta(p):
    if not Path(p).exists():
        return []
    out = []
    with open(p) as f:
        for g in [g.splitlines() for g in f.read().split(">") if g.strip()]:
            h = g[0].strip()
            s = "".join(x.strip() for x in g[1:] if x.strip()).upper()
            out.append((h, s))
    return out

# Call mafft to align ref and query
def run_mafft(ref, qry):
    key = (ref, qry)
    if key in run_mafft._cache:
        return run_mafft._cache[key]
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa') as t:
        t.write(f">R\\n{ref}\\n>Q\\n{qry}\\n")
        t.flush()
        res = subprocess.run(
            ["mafft", "--auto", "--op", "3.0", "--quiet", t.name],
            capture_output=True, text=True
        ) 
        if res.returncode != 0:
            run_mafft._cache[key] = ("", "")
            return run_mafft._cache[key]
        d = {l.strip(): "".join(lines).strip()
             for l, *lines in [g.splitlines() for g in res.stdout.split(">") if g.strip()]}
        run_mafft._cache[key] = (d.get("R", ""), d.get("Q", ""))
        return run_mafft._cache[key]
run_mafft._cache = {} 

# Extract CDS from aligned sequences, trimming leading/trailing gaps and ensuring length is a multiple of 3
def get_cds(r_aln, q_aln):
    if not r_aln or not q_aln:
        return ""
    s, e = 0, len(r_aln) - 1
    while s <= e and (r_aln[s] == '-' or q_aln[s] == '-'):
        s += 1
    while e >= s and (r_aln[e] == '-' or q_aln[e] == '-'):
        e -= 1
    if s > e:
        return ""
    cds = "".join(q for r, q in zip(r_aln[s:e+1], q_aln[s:e+1]) if r != '-')
    ungapped = sum(1 for c in cds if c != '-')
    usable = (ungapped // 3) * 3
    if usable <= 0:
        return ""
    out, seen = [], 0
    for c in cds:
        if c != '-':
            if seen >= usable:
                break
            seen += 1
        out.append(c)
    return "".join(out)
# Score alignment by counting gaps and mismatches (gaps worse than mismatches)
def aln_score(r_aln, q_aln):
    if not r_aln or not q_aln:
        return (10**12, -10**12) # Worst possible score for failed alignments
    gaps = q_aln.count('-') + r_aln.count('-')
    matches = sum(1 for r, q in zip(r_aln, q_aln) if r == q and r != '-')
    return (gaps, -matches)

# Among candidates, pick the one with the best alignment score to the query
def best_by_alignment(cands, query_seq):
    best, best_score = None, (10**12, -10**12)
    for c in cands:
        r_aln, q_aln = run_mafft(c['s'], query_seq)
        sc = aln_score(r_aln, q_aln)
        if sc < best_score:
            best, best_score = c, sc
    return best
# Choose reference sequence based on inferred subtype, segment, and variant, with fallbacks
def choose_ref(refs, inferred_subtype, seg, var, query_seq):
    h = (re.search(r"H\d+", inferred_subtype) or [None])[0]
    n = (re.search(r"N\d+", inferred_subtype) or [None])[0]
    # For H5/H7, prioritize HPAI/LPAI matches; for other segments, just pick the best by alignment among subtype/variant matches
    def hp_lp_filter(cands):
        if h in ("H5", "H7"):
            filtered = [r for r in cands if r['p'] in ("HPAI", "LPAI")]
            if filtered:
                return filtered
        return cands

    if seg == "HA" and h:
        # strict H-match for HA segment (allow HA fallback for derived HA variants)
        cands = hp_lp_filter([r for r in refs if r['h'] == h and r['v'] in ("HA", var)])
        best = best_by_alignment(cands, query_seq) if cands else None
        if best and best['h'] == h:
            return best['s'], best['id']

        # If no H-match, default fallback to H5N1 HPAI
        fb = [r for r in refs if r['sub'] == "H5N1" and r['p'] == "HPAI" and r['v'] in ("HA", var)]
        if fb:
            return fb[0]['s'], fb[0]['id']
        return None, "None"

    if seg == "NA" and n:
        # strict N-match for NA segment
        cands = hp_lp_filter([r for r in refs if r['v'] == "NA" and r['n'] == n])
        best = best_by_alignment(cands, query_seq) if cands else None
        if best and best['n'] == n:
            return best['s'], best['id']

        # If no N-match, default fallback to H5N1 HPAI
        fb = [r for r in refs if r['sub'] == "H5N1" and r['p'] == "HPAI" and r['v'] == "NA"]
        if fb:
            return fb[0]['s'], fb[0]['id']
        return None, "None"

    # Remaining segments: pick by protein variant
    if h in ("H5", "H7", "H9"):
        cands = hp_lp_filter([r for r in refs if r['v'] == var and r['h'] == h])
        best = best_by_alignment(cands, query_seq) if cands else None
        if best:
            return best['s'], best['id']

    # Final fallback: H5N1 HPAI by protein variant
    fallback = [r for r in refs if r['sub'] == "H5N1" and r['v'] == var and r['p'] == "HPAI"]
    if fallback:
        return fallback[0]['s'], fallback[0]['id']
    return None, "None"

# Main
refs_raw = parse_fasta(os.environ["REFERENCES_FASTA"])
refs = []
for h, s in refs_raw:
    p = h.split("_")
    if len(p) < 2:
        continue
    refs.append({'id': h, 'sub': p[0].upper(), 'v': p[1], 'p': p[-1].upper(),
                 'h': (re.search(r"H\d+", p[0].upper()) or [None])[0],
                 'n': (re.search(r"N\d+", p[0].upper()) or [None])[0], 's': s})

def load_subtypes(tsv_path):
    out = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not reader.fieldnames:
            return out
        sample_col = "seqName" if "seqName" in reader.fieldnames else "sample"
        subtype_col = "inferred_subtype" if "inferred_subtype" in reader.fieldnames else "subtype"
        for row in reader:
            sid = (row.get(sample_col) or "").strip()
            st  = (row.get(subtype_col) or "").strip().upper()
            if sid:
                out[sid] = st
    return out

subtypes = load_subtypes(os.environ["INFERRED_SUBTYPES"])
report = []

v_map = {
    "HA":  ["HA", "HA(-SP)", "HA1(-SP)", "HA2"],
    "PA":  ["PA", "PA-X"],
    "PB1": ["PB1", "PB1-F2"],
    "PB2": ["PB2"],
    "NA":  ["NA"],
    "NP":  ["NP"],
    "NS":  ["NS1", "NEP"],
    "MP":  ["M1", "M2"],
}

for s_dir in filter(Path.is_dir, Path("sequences").iterdir()):
    sub = subtypes.get(s_dir.name, "H5N1").upper()
    seg_dir = s_dir / "segments"
    if not seg_dir.exists():
        continue

    (seg_dir / "CDS").mkdir(exist_ok=True)

    for f in seg_dir.glob("*.fasta"):
        if "_CDS" in f.name:
            continue
        q_records = parse_fasta(f)
        if not q_records:
            continue

        variants = v_map.get(f.stem, [])
        if not variants:
            continue

        query0 = q_records[0][1]

        # HA/NA are segment-driven; remaining segments are variant-driven
        seg_ref = None
        if f.stem in ("HA", "NA"):
            seg_ref = choose_ref(refs, sub, f.stem, None, query0)

        for var in variants:
            if seg_ref is not None:
                r_s, r_id = seg_ref
            else:
                r_s, r_id = choose_ref(refs, sub, f.stem, var, query0)

            if r_s:
                res = [(f"{qh}|{var}", get_cds(*run_mafft(r_s, qs)))
                       for qh, qs in q_records]
                res = [(h, s) for h, s in res if s]
                if res:
                    with open(seg_dir / "CDS" / f"{var}_CDS.fasta", "w") as out:
                        for h, s in res:
                            out.write(f">{h}\\n{s}\\n")
                report.append({
                    "sample": s_dir.name, "segment": f.stem, "cds_variant": var,
                    "inferred_subtype": sub, "ref": r_id,
                    "status": "Success" if res else "Failed"
                })
            else:
                report.append({
                    "sample": s_dir.name, "segment": f.stem, "cds_variant": var,
                    "inferred_subtype": sub, "ref": "None", "status": "Failed-NoRef"
                })

with open("sequences/cds_report.csv", "w") as f:
    w = csv.DictWriter(f, fieldnames=["sample", "segment", "cds_variant",
                                      "inferred_subtype", "ref", "status"])
    w.writeheader()
    w.writerows(report)
PY
    """
}