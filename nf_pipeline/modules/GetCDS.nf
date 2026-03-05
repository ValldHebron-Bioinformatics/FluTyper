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
    cp -rL ${sequences_dir} sequences_work
    rm -f sequences
    mv sequences_work sequences

    export SEQUENCES_DIR="sequences"
    export INFERRED_SUBTYPES="${inferred_subtypes}"
    export REFERENCES_FASTA="${params.protocols[params.protocol].resources}/CDS_references.fasta"

    python3 - <<'PY'
import csv
import os
import re
import subprocess
import tempfile
from pathlib import Path

def parse_fasta(source):
    '''Yield (header, seq) from a file path or a string of FASTA text.'''
    if isinstance(source, Path):
        with open(source, "r") as fh:
            lines = fh.readlines()
    else:
        lines = source.splitlines()

    header, seq = None, []
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq).upper()
            header, seq = line[1:], []
        else:
            seq.append(line)
    if header is not None:
        yield header, "".join(seq).upper()


def write_fasta(path, records):
    with open(path, "w") as fh:
        for h, s in records:
            print(">" + h, file=fh)
            print(s, file=fh)

# Ref map

def get_ref_map(ref_fasta):
    '''Build refs[subtype][pathotype][variant] = seq from CDS_references.fasta.'''
    refs = {}
    for header, seq in parse_fasta(Path(ref_fasta)):
        parts = header.split("_")
        if len(parts) < 3:
            continue
        subtype = parts[0]
        variant = parts[1]
        pathotype = parts[-1]
        if not re.fullmatch(r"H[0-9]+N[0-9]+", subtype):
            continue
        if pathotype not in ("HPAI", "LPAI"):
            continue
        refs.setdefault(subtype, {}).setdefault(pathotype, {})[variant] = seq
    return refs


def load_subtype_map(tsv_path):
    '''Return {sample_id: inferred_subtype} from the TSV.'''
    out = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter=chr(9))
        sample_col = "seqName" if "seqName" in reader.fieldnames else "sample"
        subtype_col = "inferred_subtype" if "inferred_subtype" in reader.fieldnames else "subtype"
        for row in reader:
            sid = (row.get(sample_col) or "").strip()
            if sid:
                out[sid] = (row.get(subtype_col) or "").strip()
    return out


_FAMILY = {
    r"H5N[0-9]+": ("H5N1", "HPAI"),
    r"H7N[0-9]+": ("H7N9", "HPAI"),
    r"H9N[0-9]+": ("H9N2", "LPAI"),
}

def choose_family(inferred):
    for pattern, family in _FAMILY.items():
        if re.fullmatch(pattern, inferred or ""):
            return family
    return "H5N1", "HPAI"          # fallback


VARIANT_MAP = {
    "HA":  ["HA", "HA(-SP)", "HA1", "HA1(-SP)", "HA2"],
    "PA":  ["PA", "PA-X"],
    "PB1": ["PB1", "PB1-F2"],
    "PB2": ["PB2"],
    "NA":  ["NA"],
    "NP":  ["NP"],
    "NS":  ["NS1", "NEP"],
    "MP":  ["M1", "M2"],
}


def choose_ref_seq(refs, canonical, preferred_pathotype, variant):
    '''Return the best matching reference sequence, falling back broadly.'''
    alt = "LPAI" if preferred_pathotype == "HPAI" else "HPAI"

    for subtype in [canonical] + [s for s in refs if s != canonical]:
        for pathotype in [preferred_pathotype, alt]:
            seq = refs.get(subtype, {}).get(pathotype, {}).get(variant)
            if seq:
                return seq
    return ""

def extract_cds(query_seq, ref_seq):
    '''Align query to ref with MAFFT and return the CDS-trimmed query sequence.'''
    with tempfile.TemporaryDirectory() as tdir:
        pair_fa = Path(tdir) / "pair.fa"
        write_fasta(pair_fa, [("REF", ref_seq), ("QRY", query_seq)])

        proc = subprocess.run(
            ["mafft", "--globalpair", "--op", "3.0", "--ep", "0.5", "--quiet", str(pair_fa)],
            capture_output=True, text=True, check=False,
        )
        if proc.returncode != 0:
            return ""

        aln = dict(parse_fasta(proc.stdout))
        ref_aln, qry_aln = aln.get("REF", ""), aln.get("QRY", "")
        if not ref_aln or not qry_aln:
            return ""

        cds = "".join(q for r, q in zip(ref_aln, qry_aln) if r != "-")
        cds = cds.replace("-", "").upper()
        cds = re.sub(r"[^ACGTN]", "", cds)

        usable = (len(cds) // 3) * 3
        return cds[:usable] if usable > 0 else ""

def main():
    sequences_dir = Path(os.environ["SEQUENCES_DIR"])
    inferred_tsv = Path(os.environ["INFERRED_SUBTYPES"])
    references_fa = Path(os.environ["REFERENCES_FASTA"])

    refs = get_ref_map(references_fa)
    subtype_map = load_subtype_map(inferred_tsv)

    for sample_dir in filter(Path.is_dir, sequences_dir.iterdir()):
        seg_dir = sample_dir / "segments"
        if not seg_dir.is_dir():
            continue

        sample_id = sample_dir.name
        inferred = subtype_map.get(sample_id, "")
        canonical, pref = choose_family(inferred)

        if not re.fullmatch(r"H(5|7|9)N[0-9]+", inferred or ""):
            print(f"[GetCDS] WARNING: Unexpected subtype for {sample_id} ({inferred}). "
                  f"Falling back to H5N1_*_HPAI", flush=True)

        cds_dir = seg_dir / "CDS"
        cds_dir.mkdir(parents=True, exist_ok=True)

        for seg_fa in seg_dir.glob("*.fasta"):
            seg_name = seg_fa.stem
            if seg_name.endswith("_CDS"):
                continue

            variants = VARIANT_MAP.get(seg_name, [])
            query_records = list(parse_fasta(seg_fa))
            if not variants or not query_records:
                continue

            for variant in variants:
                ref_seq = choose_ref_seq(refs, canonical, pref, variant)
                out_path = cds_dir / f"{variant}_CDS.fasta"

                if not ref_seq:
                    out_path.unlink(missing_ok=True)
                    continue

                results = [
                    (f"{qh}|{variant}", extract_cds(qs, ref_seq))
                    for qh, qs in query_records
                ]
                results = [(h, s) for h, s in results if s]

                if results:
                    write_fasta(out_path, results)
                else:
                    out_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()

PY

    """
}


