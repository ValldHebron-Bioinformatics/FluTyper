process MutationsFinder {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(prot_files), val(h_tag), val(n_tag), val(pathotype)

    output:
    tuple val(sample_id), path("samples/${sample_id}/mutations/*_${sample_id}_mutations.csv"), path("samples/${sample_id}/${sample_id}_mutations.csv"), emit: results
    tuple val(sample_id), path("MFerrors.log"), optional: true, emit: errors

    script:
    """#!/usr/bin/env python3
import csv, re, subprocess
from pathlib import Path
from Bio import SeqIO

dict_path = "${params.protocols[params.protocol].resources}/AA_DICT.csv"
ref_prot_path = "${params.protocols[params.protocol].resources}/PROT_references.fasta"
markers_dir = Path("${params.protocols[params.protocol].resources}/markers")
mutations_prog = "${params.programs.MutationsDictionary}"
target_h = "${h_tag}" if "${h_tag}" in {"H1", "H3", "H5", "H7", "H9"} else "H5"

# Force creation of the mutations directory so Nextflow publishes it
mut_dir = Path("samples/${sample_id}/mutations")
mut_dir.mkdir(parents=True, exist_ok=True)

refs = {rec.id: str(rec.seq) for rec in SeqIO.parse(ref_prot_path, "fasta")}
output_files = []

files = "${prot_files.join(' ')}".split()

for file in files:
    file_path = Path(file)
    if not file_path.exists(): continue
    prot_name = file_path.name.split('_')[1]

    # Find the matching reference sequence dynamically
    req_tag = "${n_tag}" if prot_name == "NA" else "${h_tag}"
    base_pattern = f"^{req_tag}_{prot_name}_"
    matching_refs = [r_id for r_id in refs.keys() if re.search(base_pattern, r_id)]
    
    ref_id_to_use = None
    if matching_refs:
        if "${pathotype}":
            patho_matches = [r for r in matching_refs if "${pathotype}" in r]
            ref_id_to_use = patho_matches[0] if patho_matches else matching_refs[0]
        else:
            ref_id_to_use = matching_refs[0]
            
    if not ref_id_to_use and not prot_name.startswith("HA") and prot_name != "NA":
        fallback_pattern = f"^H5_{prot_name}_"
        fallback_refs = [r_id for r_id in refs.keys() if re.search(fallback_pattern, r_id)]
        if fallback_refs:
            patho_matches = [r for r in fallback_refs if "HPAI" in r]
            ref_id_to_use = patho_matches[0] if patho_matches else fallback_refs[0]
            req_tag = "H5"
            
    if not ref_id_to_use: continue
    ref_seq = refs[ref_id_to_use]

    # Load known markers
    markers = {}
    m_file = markers_dir / f"{prot_name}_markers.csv"
    if m_file.exists():
        with open(m_file) as f:
            for r in csv.DictReader(f):
                markers[(r['POSITION'].strip(), r['AA'].strip())] = (r.get('ORIGIN',''), r['EFFECT'], r.get('REFERENCE',''))

    # Set the precise naming convention you requested
    output_csv = mut_dir / f"{prot_name}_${sample_id}_mutations.csv"
    output_files.append(output_csv)
    query_seq = str(next(SeqIO.parse(file_path, "fasta")).seq)

    with open(output_csv, 'w', newline='') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(["SAMPLE_ID", "SUBTYPE", "PROTEIN", "REF_SUBTYPE", "POSITION", "REFERENCE_AA", "QUERY_AA", "MARKER", "ORIGIN", "EFFECT", "REFERENCE"])
        
        # Direct comparison with strict positive numbering starting at 1
        pos = 0
        for r_aa, q_aa in zip(ref_seq, query_seq):
            if r_aa != '-': 
                pos += 1
                
            pos_str = str(pos)
            is_marker = (pos_str, q_aa) in markers
            is_mutation = r_aa != q_aa
            
            if is_marker:
                m_info = markers[(pos_str, q_aa)]
                writer.writerow(["${sample_id}", "${h_tag}${n_tag}", prot_name, req_tag, pos_str, r_aa, q_aa, "TRUE", m_info[0], m_info[1], m_info[2]])
            elif is_mutation:
                writer.writerow(["${sample_id}", "${h_tag}${n_tag}", prot_name, req_tag, pos_str, r_aa, q_aa, "FALSE", "", "", ""])

    # Delegate external translation for non-H5 HA sequences
    if prot_name.startswith("HA") and target_h != "H5":
        subprocess.run(["python3", mutations_prog, "--subtype", target_h, "--input", str(output_csv), "--dictionary", dict_path, "--output", str(output_csv)], check=True)

# Compile everything into the master CSV
master_csv = Path(f"samples/${sample_id}/${sample_id}_mutations.csv")
if output_files:
    with open(master_csv, 'w') as master_f:
        for i, f_path in enumerate(output_files):
            with open(f_path) as f:
                if i > 0: next(f)
                master_f.write(f.read())
"""
}