#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsFinder {
    errorStrategy 'ignore'
    debug true

    input:
    tuple val(sample_id), path(prot_files), val(h_tag), val(n_tag), val(pathotype)

    output:
    tuple val(sample_id), path("samples/${sample_id}/mutations/${sample_id}_*_mutations.csv"), path("samples/${sample_id}/${sample_id}_mutations.csv"), emit: results
    tuple val(sample_id), path("MFerrors.log"), optional: true, emit: errors

    script:
    """#!/usr/bin/env python3
import os, csv
from pathlib import Path
from Bio import SeqIO

markers_dir = Path("${params.protocols[params.protocol].resources}/markers")
ha_dict = "${params.protocols[params.protocol].resources}/HA_DICT.csv"
na_dict = "${params.protocols[params.protocol].resources}/NA_DICT.csv"
out_dir = Path("samples/${sample_id}/mutations")
out_dir.mkdir(parents=True, exist_ok=True)

target_H = "${h_tag}" if "${h_tag}".startswith("H") and "${h_tag}"[1:].isdigit() else "H5"
target_N = "${n_tag}" if "${n_tag}".startswith("N") and "${n_tag}"[1:].isdigit() else "N1"

def build_pos_lookup(dict_path, from_col, to_col, prot_filter=None):
    if not os.path.exists(dict_path):
        return {}
    with open(dict_path, 'r') as f:
        reader = csv.DictReader(f) 
        # Clean up header names
        headers = [h.strip() for h in (reader.fieldnames or [])]
        # Validation
        if from_col not in headers or to_col not in headers:
            return {}
        # Build the mapping with additional filtering for the target protein
        return {
            row[from_col].strip(): row[to_col].strip() 
            for row in reader
            if row.get(from_col, '').strip() and row.get(to_col, '').strip()
            and (not prot_filter or row.get('PROTEIN', '').strip().upper() == prot_filter.upper())
        }

def get_unique(seq):
    already_recorded = set()
    return [x for x in seq if not (x in already_recorded or already_recorded.add(x))]

output_files = []
for aligned_prot in "${prot_files}".split():
    file_path = Path(aligned_prot)
    prot_name = file_path.name.replace("_PROT_aligned.fasta", "").split('_')[1]

    records = list(SeqIO.parse(file_path, "fasta"))
    if len(records) < 2:
        with open("MFerrors.log", 'a') as f: f.write(f"ERROR: Less than 2 sequences in ${sample_id} {prot_name}, error in the alignment\\n")
        continue

    ref_seq, query_seq = str(records[0].seq), str(records[1].seq)
    ref_tag = str(records[0].description).split('_')[0]
    subtype_val = "${h_tag}${n_tag}(${pathotype})" if "${pathotype}" else "${h_tag}${n_tag}"

    pos_to_target, pos_to_base = {}, {}
    if prot_name.startswith("HA"):
        ref_H = ref_tag if ref_tag.startswith("H") and ref_tag[1:].isdigit() else "H5"
        # Build position lookups for both target and reference numbering systems
        pos_to_target = build_pos_lookup(ha_dict, f"{ref_H}_numbering", f"{target_H}_numbering", prot_name)
        pos_to_base = build_pos_lookup(ha_dict, f"{ref_H}_numbering", "H5_numbering", prot_name)
    elif prot_name.startswith("NA"):
        ref_N = ref_tag if ref_tag.startswith("N") and ref_tag[1:].isdigit() else "N1"
        # Build position lookups for both target and reference numbering systems
        pos_to_target = build_pos_lookup(na_dict, f"{ref_N}_pos", f"{target_N}_pos")
        pos_to_base = build_pos_lookup(na_dict, f"{ref_N}_pos", "N1_pos")

    markers_by_id, marker_info = {}, {}
    m_file = markers_dir / f"{prot_name}_markers.csv"
    if m_file.exists():
        with open(m_file) as f:
            for row in csv.DictReader(f):
                m_id = row['MARKER_ID'].strip()
                # Store the mutation combination for each marker ID, setdefault is used to handle multiple entries for the same marker ID without overwriting
                markers_by_id.setdefault(m_id, set()).add((row['POSITION'].strip(), row['AA'].strip()))
                
                info = (row.get('EFFECT', '').strip(), row.get('FOUND_IN', '').strip(), row.get('REFERENCE', '').strip())
                if info not in marker_info.setdefault(m_id, []): marker_info[m_id].append(info)

    observed_mutations, protein_pos = set(), 0
    for r_aa, q_aa in zip(ref_seq, query_seq):
        if r_aa != "-": protein_pos += 1
        # Translate the current position into the standardized base numbering (like H5 or N1)
        standard_pos = pos_to_base.get(str(protein_pos), str(protein_pos))    
        # Record the combination of the standardized position and the amino acid found
        observed_mutations.add((standard_pos, q_aa))

    active_markers = {}
    for m_id, mut_set in markers_by_id.items():
        # Check if all the mutations in the marker's combination are present in the observed mutations
        if mut_set.issubset(observed_mutations):
            # If the marker is active, associate it with all its mutations for later reference in the output
            for mut in mut_set: active_markers.setdefault(mut, []).append(m_id)

    out_csv = out_dir / f"${sample_id}_{prot_name}_mutations.csv"
    output_files.append(out_csv)
    
    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        # Added FOUND_IN to the master header
        writer.writerow(["SAMPLE_ID", "SUBTYPE", "PROTEIN", "REF_SUBTYPE", "POSITION", "POSITION_REF", "REFERENCE_AA", "QUERY_AA", "AA_MUTATION", "MUTATION_TYPE", "MARKER", "MARKER_ID", "IS_COMBINATION", "EFFECT", "FOUND_IN", "REFERENCE"])
        
        pos = 0
        for r_aa, q_aa in zip(ref_seq, query_seq):
            if r_aa != "-": pos += 1
            pos_raw = str(pos)
            # Translate the raw position to the target numbering system and also get the standardized base position for reference
            # If no translation is found, it defaults to the raw position
            position, pos_ref = pos_to_target.get(pos_raw, pos_raw), pos_to_base.get(pos_raw, pos_raw)

            is_marker, is_mutation = (pos_ref, q_aa) in active_markers, r_aa != q_aa and "X" not in (r_aa, q_aa)
            if not (is_marker or is_mutation): continue

            mut_type = "Marker" if not is_mutation else "Insertion" if r_aa == "-" else "Deletion" if q_aa == "-" else "Substitution"
            aa_mut = f"{position}{q_aa}" if is_marker else f"{r_aa}{position}{q_aa}"

            if is_marker:
                m_ids = get_unique(active_markers[(pos_ref, q_aa)])
                is_combo = " | ".join(get_unique("Yes" if len(set(m[0] for m in markers_by_id[mid])) > 1 else "No" for mid in m_ids))
                
                # Extract and join effect, found_in, and reference cleanly
                effect = " | ".join(get_unique(eff for mid in m_ids for eff, _, _ in marker_info[mid] if eff))
                found = " | ".join(get_unique(fnd for mid in m_ids for _, fnd, _ in marker_info[mid] if fnd))
                ref = " | ".join(get_unique(ref for mid in m_ids for _, _, ref in marker_info[mid] if ref))
                
                writer.writerow(["${sample_id}", subtype_val, prot_name, ref_tag, position, pos_ref, r_aa, q_aa, aa_mut, mut_type, "Yes", " | ".join(m_ids), is_combo, effect, found, ref])
            else:
                writer.writerow(["${sample_id}", subtype_val, prot_name, ref_tag, position, pos_ref, r_aa, q_aa, aa_mut, mut_type, "No", "", "No", "", "", ""])

master_csv = Path(f"samples/${sample_id}/${sample_id}_mutations.csv")
if output_files:
    with open(master_csv, 'w', newline='') as mf:
        writer = csv.writer(mf)
        # Added FOUND_IN to the master compiler header
        writer.writerow(["SAMPLE_ID", "SUBTYPE", "PROTEIN", "REF_SUBTYPE", "POSITION", "POSITION_REF", "REFERENCE_AA", "QUERY_AA", "AA_MUTATION", "MUTATION_TYPE", "MARKER", "MARKER_ID", "IS_COMBINATION", "EFFECT", "FOUND_IN", "REFERENCE"])
        for f_path in output_files:
            with open(f_path, 'r') as f:
                next(f)
                writer.writerows(csv.reader(f))
"""
}