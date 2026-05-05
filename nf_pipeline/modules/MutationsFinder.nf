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
    with open(dict_path, 'r', encoding='utf-8-sig') as f:
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
            if row.get(from_col, '').strip() 
            and row.get(to_col, '').strip()
            and (not prot_filter or row.get('PROTEIN', '').strip().upper() == prot_filter.upper())
        }

output_files = []
for aligned_prot in "${prot_files}".split():
    file_path = Path(aligned_prot)
    prot_name = file_path.name.split('_')[1]

    records = list(SeqIO.parse(file_path, "fasta"))
    if len(records) < 2:
        with open("MFerrors.log", 'a') as f: f.write(f"ERROR: Less than 2 sequences in ${sample_id} {prot_name}, error in the alignment\\n")
        continue

    ref_seq, query_seq = str(records[0].seq), str(records[1].seq)
    ref_tag = str(records[0].description).split('_')[0]
    subtype_val = "${h_tag}${n_tag}(${pathotype})" if "${pathotype}" else "${h_tag}${n_tag}"

    pos_to_base = {}
    if "${params.protocol}" == "AVIAN":
        if prot_name.startswith("HA"):
            ref_H = ref_tag if ref_tag.startswith("H") and ref_tag[1:].isdigit() else "H5"
            pos_to_base = build_pos_lookup(ha_dict, f"{ref_H}_numbering", "H5_numbering", prot_name)
        elif prot_name.startswith("NA"):
            ref_N = ref_tag if ref_tag.startswith("N") and ref_tag[1:].isdigit() else "N1"
            pos_to_base = build_pos_lookup(na_dict, f"{ref_N}_pos", "N1_pos")
    elif "${params.protocol}" == "HUMAN":
        if prot_name.startswith("HA"):
            ref_H = "${h_tag}" if "${h_tag}".startswith("H") else "H1"
            pos_to_base = build_pos_lookup(ha_dict, f"{ref_H}_numbering", "H1_numbering", prot_name)
        elif prot_name.startswith("NA"):
            ref_N = "${n_tag}" if "${n_tag}".startswith("N") else "N1"
            pos_to_base = build_pos_lookup(na_dict, f"{ref_N}_pos", "N1_pos")

    markers_by_id, marker_info = {}, {}
    m_file = markers_dir / f"{prot_name}_markers.csv"
    if m_file.exists():
        with open(m_file, encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                # Enforce subtype matching for HUMAN protocol
                if "${params.protocol}" == "HUMAN" and "${h_tag}${n_tag}" not in row.get('FOUND_IN', ''):
                    continue
                    
                m_id = row['MARKER_ID'].strip()
                # setdefault allows us to merge multiple position/AA pairs for the same marker ID
                markers_by_id.setdefault(m_id, set()).add((row['POSITION'].strip(), row['AA'].strip()))
                
                info = (row.get('EFFECT', '').strip(), row.get('FOUND_IN', '').strip(), row.get('REFERENCE', '').strip())
                # Avoid duplicate info entries for the same marker ID
                if info not in marker_info.setdefault(m_id, []): marker_info[m_id].append(info)

    # Add "X" mutations for all observed positions to allow marker detection ONLY if it is an actual mutation
    observed_mutations, protein_pos = set(), 0
    for r_aa, q_aa in zip(ref_seq, query_seq):
        if r_aa != "-": protein_pos += 1
        standard_pos = pos_to_base.get(str(protein_pos), str(protein_pos))
        observed_mutations.add((standard_pos, q_aa))
        
        # X acts as a wildcard only if the amino acid changed
        if r_aa != q_aa and "X" not in (r_aa, q_aa):
            observed_mutations.add((standard_pos, "X"))
            
    # Check which observed mutations are part of any marker sets
    active_markers = {}
    for m_id, mut_set in markers_by_id.items():
        if mut_set.issubset(observed_mutations):
            # Mutations can have multiple marker ids, so we append to a list for each mutation
            for mut in mut_set: active_markers.setdefault(mut, []).append(m_id)
            
    # Build list to store indels
    events = []
    pos = 0
    for r_aa, q_aa in zip(ref_seq, query_seq):
        if r_aa != "-": pos += 1
        pos_raw = str(pos)
        pos_ref = pos_to_base.get(pos_raw, pos_raw)

        is_mutation = r_aa != q_aa and "X" not in (r_aa, q_aa)

        # Look first for exact matches, then for "X" matches only if it is a true mutation
        m_ids_exact = active_markers.get((pos_ref, q_aa), [])
        m_ids_x = active_markers.get((pos_ref, "X"), []) if is_mutation else []
        combined_m_ids = list(dict.fromkeys(m_ids_exact + m_ids_x)) 
        
        is_marker = len(combined_m_ids) > 0
        
        if not (is_marker or is_mutation):
            continue

        if is_marker:
            m_ids = combined_m_ids
            is_combo = " | ".join(dict.fromkeys("Yes" if len(set(m[0] for m in markers_by_id[mid])) > 1 else "No" for mid in m_ids))
            # Clean and deduplicate by splitting on pipes inline
            effect = " | ".join(dict.fromkeys(item.strip() for mid in m_ids for eff, _, _ in marker_info[mid] if eff for item in eff.split('|') if item.strip()))
            found  = " | ".join(dict.fromkeys(item.strip() for mid in m_ids for _, fnd, _ in marker_info[mid] if fnd for item in fnd.split('|') if item.strip()))
            ref    = " | ".join(dict.fromkeys(item.strip() for mid in m_ids for _, _, r   in marker_info[mid] if r for item in r.split('|') if item.strip()))
            mut_type = "Marker"
            # Included the reference amino acid conditionally if a mutation actually occurred
            aa_mut = f"{r_aa}{pos_raw}{q_aa}" if is_mutation else f"{pos_raw}{q_aa}"
            
            events.append({"pos": pos_raw, "pos_ref": pos_ref, "r_aa": r_aa, "q_aa": q_aa,
                            "aa_mut": aa_mut, "mut_type": mut_type, "marker": "Yes",
                            "m_ids": " | ".join(m_ids), "is_combo": is_combo,
                            "effect": effect, "found": found, "ref": ref})
        else:
            mut_type = "Insertion" if r_aa == "-" else "Deletion" if q_aa == "-" else "Substitution"
            aa_mut = f"{r_aa}{pos_raw}{q_aa}"
            events.append({"pos": pos_raw, "pos_ref": pos_ref, "r_aa": r_aa, "q_aa": q_aa,
                            "aa_mut": aa_mut, "mut_type": mut_type, "marker": "No",
                            "m_ids": "", "is_combo": "No", "effect": "", "found": "", "ref": ""})

    # Collapse consecutive indels into single events to avoid multiple rows
    def collapse_indels(events):
        result = []
        i = 0
        while i < len(events):
            indel = events[i]
            if indel["mut_type"] in ("Insertion", "Deletion"):
                # Collect consecutive events of the same indel type
                group = [indel]
                j = i + 1
                while j < len(events) and events[j]["mut_type"] == indel["mut_type"]:
                    group.append(events[j])
                    j += 1
                # Merge: position is the first one; AAs are concatenated
                first = group[0]
                ref_aas   = "".join(e["r_aa"] for e in group)
                query_aas = "".join(e["q_aa"] for e in group)
                start_pos = first["pos"]
                end_pos   = group[-1]["pos"]
                pos_label = start_pos if start_pos == end_pos else f"{start_pos}-{end_pos}"
                pos_ref_label = first["pos_ref"] if first["pos_ref"] == group[-1]["pos_ref"] \
                                else f"{first['pos_ref']}-{group[-1]['pos_ref']}"
                if indel["mut_type"] == "Deletion":
                    aa_mut = f"del{pos_label}"
                    query_aas = "-"
                else:
                    aa_mut = f"ins{pos_label}{query_aas}"
                    ref_aas = "-"
                result.append({**first,
                                "pos": pos_label, "pos_ref": pos_ref_label,
                                "r_aa": ref_aas, "q_aa": query_aas, "aa_mut": aa_mut})
                i = j
            else:
                result.append(indel)
                i += 1
        return result

    events = collapse_indels(events)

    # --- Write output CSV ---
    out_csv = out_dir / f"${sample_id}_{prot_name}_mutations.csv"
    output_files.append(out_csv)

    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["SAMPLE_ID", "SUBTYPE", "PROTEIN", "REF_SUBTYPE", "POSITION", "POSITION_REF",
                         "REFERENCE_AA", "QUERY_AA", "AA_MUTATION", "MUTATION_TYPE", "MARKER",
                         "MARKER_ID", "IS_COMBINATION", "EFFECT", "FOUND_IN", "REFERENCE"])
        for indel in events:
            writer.writerow(["${sample_id}", subtype_val, prot_name, ref_tag,
                             indel["pos"], indel["pos_ref"], indel["r_aa"], indel["q_aa"],
                             indel["aa_mut"], indel["mut_type"], indel["marker"],
                             indel["m_ids"], indel["is_combo"], indel["effect"], indel["found"], indel["ref"]])

master_csv = Path(f"samples/${sample_id}/${sample_id}_mutations.csv")
if output_files:
    with open(master_csv, 'w', newline='') as mf:
        writer = csv.writer(mf)
        writer.writerow(["SAMPLE_ID", "SUBTYPE", "PROTEIN", "REF_SUBTYPE", "POSITION", "POSITION_REF", "REFERENCE_AA", "QUERY_AA", "AA_MUTATION", "MUTATION_TYPE", "MARKER", "MARKER_ID", "IS_COMBINATION", "EFFECT", "FOUND_IN", "REFERENCE"])
        for f_path in output_files:
            with open(f_path, 'r') as f:
                next(f)
                writer.writerows(csv.reader(f))
"""
}