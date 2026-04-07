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
import os, csv, re, subprocess
from pathlib import Path
from Bio import SeqIO

markers_dir = Path("${params.protocols[params.protocol].resources}/markers")
dictionary = "${params.protocols[params.protocol].resources}/AA_DICT.csv"
mutations_prog = "${params.programs.MutationsDictionary}"
log_file = "MFerrors.log"
out_dir = Path("samples/${sample_id}/mutations")
out_dir.mkdir(parents=True, exist_ok=True)
output_files = []
aligned_prots = "${prot_files}".split()

for aligned_prot in aligned_prots:
    file_path = Path(aligned_prot)
    prot_name = file_path.name.replace("_PROT_aligned.fasta", "").split('_')[1]

    records = list(SeqIO.parse(file_path, "fasta"))
    if len(records) < 2:
        with open(log_file, 'a') as log_f:
            log_f.write(f"ERROR: Less than 2 sequences found in ${sample_id} {prot_name} alignment file.\\n")
        continue

    ref_header = str(records[0].description)
    ref_seq = str(records[0].seq)
    query_header = str(records[1].id)
    query_seq = str(records[1].seq)
    
    ref_tag = ref_header.split('_')[0]
    subtype_val = "${h_tag}${n_tag}(${pathotype})" if "${pathotype}" != "" else "${h_tag}${n_tag}"
    target_H = "${h_tag}" if "${h_tag}" in {"H1", "H3", "H5", "H7", "H9"} else "H5"
    markers = {}
    
    m_file = markers_dir / f"{prot_name}_markers.csv"
    if m_file.exists():
        with open(m_file) as f:
            for row in csv.DictReader(f):
                markers[(row['POSITION'].strip(), row['AA'].strip())] = (row['EFFECT'], row.get('REFERENCE',''))

    output_csv = out_dir / f"${sample_id}_{prot_name}_mutations.csv"
    output_files.append(output_csv)
    
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        # HEADER
        writer.writerow(["SAMPLE_ID", "SUBTYPE", "PROTEIN", "REF_SUBTYPE", "POSITION", "POSITION_SUBTYPE", "REFERENCE_AA", "QUERY_AA", "AA_MUTATION", "MUTATION_TYPE", "MARKER", "EFFECT", "REFERENCE"])
        
        pos = 0
        for ref_aa, query_aa in zip(ref_seq, query_seq):
            # Only increment position if reference amino acid is not a gap
            if ref_aa != "-":
                pos += 1

            pos_str = str(pos)
            is_marker = (pos_str, query_aa) in markers
            is_mutation = ref_aa != query_aa

            pos_subtype = pos_str if prot_name.startswith("HA") else ""

            if is_marker:
                aa_mutation = f"{pos_str}{query_aa}"
                m_info = markers[(pos_str, query_aa)]
                if not is_mutation:
                    mutation_type = "None"
                elif ref_aa == "-" and query_aa != "-":
                    mutation_type = "Insertion"
                elif ref_aa != "-" and query_aa == "-":
                    mutation_type = "Deletion"
                else:
                    mutation_type = "Substitution"
                    
                writer.writerow(["${sample_id}", subtype_val, prot_name, ref_tag, pos_str, pos_subtype, ref_aa, query_aa, aa_mutation, mutation_type, "Yes" , m_info[0], m_info[1]])

            elif is_mutation:
                aa_mutation = f"{ref_aa}{pos_str}{query_aa}"
                if ref_aa == "-" and query_aa != "-":
                    mutation_type = "Insertion"
                elif ref_aa != "-" and query_aa == "-":
                    mutation_type = "Deletion"
                else:
                    mutation_type = "Substitution"
                writer.writerow(["${sample_id}", subtype_val, prot_name, ref_tag, pos_str, pos_subtype, ref_aa, query_aa, aa_mutation, mutation_type, "No", "", ""])
            else:
                continue
                
    # Run external dictionary script safely via subprocess
    if prot_name.startswith("HA") and target_H != "H5":
        subprocess.run(["python3", mutations_prog, "--subtype", target_H, "--input", str(output_csv), "--dictionary", dictionary, "--output", str(output_csv)], check=True)

# Compile master CSV for the sample with all mutations from individual protein files
master_csv = Path(f"samples/${sample_id}/${sample_id}_mutations.csv")
if output_files:
    with open(master_csv, 'w', newline='') as master_f:
        writer = csv.writer(master_f)
        writer.writerow(["SAMPLE_ID", "SUBTYPE", "PROTEIN", "REF_SUBTYPE", "POSITION", "POSITION_SUBTYPE", "REFERENCE_AA", "QUERY_AA", "AA_MUTATION", "MUTATION_TYPE", "MARKER", "EFFECT", "REFERENCE"]) # Write header once

        for f_path in output_files:
            with open(f_path, 'r') as f:
                reader = csv.reader(f)
                next(reader) # Skip header
                writer.writerows(reader) 
"""
}