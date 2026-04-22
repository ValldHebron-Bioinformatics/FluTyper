#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GetCDS {
    errorStrategy 'ignore'

    input:
    tuple val(h_tag), val(n_tag), val(sample_id), val(pathotype), path(sample_dir)

    output:
    tuple val(sample_id), path("samples/${sample_id}/CDS/*_CDS.fasta"), emit: results
    tuple val(sample_id), path("${sample_id}_*_CDS_aligned.fasta"), emit: aligned
    tuple val(sample_id), path("CDSerrors.log"), optional: true, emit: errors

    script:
    """
#!/usr/bin/env python3
import os, subprocess, io, csv
from Bio import SeqIO

ref_fasta = "${params.protocols[params.protocol].resources}/CDS_references.fasta"
threshold_csv = "${params.protocols[params.protocol].resources}/Identity_thresholds.csv"
cds_dir = "samples/${sample_id}/CDS"
os.makedirs(cds_dir, exist_ok=True)
log_file = "CDSerrors.log"

# Load identity thresholds from CSV into a dictionary
identity_thresholds = {}
if os.path.isfile(threshold_csv):
    with open(threshold_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            identity_thresholds[row['Target_Group']] = float(row['Calculated_Threshold']) / 100.0

prot_dict = {
    "HA":  ["HA1", "HA2"], "NA":  ["NA"], "PB2": ["PB2"],
    "PB1": ["PB1", "PB1-F2"], "PA":  ["PA", "PA-X"], "NP":  ["NP"],
    "MP":  ["M1", "M2"], "NS":  ["NS1", "NS2"]
}

def TrimCDS(ref_seq, aligned_seq, gap_threshold):
    start = len(ref_seq) - len(ref_seq.lstrip('-'))
    end   = len(ref_seq.rstrip('-'))
    
    final_query, final_ref = [], []
    tmp_query, tmp_ref = [], []
    gap_len = 0
    
    for ref_ch, query_ch in zip(ref_seq[start:end], aligned_seq[start:end]):
        if ref_ch == '-':
            gap_len += 1
            tmp_query.append(query_ch)
            tmp_ref.append(ref_ch)
        else:
            if gap_len < gap_threshold:
                final_query.extend(tmp_query)
                final_ref.extend(tmp_ref)
            gap_len, tmp_query, tmp_ref = 0, [], []
            
            final_query.append(query_ch)
            final_ref.append(ref_ch)
            
    if gap_len < gap_threshold:
        final_query.extend(tmp_query)
        final_ref.extend(tmp_ref)
        
    return "".join(final_ref), "".join(final_query)

PATHO_SUBTYPES = {"H5", "H7", "H9"}

for seg, prots in prot_dict.items():
    if seg == "NA":
        ref_tag, ref_patho = "${n_tag}", "${pathotype}"
    elif seg == "HA":
        ref_tag, ref_patho = "${h_tag}", "${pathotype}"
    else:
        ref_tag = "${h_tag}" if "${h_tag}" in PATHO_SUBTYPES else "H5"
        ref_patho = "${pathotype}" if "${h_tag}" in PATHO_SUBTYPES else "HPAI"

    seg_fasta = f"${sample_dir}/segments/${sample_id}_{seg}.fasta"
    
    if not (os.path.isfile(seg_fasta) and os.path.getsize(seg_fasta) > 0):
        with open(log_file, 'a') as f:
            f.write(f"GetCDS: No valid FASTA found for sample ${sample_id} segment {seg}. Skipping all proteins for this segment.\\n")
        continue 
        
    for prot in prots:
        pattern = f"^{ref_tag}_{prot}_.*{ref_patho}"
        
        # Determine the target group for identity threshold lookup
        if seg == "HA":
            target_group = f"{prot}_${h_tag}"
        elif seg == "NA":
            target_group = f"{prot}_${n_tag}"
        else:
            target_group = prot
            
        # Obtain the specific threshold (or 60% by default if something fails)
        min_identity = identity_thresholds.get(target_group, 0.60)
        
        cmd = f"seqkit grep -r -p '{pattern}' {ref_fasta} | cat - '{seg_fasta}' | mafft --localpair --maxiterate 1000 --op 3 --ep 0.123 --quiet -"
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            
            # Check for empty output from MAFFT
            if not result.stdout or not result.stdout.strip():
                with open(log_file, 'a') as f:
                    f.write(f"GetCDS: Alignment for ${sample_id} produced no data for {prot}. Protein likely missing in segment.\\n")
                continue

            sequences = list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))
            
            if len(sequences) < 2:
                with open(log_file, 'a') as f:
                    f.write(f"GetCDS: No homologous sequence found in ${sample_id} for protein {prot}. Skipping.\\n")
                continue

            ref_header = sequences[0].id
            ref_seq = str(sequences[0].seq)
            aligned_header = sequences[1].id
            aligned_seq = str(sequences[1].seq)

            # Homology check: Calculate coverage, identity, and N ratio
            matches = 0
            aligned_informative_positions = 0
            n_count = 0
            query_bases = 0

            # Count position by position
            for r, q in zip(ref_seq, aligned_seq):
                q_upper = q.upper()
                
                # Track query bases and Ns (ignoring gaps)
                if q_upper != '-':
                    query_bases += 1
                    if q_upper == 'N':
                        n_count += 1

                # Track alignment matches and coverage (ignoring gaps AND ignoring Ns)
                if r != '-' and q != '-':
                    if q_upper != 'N':
                        aligned_informative_positions += 1
                        if r.upper() == q_upper:
                            matches += 1
            
            # Check for bad reads based on N ratio
            n_ratio = n_count / query_bases if query_bases > 0 else 0
            
            if n_ratio > 0.5:
                with open(log_file, 'a') as f:
                    f.write(f"GetCDS: ${sample_id} {prot} ignored. Bad read with >50% Ns (N ratio: {n_ratio:.1%}).\\n")
                continue

            # Calculate coverage and identity based on informative sites
            ref_length = len(ref_seq.replace('-', ''))
            
            coverage_ratio = aligned_informative_positions / ref_length if ref_length > 0 else 0
            identity_ratio = matches / ref_length if ref_length > 0 else 0

            # Minimum 50% coverage AND dynamic minimum identity
            if coverage_ratio < 0.5 or identity_ratio < min_identity:
                with open(log_file, 'a') as f:
                    f.write(f"GetCDS: ${sample_id} {prot} ignored. No real homology (Coverage: {coverage_ratio:.1%}, Identity: {identity_ratio:.1%}, Required Identity: {min_identity:.1%}).\\n")
                continue
                
            if prot == "NS2":
                current_threshold = 400
            elif prot == "M2":
                current_threshold = 600
            else:
                current_threshold = float('inf')
                
            clean_ref, clean_query = TrimCDS(ref_seq, aligned_seq, current_threshold)
            
            # Final check if trimming resulted in an empty sequence
            if not clean_query.replace('-', '').strip():
                with open(log_file, 'a') as f:
                    f.write(f"GetCDS: Trimmed sequence for ${sample_id} {prot} is empty. Skipping.\\n")
                continue

            with open(f"{cds_dir}/${sample_id}_{prot}_CDS.fasta", 'w') as f_out:
                f_out.write(f">{aligned_header}\\n")
                for i in range(0, len(clean_query), 80):
                    f_out.write(clean_query[i:i+80] + '\\n')
                    
            aligned_file = f"${sample_id}_{prot}_CDS_aligned.fasta"
            with open(aligned_file, 'w') as f_aln:
                f_aln.write(f">{ref_header}_trimmed\\n{clean_ref}\\n")
                f_aln.write(f">{aligned_header}_trimmed\\n{clean_query}\\n")
                            
        except subprocess.CalledProcessError as e:
            with open(log_file, 'a') as f:
                f.write(f"GetCDS: MAFFT alignment failed for ${sample_id} {prot}: {e}\\n")
    """
}