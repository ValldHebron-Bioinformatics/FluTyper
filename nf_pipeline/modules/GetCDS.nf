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
    import os, subprocess, io
    from Bio import SeqIO

    ref_fasta = "${params.protocols[params.protocol].resources}/CDS_references.fasta"
    cds_dir = "samples/${sample_id}/CDS"
    os.makedirs(cds_dir, exist_ok=True)
    log_file = "CDSerrors.log"

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

                # Homology check: Calculate coverage and identity to ensure we have a real match before proceeding with trimming and CDS extraction.
                matches = 0
                aligned_positions = 0

                # Count position by position (ignoring spaces where either sequence has a gap)
                for r, q in zip(ref_seq, aligned_seq):
                    if r != '-' and q != '-':
                        aligned_positions += 1
                        if r.upper() == q.upper():
                            matches += 1
                
                ref_length = len(ref_seq.replace('-', ''))
                
                coverage_ratio = aligned_positions / ref_length if ref_length > 0 else 0
                identity_ratio = matches / aligned_positions if aligned_positions > 0 else 0

                # Minimum 30% coverage AND minimum 60% identity
                if coverage_ratio < 0.3 or identity_ratio < 0.6:
                    with open(log_file, 'a') as f:
                        f.write(f"GetCDS: ${sample_id} {prot} ignored. No real homology (Coverage: {coverage_ratio:.1%}, Identity: {identity_ratio:.1%}).\\n")
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