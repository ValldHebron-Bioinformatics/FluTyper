#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GetCDS {
    errorStrategy 'ignore'

    input:
    tuple val(h_tag), val(n_tag), val(sample_id), val(pathotype), path(sample_dir)

    output:
    tuple val(sample_id), path("samples/${sample_id}/segments/CDS/*_CDS.fasta")

    script:
    """
    #!/usr/bin/env python3
    import os, subprocess, re

    h_tag = "${h_tag}"
    n_tag = "${n_tag}"
    sample_id = "${sample_id}"
    pathotype = "${pathotype}"
    sample_dir = "${sample_dir}"
    references_fasta = "${params.protocols[params.protocol].resources}/CDS_references.fasta"
    
    cds_dir = f"samples/{sample_id}/segments/CDS"
    os.makedirs(cds_dir, exist_ok=True)


    # Python dictionary
    prot_dict = {
        "HA":  ["HA", "HA-SP", "HA1-SP", "HA2"],
        "NA":  ["NA"],
        "PB2": ["PB2"],
        "PB1": ["PB1", "PB1-F2"],
        "PA":  ["PA", "PA-X"],
        "NP":  ["NP"],
        "MP":  ["M1", "M2"],
        "NS":  ["NS1", "NEP"]
    }

    for seg, prots in prot_dict.items():
        if seg == "NA":
            ref_tag = n_tag
            ref_patho = ""
        elif seg == "HA":
            ref_tag = h_tag
            ref_patho = pathotype
        else:
            # Internal segments
            if h_tag not in ["H5", "H7", "H9"]:
                ref_tag = "H5"
                ref_patho = "HPAI"
            else:
                ref_tag = h_tag
                ref_patho = pathotype

        seg_fasta = f"{sample_dir}/segments/{seg}/{sample_id}_{seg}.fasta"
        if not os.path.isfile(seg_fasta):
            continue
            
        for prot in prots:
            pattern = f"^{ref_tag}_{prot}_.*_{ref_patho}" if (ref_patho and seg != "NA") else f"^{ref_tag}_{prot}_"
            ref_out = f"{cds_dir}/{sample_id}_{prot}_ref.fasta"
            mafft_in = f"{cds_dir}/{sample_id}_{prot}_mafft_in.fasta"
            mafft_out = f"{cds_dir}/{sample_id}_{prot}_aligned.fasta"
            
            subprocess.run(f"seqkit grep -r -p '{pattern}' {references_fasta} > {ref_out}", shell=True)
            
            if os.path.isfile(ref_out) and os.path.getsize(ref_out) > 0:
                os.system(f'cat "{ref_out}" "{seg_fasta}" > "{mafft_in}"')
                os.system(f'mafft --auto --quiet "{mafft_in}" > "{mafft_out}"')
                
                if os.path.isfile(mafft_out):
                    with open(mafft_out, 'r') as f:
                        parts = [p.strip().split('\\n', 1) for p in f.read().split('>') if p.strip()]
                    
                    if len(parts) >= 2:
                        ref_seq = parts[0][1].replace('\\n', '')
                        aligned_header, aligned_seq = parts[1][0], parts[1][1].replace('\\n', '')
                        
                        start = len(ref_seq) - len(ref_seq.lstrip('-'))
                        end = len(ref_seq.rstrip('-'))
                        
                        final_seq, temp_chunk = [], []
                        gap_len = 0
                        
                        for r, s in zip(ref_seq[start:end], aligned_seq[start:end]):
                            if r == '-':
                                gap_len += 1
                                if s != '-': temp_chunk.append(s)
                            else:
                                if gap_len < 50: final_seq.extend(temp_chunk)
                                gap_len, temp_chunk = 0, []
                                if s != '-': final_seq.append(s)
                        
                        if gap_len < 50: final_seq.extend(temp_chunk)
                        
                        clean_seq = "".join(final_seq)
                        with open(f"{cds_dir}/{sample_id}_{prot}_CDS.fasta", 'w') as f_out:
                            f_out.write(f">{aligned_header}\\n")
                            f_out.write('\\n'.join(clean_seq[i:i+80] for i in range(0, len(clean_seq), 80)) + '\\n') # https://softwareengineering.stackexchange.com/questions/148677/why-is-80-characters-the-standard-limit-for-code-width 
    
    """
}