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

    references_fasta = "${params.protocols[params.protocol].resources}/CDS_references.fasta"
    
    # Create output directory for CDS if it doesn't exist
    cds_dir = "samples/${sample_id}/segments/CDS"
    os.makedirs(cds_dir, exist_ok=True)

    # Define the map here, I wanted to keep it in params but it was too complex for Groovy interpolation
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

    # Iterate over segments and their associated proteins
    for seg, prots in prot_dict.items():
        if seg == "NA":
            ref_tag = "${n_tag}"
            ref_patho = ""
        elif seg == "HA":
            ref_tag = "${h_tag}"
            ref_patho = "${pathotype}"
        else:
            # Internal segments
            if "${h_tag}" not in ["H5", "H7", "H9"]:
                ref_tag = "H5"
                ref_patho = "HPAI"
            else:
                ref_tag = "${h_tag}"
                ref_patho = "${pathotype}"

        # Groovy interpolates \${sample_dir} and \${sample_id}, Python interpolates {seg}
        seg_fasta = f"${sample_dir}/segments/{seg}/${sample_id}_{seg}.fasta"
        if not os.path.isfile(seg_fasta):
            continue
        # For each protein associated with the segment, extract the reference, align, and trim to get the CDS    
        for prot in prots:
            pattern = f"^{ref_tag}_{prot}_.*_{ref_patho}" if (ref_patho and seg != "NA") else f"^{ref_tag}_{prot}_"
            ref_out = f"{cds_dir}/${sample_id}_{prot}_ref.fasta"
            mafft_in = f"{cds_dir}/${sample_id}_{prot}_mafft_in.fasta"
            mafft_out = f"{cds_dir}/${sample_id}_{prot}_aligned.fasta"
            # Use seqkit to extract the reference sequence based on the pattern
            subprocess.run(f"seqkit grep -r -p '{pattern}' {references_fasta} > {ref_out}", shell=True)
            
            # If we got a reference sequence, proceed with alignment
            if os.path.isfile(ref_out) and os.path.getsize(ref_out) > 0:
                os.system(f'cat "{ref_out}" "{seg_fasta}" > "{mafft_in}"')
                os.system(f'mafft --auto --quiet "{mafft_in}" > "{mafft_out}"')
                
                # If the alignment was successful, extract the CDS sequence by removing gaps
                # from the reference sequence and applying the same trimming to the aligned sequence
                if os.path.isfile(mafft_out):
                    
                    with open(mafft_out, 'r') as f: 
                        sequences = []
                        for parts in f.read().split('>'): # Split the FASTA file into entries
                            parts = parts.strip() 
                            if parts:
                                sequences.append(parts.split('\\n', 1)) # Split header and sequence, maybe use SeqIO here?

                    if len(sequences) >= 2: # We expect exactly 2 sequences, but in case there are more (e.g. multiple references matching the pattern)
                        ref_seq = sequences[0][1].replace('\\n', '')
                        aligned_header, aligned_seq = sequences[1][0], sequences[1][1].replace('\\n', '')
                        
                        start = len(ref_seq) - len(ref_seq.lstrip('-')) # Count leading gaps to find the true start of the sequence
                        end = len(ref_seq.rstrip('-')) # Count trailing gaps to find the true end of the sequence
                        
                        final_seq, tmp_chunk = [], []
                        gap_len = 0
                        
                        for ref_ch, query_ch in zip(ref_seq[start:end], aligned_seq[start:end]): # Here we already trim the reference to the expected CDS region, so we only consider that part for gap removal
                            if ref_ch == '-':
                                gap_len += 1
                                tmp_chunk.append(query_ch)
                            else:
                                if gap_len < 50: # This threshold can be adjusted based on expected indel sizes, it allows small gaps to be included in the final sequence but filters out large gaps like the alternative splicing in NEP or M2
                                    final_seq.extend(tmp_chunk)
                                gap_len, tmp_chunk = 0, []
                                final_seq.append(query_ch)

                        if gap_len < 50: 
                            final_seq.extend(tmp_chunk)

                        clean_seq = "".join(final_seq)
                        with open(f"{cds_dir}/${sample_id}_{prot}_CDS.fasta", 'w') as f_out:
                            f_out.write(f">{aligned_header}\\n")
                            f_out.write('\\n'.join(clean_seq[i:i+80] for i in range(0, len(clean_seq), 80)) + '\\n')
    """
}
