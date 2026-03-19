process MutationsFinder {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(prot_files), val(h_tag), val(n_tag), val(pathotype)

    output:
    tuple val(sample_id), path("samples/${sample_id}/mutations/${sample_id}_*_mutations.csv"), path("samples/${sample_id}/${sample_id}_mutations.csv"), optional: true, emit: results
    tuple val(sample_id), path("MFerrors.log"), optional: true, emit: errors

    script:
    """
    DICTIONARY="${params.protocols[params.protocol].resources}/AA_Sites.csv"
    MARKERS_DIR="${params.protocols[params.protocol].resources}/MARKERS"
    REFERENCE_PROT="${params.protocols[params.protocol].resources}/PROT_references.fasta"
    target_H=\$(echo "${h_tag}" | grep -E '^H[13579]\$' || echo "H5")
    FINAL_MARKERS="HA-SP_\${target_H}.csv"
    if [[ "\${target_H}" != "H5" ]]; then
        python3 "${params.programs.MutationsDictionary}" --subtype \${target_H} --input "\$MARKERS_DIR/HA-SP.csv" --dictionary "\$DICTIONARY" --output \$FINAL_MARKERS
    fi
    
    mkdir -p samples/${sample_id}/mutations 

    # PROCESS EACH PROTEIN FILE
    for file in ${prot_files}; do
        if [[ ! -f "\$file" ]];then
            echo "Warning: Protein file '\$file' not found. Skipping." >> "MFerrors.log"
            continue
        fi 
        
        prot_name=\$(basename "\$file" | cut -d'_' -f2)
        
        case "\$prot_name" in
            NA) ref_tag="${n_tag}"; ref_patho="" ;;
            HA*) ref_tag="${h_tag}"; ref_patho="${pathotype}" ;;
            *) 
                if [[ "${h_tag}" =~ ^(H5|H7|H9)\$ ]]; then
                    ref_tag="${h_tag}"
                    ref_patho="${pathotype}"
                else
                    ref_tag="H5"
                    ref_patho="HPAI"
                fi
                ;;
        esac

        pattern="^\${ref_tag}_\${prot_name}_.*\${ref_patho}"        
        seqkit grep -r -p "\$pattern" "\$REFERENCE_PROT" > ref_prot.fasta || continue
        
        ref_seq=\$(grep -v ">" "ref_prot.fasta" | tr -d '\\n')
        query_seq=\$(grep -v ">" "\$file" | tr -d '\\n')
        
        if [[ "${h_tag}" =~ ^(H5|H7|H9)\$ ]] && [[ -n "${pathotype}" ]]; then
            subtype_val="${h_tag}${n_tag}(${pathotype})"
        else
            subtype_val="${h_tag}${n_tag}"
        fi
        ref_pattern="\${ref_tag}\${ref_patho:+(\$ref_patho)}"
        
        output_csv="samples/${sample_id}/mutations/${sample_id}_\${prot_name}_mutations.csv"
        # HEADER
        echo "SAMPLE_ID,SUBTYPE,PROTEIN,REF_SUBTYPE,POSITION,REFERENCE_AA,QUERY_AA,MARKER,ORIGIN,EFFECT,REFERENCE" > "\$output_csv"
        
        if [[ "\$prot_name" == "HA-SP" && -f "\$FINAL_MARKERS" ]]; then
            ref_file="\$FINAL_MARKERS"
        else
            ref_file="\$MARKERS_DIR/\${prot_name}.csv"
        fi
        
        # FIND MUTATIONS
        pos=1
        paste <(echo "\$ref_seq" | fold -w1) <(echo "\$query_seq" | fold -w1) | while read -r ref_aa query_aa; do
            if [[ -n "\$ref_aa" && -n "\$query_aa" && "\$ref_aa" != "\$query_aa" ]]; then
                marker_found=false
                if [[ -f "\$ref_file" ]]; then
                    while IFS=, read -r m_pos m_aa m_origin m_effect m_ref; do
                        if [[ "\$m_pos" == "\$pos" && "\$m_aa" == "\$query_aa" ]]; then
                            echo "${sample_id},\${subtype_val},\${prot_name},\${ref_pattern},\${pos},\${ref_aa},\${query_aa},TRUE,\${m_origin},\${m_effect},\${m_ref}" >> "\$output_csv"
                            marker_found=true
                            break
                        fi
                    done < <(tail -n +2 "\$ref_file" | tr -d '\\r') # Skip header and sanitize line endings
                fi
                
                if [[ "\$marker_found" == "false" ]]; then
                    echo "${sample_id},\${subtype_val},\${prot_name},\${ref_pattern},\${pos},\${ref_aa},\${query_aa},FALSE,,," >> "\$output_csv"
                fi
            fi
            ((pos++))
        done
        
        # CHECK NON-MUTATION MARKERS
        if [[ -f "\$ref_file" ]]; then
            while IFS=, read -r m_pos m_aa m_origin m_effect m_ref; do
                if [[ "\$m_pos" =~ ^[0-9]+\$ ]]; then
                    q_aa="\${query_seq:m_pos-1:1}"
                    r_aa="\${ref_seq:m_pos-1:1}"

                    if [[ "\$m_aa" == "\$q_aa" && "\$r_aa" == "\$q_aa" ]]; then
                        echo "${sample_id},\${subtype_val},\${prot_name},\${ref_pattern},\$m_pos,\$r_aa,\$q_aa,TRUE,\$m_origin,\$m_effect,\$m_ref" >> "\$output_csv"
                    fi
                fi
            done < <(tail -n +2 "\$ref_file" | tr -d '\\r')
        fi
        # CONVERT HA TO H5 NUMBERING
        if [[ "\$prot_name" == HA* && "\${target_H}" != "H5" ]]; then            
            python3 "${params.programs.MutationsDictionary}" \
                --base "\${target_H}" \
                --subtype "H5" \
                --input "\${output_csv}" \
                --dictionary "\$DICTIONARY" \
                --output "\$output_csv"
        fi
    done
    # COMPILE ALL PROTEIN MUTATION FILES INTO ONE MASTER CSV
    # Grab the header
    head -n 1 \$(ls samples/${sample_id}/mutations/${sample_id}_*_mutations.csv | head -n 1) > samples/${sample_id}/${sample_id}_mutations.csv
    # Append all lines EXCEPT the header -q to avoid printing file names, -n +2 to skip the header
    tail -q -n +2 samples/${sample_id}/mutations/${sample_id}_*_mutations.csv >> samples/${sample_id}/${sample_id}_mutations.csv
    """ 
}