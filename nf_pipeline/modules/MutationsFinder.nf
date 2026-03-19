#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsFinder {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(prot_files), val(h_tag), val(n_tag), val(pathotype)

    output:
    tuple val(sample_id), path("samples/${sample_id}/mutations/${sample_id}_*_mutations.csv"), optional: true, emit: results
    tuple val(sample_id), path("MFerrors.log"), optional: true, emit: errors

    script:
    """
    DICTIONARY="${params.protocols[params.protocol].resources}/AA_Sites.csv"
    MARKERS_DIR="${params.protocols[params.protocol].resources}/MARKERS"
    REFERENCE_PROT="${params.protocols[params.protocol].resources}/PROT_references.fasta"

    mkdir -p "samples/${sample_id}/mutations"

    for file in ${prot_files}; do
        [[ ! -f "\$file" ]] && continue
        
        prot_name=\$(basename "\$file" | cut -d'_' -f2)
        
        # Determine reference pattern
        case "\$prot_name" in
            "NA") ref_tag="${n_tag}"; ref_patho="" ;;
            HA*) ref_tag="${h_tag}"; ref_patho="${pathotype}" ;;
            *) if [[ "${h_tag}" =~ ^(H5|H7|H9)\$ ]]; then
                ref_tag="${h_tag}"
                ref_patho="${pathotype}"
                else
                ref_tag="H5"
                ref_patho="HPAI"
                fi
        esac

        pattern="^\${ref_tag}_\${prot_name}_.*\${ref_patho}"        
        seqkit grep -r -p "\$pattern" "\$REFERENCE_PROT" > ref_prot.fasta || continue
        
        ref_seq=\$(grep -v ">" "ref_prot.fasta" | tr -d '\\n')
        query_seq=\$(grep -v ">" "\$file" | tr -d '\\n')
        
        # Set subtype and reference values
        if [[ "${h_tag}" =~ ^(H5|H7|H9)\$ ]] && [[ -n "${pathotype}" ]]; then
            subtype_val="${h_tag}${n_tag}(${pathotype})"
        else
            subtype_val="${h_tag}${n_tag}"
        fi
        ref_pattern="\${ref_tag}\${ref_patho:+(\$ref_patho)}"
        
        output_csv="samples/${sample_id}/mutations/${sample_id}_\${prot_name}_mutations.csv"
        echo "SAMPLE_ID,SUBTYPE,PROTEIN,REF_SUBTYPE,POSITION,REFERENCE_AA,QUERY_AA,MARKER,ORIGIN,EFFECT,REFERENCE" > "\$output_csv"
        
        # Get marker file
        target_H=\$(echo "${h_tag}" | grep -E '^H[1357912]\$' || echo "H5")
        FINAL_MARKERS="HA-SP_\${target_H}.csv"
        
        if [[ "\${target_H}" != "H5" ]]; then
            python3 "${params.programs.MutationsDictionary}" --subtype \${target_H} --markers "\$MARKERS_DIR/HA-SP.csv" --dictionary "\$DICTIONARY" --output \$FINAL_MARKERS
        fi
        
        if [[ "\$prot_name" == "HA-SP" && -f "\$FINAL_MARKERS" ]]; then
            ref_file="\$FINAL_MARKERS"
        else
            ref_file="\$MARKERS_DIR/\${prot_name}.csv"
        fi
        
        # Find mutations and check markers
        pos=1
        paste <(echo "\$ref_seq" | fold -w1) <(echo "\$query_seq" | fold -w1) | while read -r ref_aa query_aa; do
            if [[ -n "\$ref_aa" && -n "\$query_aa" && "\$ref_aa" != "\$query_aa" ]]; then
                # Check if mutation matches a marker
                marker_found=false
                if [[ -f "\$ref_file" ]]; then
                    while IFS=, read -r marker_pos marker_aa marker_origin marker_effect marker_ref; do
                        marker_pos=\$(echo "\$marker_pos" | tr -d '\\r')
                        marker_aa=\$(echo "\$marker_aa" | tr -d '\\r')
                        if [[ "\$marker_pos" == "\$pos" && "\$marker_aa" == "\$query_aa" ]]; then
                            echo "${sample_id},\${subtype_val},\${prot_name},\${ref_pattern},\${pos},\${ref_aa},\${query_aa},TRUE,\${marker_origin},\${marker_effect},\${marker_ref}" >> "\$output_csv"
                            marker_found=true
                            break
                        fi
                    done < <(tail -n +2 "\$ref_file")
                fi
                
                if [[ "\$marker_found" == "false" ]]; then
                    echo "${sample_id},\${subtype_val},\${prot_name},\${ref_pattern},\${pos},\${ref_aa},\${query_aa},FALSE,,," >> "\$output_csv"
                fi
            fi
            ((pos++))
        done
        
        # Check for non-mutation markers
        if [[ -f "\$ref_file" ]]; then
            tail -n +2 "\$ref_file" | while IFS=, read -r marker_pos marker_aa marker_origin marker_effect marker_ref; do
                marker_pos=\$(echo "\$marker_pos" | tr -d '\\r')
                marker_aa=\$(echo "\$marker_aa" | tr -d '\\r')
                if [[ "\$marker_pos" =~ ^[0-9]+\$ ]]; then
                    q_aa=\$(echo "\$query_seq" | cut -c "\$marker_pos")
                    r_aa=\$(echo "\$ref_seq" | cut -c "\$marker_pos")
                    if [[ "\$marker_aa" == "\$q_aa" && "\$r_aa" == "\$q_aa" ]]; then
                        echo "${sample_id},\${subtype_val},\${prot_name},\${ref_pattern},\${marker_pos},\${r_aa},\${q_aa},TRUE,\${marker_origin},\${marker_effect},\${marker_ref}" >> "\$output_csv"
                    fi
                fi
            done
        fi
    done
    """ 
}
