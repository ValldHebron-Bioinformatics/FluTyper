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
        echo "---"
        echo "Inspecting file: \$file for sample: ${sample_id}"

        if [[ ! -f "\$file" ]]; then
            echo "FAIL: Protein file \$file not found."
            echo "MutationsFinder: Protein file \$file not found for sample ${sample_id}." >> "MFerrors.log"
            continue
        fi

        prot_name=\$(basename "\$file" | cut -d'_' -f2)
        echo "Identified protein: \$prot_name"

        case "\$prot_name" in
            "NA")
                ref_tag="${n_tag}"
                ref_patho=""
                ;;
            HA*)
                ref_tag="${h_tag}"
                ref_patho="${pathotype}"
                ;;
            *)
                case "${h_tag}" in
                    H5|H7|H9)
                        ref_tag="${h_tag}"
                        ref_patho="${pathotype}"
                        ;;
                    *)
                        ref_tag="H5"
                        ref_patho="HPAI"
                        ;;
                esac
                ;;
        esac

        pattern="^\${ref_tag}_\${prot_name}_.*\${ref_patho}"        
        seqkit grep -r -p "\$pattern" "\$REFERENCE_PROT" > ref_prot.fasta
        
        if [[ ! -s ref_prot.fasta ]]; then
            echo "MutationsFinder: No reference sequence found for pattern \$pattern in sample ${sample_id}." >> "MFerrors.log"
            continue
        fi

        ref_seq=\$(grep -v ">" "ref_prot.fasta" | tr -d '\\n')
        query_seq=\$(grep -v ">" "\$file" | tr -d '\\n')
        
        pos=1
        
        if [[ "${h_tag}" =~ ^(H5|H7|H9)\$ ]] && [[ -n "${pathotype}" ]]; then
            subtype_val="${h_tag}${n_tag}(${pathotype})"
        else
            subtype_val="${h_tag}${n_tag}"
        fi

        if [[ -n "\${ref_patho}" ]]; then
            ref_pattern="\${ref_tag}(\${ref_patho})"
        else
            ref_pattern="\${ref_tag}"
        fi
        output_csv="samples/${sample_id}/mutations/${sample_id}_\${prot_name}_mutations.csv"
        echo "SAMPLE_ID,SUBTYPE,PROTEIN,REF_SUBTYPE,POSITION,REFERENCE_AA,QUERY_AA,MARKER,ORIGIN,EFFECT,REFERENCE" > "\$output_csv"
        
        # Assemble the data first, then pipe it directly into the loop
        paste <(echo "\$ref_seq" | fold -w1) <(echo "\$query_seq" | fold -w1) | while read -r ref_aa query_aa; do
            if [[ -n "\$ref_aa" && -n "\$query_aa" && "\$ref_aa" != "\$query_aa" ]]; then
                echo "Mutation found at position \$pos: reference \${ref_aa} -> query \${query_aa}"
                echo "${sample_id},\${subtype_val},\${prot_name},\${ref_pattern},\${pos},\${ref_aa},\${query_aa}" >> "\$output_csv"
            fi
            ((pos++))
        done
        case "${h_tag}" in
            H1|H3|H5|H7|H9)
                target_H="${h_tag}"
            *)
                target_H="H5"
                echo "MutationsFinder: Subtype ${h_tag} not found in dictionary for sample ${sample_id}. Defaulting to H5 numbering." >> "MFerrors.log"
        esac
        FINAL_MARKERS="HA-SP_\${target_H}.csv"
        if [[ "\${target_H}"!="H5" ]]; then
            python3 "${params.programs.MutationsDictionary}" \\
                --subtype \${target_H} \\
                --markers "\$MARKERS_DIR/HA-SP.csv" \\
                --dictionary "\$DICTIONARY" \\
                --output \$FINAL_MARKERS
        fi
        
        if [[ "\$prot_name" == "HA-SP" && -f "\$FINAL_MARKERS" ]]; then
            ref_file="\$FINAL_MARKERS"
        else
            ref_file="\$MARKERS_DIR/\${prot_name}.csv"
        fi

    done
    """ 
}
   