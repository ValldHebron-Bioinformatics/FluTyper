process MutationsFinder {
    errorStrategy 'ignore'
    debug true

    input:
    tuple val(sample_id), path(prot_files), val(h_tag), val(n_tag), val(pathotype)

    output:
    tuple val(sample_id), path("samples/${sample_id}/mutations/${sample_id}_*_mutations.csv"), optional: true, emit: results
    tuple val(sample_id), path("MFerrors.log"), optional: true, emit: errors

    script:
    """
    DICTIONARY="${params.protocols[params.protocol].resources}/AA_Sites.csv"
    MARKERS_dir="${params.protocols[params.protocol].resources}/MARKERS"
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
        echo "Searching REFERENCE_PROT for pattern: \$pattern"
        
        seqkit grep -r -p "\$pattern" "\$REFERENCE_PROT" > ref_prot.fasta
        
        if [[ ! -s ref_prot.fasta ]]; then
            echo "FAIL: No reference sequence found for pattern \$pattern."
            echo "MutationsFinder: No reference sequence found for pattern \$pattern in sample ${sample_id}." >> "MFerrors.log"
            continue
        fi

        echo "SUCCESS: Reference sequence found and extracted."

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
        
        echo "SAMPLE_ID,SUBTYPE,PROTEIN,REFERENCE,POSITION,REFERENCE_AA,QUERY_AA" > "samples/${sample_id}/mutations/${sample_id}_\${prot_name}_mutations.csv"
        
        # Assemble the data first, then pipe it directly into the loop
        paste <(echo "\$ref_seq" | fold -w1) <(echo "\$query_seq" | fold -w1) | while read -r ref_aa query_aa; do
            if [[ -n "\$ref_aa" && -n "\$query_aa" && "\$ref_aa" != "\$query_aa" ]]; then
                echo "Mutation found at position \$pos: reference \${ref_aa} -> query \${query_aa}"
                echo "${sample_id},\${subtype_val},\${prot_name},\${ref_pattern},\${pos},\${ref_aa},\${query_aa}" >> "samples/${sample_id}/mutations/${sample_id}_\${prot_name}_mutations.csv"
            fi
            ((pos++))
        done
    done
    """ 
}
   


//    if [[ "${h_tag}" == "H1" || "${h_tag}" == "H3" || "${h_tag}" == "H5" || "${h_tag}" == "H7" || "${h_tag}" == "H9" ]]; then
//        TARGET_H="${h_tag}"
//    else
//        TARGET_H="H5"
//        echo "MutationsFinder: Subtype ${h_tag} not found in dictionary for sample ${sample_id}. Defaulting to H5 numbering." >> "MFerrors.log"
//    fi
//
//    FINAL_MARKERS="HA-SP_\${TARGET_H}.csv"
//    if [[ ${h_tag} != "H5" ]]; then
//      python3 "${params.programs.MutationsDictionary}" \\
//        --subtype \${TARGET_H} \\
//        --markers "\$MARKERS_dir/HA-SP.csv" \\
//        --dictionary "\$DICTIONARY" \\
//        --output \$FINAL_MARKERS
//    fi
//    
//    mkdir -p "samples/${sample_id}/mutations"
//
//    for prot_file in ${prot_files}; do
//            prot_name=\$(basename "\$prot_file" | grep -oE "HA-SP|NA|PB1|PB1-F2|PB2|PA|NP|M1|M2|NS1|NS2|PA-X" | head -n 1)
//            if [[ -z "\$prot_name" ]]; then
//                continue
//            fi
//        output_csv="samples/${sample_id}/mutations/${sample_id}_\${prot_name}_mutations.csv"
//        
//        # Determine the correct reference marker file
//        # Use the translated HA markers if the protein is HA-SP and translation was performed
//        if [[ "\$prot_name" == "HA-SP" && -f "\$FINAL_MARKERS" ]]; then
//            ref_file="\$FINAL_MARKERS"
//        else
//            ref_file="\$MARKERS_dir/\${prot_name}.csv"
//        fi
//
//        # Extract sequence into a clean string
//        sequence=\$(grep -v ">" "\$prot_file" | tr -d '\\n')
//        
//        if [[ -f "\$ref_file" ]]; then
//            # Initialize the output file with the header from the reference
//            head -n 1 "\$ref_file" > "\$output_csv"
//            
//            # Read markers line by line (skipping header)
//            tail -n +2 "\$ref_file" | while IFS=, read -r pos aa rest; do
//                # Check if position is within sequence bounds
//                if [[ \$pos -le \${#sequence} && \$pos -gt 0 ]]; then
//                    # Extract the amino acid at the specified position (0-indexed adjustment)
//                    actual_aa=\${sequence:\$((pos-1)):1}
//                    
//                    # If match is found, record the marker data
//                    if [[ "\$actual_aa" == "\$aa" ]]; then
//                        echo "\$pos,\$aa,\$rest" >> "\$output_csv"
//                    fi
//                fi
//            done
//        else
//            echo "MutationsFinder: Reference file for \${prot_name} not found at \$ref_file." >> "MFerrors.log"
//            # Ensure the file exists to satisfy the output pattern
//            touch "\$output_csv"
//        fi
//    done
//   