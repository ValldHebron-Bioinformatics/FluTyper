process MutationsFinder {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(prot_files), val(h_tag)

    output:
    tuple val(sample_id), path("samples/${sample_id}/mutations/${sample_id}_*_mutations.csv")

    script:
    def logDir = file(params.outDir)
    
    """
    DICTIONARY="${params.protocols[params.protocol].resources}/AA_Sites.csv"
    MARKERS_dir="${params.protocols[params.protocol].resources}/MARKERS"

    if [[ "${h_tag}" == "H1" || "${h_tag}" == "H3" || "${h_tag}" == "H5" || "${h_tag}" == "H7" || "${h_tag}" == "H9" ]]; then
        TARGET_H="${h_tag}"
    else
        TARGET_H="H5"
        echo "Subtype ${h_tag} not found in dictionary for sample ${sample_id}. Defaulting to H5 numbering." >> "${logDir}/errors.log"
    fi

    FINAL_MARKERS="HA-SP_\${TARGET_H}.csv"
    if [[ ${h_tag} != "H5" ]]; then
      python3 "${params.programs.MutationsDictionary}" \\
        --subtype \${TARGET_H} \\
        --markers "\$MARKERS_dir/HA-SP.csv" \\
        --dictionary "\$DICTIONARY" \\
        --output \$FINAL_MARKERS
    fi
    
    mkdir -p "samples/${sample_id}/mutations"

    for prot_file in ${prot_files}; do
            prot_name=\$(basename "\$prot_file" | grep -oE "HA-SP|NA|PB1|PB1-F2|PB2|PA|NP|M1|M2|NS1|NS2|PA-X" | head -n 1)
            if [[ -z "\$prot_name" ]]; then
                continue
            fi
        output_csv="samples/${sample_id}/mutations/${sample_id}_\${prot_name}_mutations.csv"
        
        # Determine the correct reference marker file
        # Use the translated HA markers if the protein is HA-SP and translation was performed
        if [[ "\$prot_name" == "HA-SP" && -f "\$FINAL_MARKERS" ]]; then
            ref_file="\$FINAL_MARKERS"
        else
            ref_file="\$MARKERS_dir/\${prot_name}.csv"
        fi

        # Extract sequence into a clean string
        sequence=\$(grep -v ">" "\$prot_file" | tr -d '\\n')
        
        if [[ -f "\$ref_file" ]]; then
            # Initialize the output file with the header from the reference
            head -n 1 "\$ref_file" > "\$output_csv"
            
            # Read markers line by line (skipping header)
            tail -n +2 "\$ref_file" | while IFS=, read -r pos aa rest; do
                # Check if position is within sequence bounds
                if [[ \$pos -le \${#sequence} && \$pos -gt 0 ]]; then
                    # Extract the amino acid at the specified position (0-indexed adjustment)
                    actual_aa=\${sequence:\$((pos-1)):1}
                    
                    # If match is found, record the marker data
                    if [[ "\$actual_aa" == "\$aa" ]]; then
                        echo "\$pos,\$aa,\$rest" >> "\$output_csv"
                    fi
                fi
            done
        else
            echo "Reference file for \${prot_name} not found at \$ref_file." >> "${logDir}/errors.log"
            # Ensure the file exists to satisfy the output pattern
            touch "\$output_csv"
        fi
    done
    """
}