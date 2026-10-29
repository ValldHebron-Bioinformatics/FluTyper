#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process OrganizeBySample {
    errorStrategy 'ignore'
    debug true 
    
    input:
    val(sample_id)

    output:
    tuple val(sample_id), path("samples/${sample_id}"), emit: results
    tuple val(sample_id), path("OSerrors.log"), optional: true, emit: errors
    tuple val(sample_id), path("${sample_id}_orientation.tsv"), optional: true, emit: orientation

    script:
    // Casting the parameter to a file() object forces Nextflow to stage it into the work directory
    def staged_fasta = file(params.inputFasta)
    """
    mkdir -p "samples/${sample_id}/segments"
    
    raw_sample="samples/${sample_id}/${sample_id}.fasta"
    combined_fasta="${sample_id}_combined.fasta"
    # Use seqkit to extract sequences for the given sample ID from the staged FASTA file
    # The ID must be followed by a header separator (or end the ID) so that e.g. sample 1913 does not also match 19135815
    seqkit grep -r -p "^\\Q${sample_id}\\E([|_\\[\\]]|\$)" "${staged_fasta}" > "\$raw_sample"

    if [ ! -s "\$raw_sample" ]; then
        echo "OrganizeBySample: No records found for sample ${sample_id}." >> "OSerrors.log"
        exit 0
    fi

    # Generate reverse complements and rename headers with seqkit
    seqkit seq --reverse --complement --validate-seq "\$raw_sample" | sed 's/^>/>rev_/' > rev_comp.fasta
    if [ "\${PIPESTATUS[0]}" -ne 0 ]; then
        echo "OrganizeBySample: seqkit failed to reverse-complement sample ${sample_id} (invalid sequence characters?). Orientation check used forward sequences only." >> "OSerrors.log"
    fi
    cat "\$raw_sample" rev_comp.fasta > "\$combined_fasta"

    # Orientation check using Nextclade with the appropriate minimizer index based on the protocol
    minimizer_index="${params.protocols[params.protocol].resources}/Segments_minimizers.json"
    nextclade sort -m "\${minimizer_index}" -r "${sample_id}_orientation.tsv" "\$combined_fasta"

    # Rescue and Split segments based on the orientation check results
    while read -r seq_name; do
        
        # Identify which segment each sequence belongs to
        seg_type=""
        for s in ${params.segments.join(' ')}; do
            if [[ "\$seq_name" =~ [_|]\${s}([_|]|\$) ]]; then
                seg_type="\$s"
                break
            fi
        done

        target_file="samples/${sample_id}/segments/${sample_id}_\${seg_type}.fasta"

        # Check the orientation score
        if grep -F "rev_\$seq_name" "${sample_id}_orientation.tsv" | cut -f 4 | grep -q "[0-9]"; then
            # The reverse complement is correct: extract it and remove the '_rev' tag
            seqkit grep -p "rev_\$seq_name" "\$combined_fasta" | sed 's/^>rev_/>/' >> "\$target_file"
        else
            # The original orientation is correct
            seqkit grep -p "\$seq_name" "\$combined_fasta" >> "\$target_file"
        fi
    done < <(grep "^>" "\$raw_sample" | tr -d '>')

    # Missing and Duplicate segment logging
    for seg in ${params.segments.join(' ')}; do
        target_fasta="samples/${sample_id}/segments/${sample_id}_\${seg}.fasta"
        
        if [ ! -s "\$target_fasta" ]; then
            echo "OrganizeBySample: No records found for sample ${sample_id} segment \${seg}, skipping." >> "OSerrors.log"
        else
            # Count the fasta headers to determine the number of sequences
            seq_count=\$(grep -c "^>" "\$target_fasta")
            
            if [ "\$seq_count" -gt 1 ]; then
                echo "OrganizeBySample: Multiple records (\$seq_count) found for sample ${sample_id} segment \${seg}. Check \$target_fasta for possible coinfection."
            fi
        fi
    done
    """
}