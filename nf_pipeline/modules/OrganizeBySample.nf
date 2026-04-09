#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process OrganizeBySample {
    errorStrategy 'ignore' 
    
    input:
    val(sample_id)

    output:
    tuple val(sample_id), path("samples/${sample_id}"), emit: results
    tuple val(sample_id), path("OSerrors.log"), optional: true, emit: errors
    tuple val(sample_id), path("${sample_id}_orientation.tsv"), emit: orientation

    script:
    // Casting the parameter to a file() object forces Nextflow to stage it into the work directory
    def staged_fasta = file(params.inputFasta)
    """
    mkdir -p "samples/${sample_id}/segments"
    
    raw_sample="samples/${sample_id}/${sample_id}.fasta"
    combined_fasta="${sample_id}_combined.fasta"

    seqkit grep -r -p "^${sample_id}" "${staged_fasta}" > "\$raw_sample"

    if [ ! -s "\$raw_sample" ]; then
        echo "OrganizeBySample: No records found for sample ${sample_id}." >> "OSerrors.log"
        exit 0
    fi

    # Build the combined file (Original + RevComp)
    rev_comp() {
        if [[ -n "\$header" ]]; then
            echo -e "\${header}\\n\${seq}\\n\${header}_rev" >> "\${combined_fasta}"
            echo "\$seq" | rev | tr 'ACGTRYSWKMBDHVNacgtryswkmbdhvn' 'TGCAYRSWMKVHDBNtgcayrswmkvhdbn' >> "\${combined_fasta}"
        fi
    }
    header=""; seq=""; > "\${combined_fasta}"
    while read -r line || [[ -n "\$line" ]]; do
        line=\$(echo "\$line" | tr -d '\\r')
        if [[ "\$line" == ">"* ]]; then rev_comp; header="\$line"; seq=""; else seq+="\$line"; fi
    done < "\$raw_sample"; rev_comp

    # Nextclade Orientation Check
    minimizer_index="${params.protocols[params.protocol].resources}/Segments_minimizers.json"
    nextclade sort -m "\${minimizer_index}" -r ${sample_id}_orientation.tsv "\${combined_fasta}"

    # Rescue and Split segments based on the orientation check results
    while read -r full_head; do
        clean_name=\${full_head#>} # Remove FASTA header symbol
        
        seg_type=""
        for s in ${params.segments.join(' ')}; do
            if [[ "\${full_head}" =~ [_|]\${s}[_|] ]]; then
                seg_type="\$s"
                break
            fi
        done

        target_file="samples/${sample_id}/segments/${sample_id}_\${seg_type}.fasta"

        is_rev=false
        while IFS=\$'\\t' read -r col_idx col_seqName col_dataset col_score col_hits; do
            # If the name matches the reverse sequence AND the score column is not empty
            if [[ "\$col_seqName" == "\${clean_name}_rev" && -n "\$col_score" ]]; then
                is_rev=true
                break
            fi
        done < "${sample_id}_orientation.tsv"

        if [ "\$is_rev" = true ]; then
            seqkit grep -p "\${clean_name}_rev" "\${combined_fasta}" | 
            seqkit replace -p "_rev\$" -r "" > "\$target_file"
        else
            seqkit grep -p "\${clean_name}" "\${combined_fasta}" > "\$target_file"
        fi
    done < <(grep ">" "\$raw_sample")

    # Missing segment logging
    for seg in ${params.segments.join(' ')}; do
        if [ ! -s "samples/${sample_id}/segments/${sample_id}_\${seg}.fasta" ]; then
            echo "OrganizeBySample: No records found for sample ${sample_id} segment \${seg}, skipping." >> "OSerrors.log"
        fi
    done
    """
}