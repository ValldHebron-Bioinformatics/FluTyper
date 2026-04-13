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

    # Generate reverse complements and rename headers with seqkit
    seqkit seq --reverse --complement --validate-seq "\$raw_sample" | sed 's/^>/>rev_/' > rev_comp.fasta
    cat "\$raw_sample" rev_comp.fasta > "\$combined_fasta"

    # Nextclade Orientation Check
    minimizer_index="${params.protocols[params.protocol].resources}/Segments_minimizers.json"
    nextclade sort -m "\${minimizer_index}" -r "orientation.tsv" "\$combined_fasta"

    # Rescue and Split segments based on the orientation check results
    while read -r seq_name; do
        
        # Identify which segment this sequence belongs to
        seg_type=""
        for s in ${params.segments.join(' ')}; do
            if [[ "\$seq_name" =~ [_|]\${s}[_|] ]]; then
                seg_type="\$s"
                break
            fi
        done

        target_file="samples/${sample_id}/segments/${sample_id}_\${seg_type}.fasta"

        # Check the orientation score
        if grep -F "rev_\$seq_name" "\$orientation_tsv" | cut -f 4 | grep -q "[0-9]"; then
            # The reverse complement is correct: extract it and remove the '_rev' tag
            seqkit grep -p "rev_\$seq_name" "\$combined_fasta" | sed 's/^>rev_/>/' > "\$target_file"
        else
            # The original orientation is correct
            seqkit grep -p "\$seq_name" "\$combined_fasta" > "\$target_file"
        fi
    done < <(grep "^>" "\$raw_sample" | tr -d '>')

    # Missing segment logging
    for seg in ${params.segments.join(' ')}; do
        if [ ! -s "samples/${sample_id}/segments/${sample_id}_\${seg}.fasta" ]; then
            echo "OrganizeBySample: No records found for sample ${sample_id} segment \${seg}, skipping." >> "OSerrors.log"
        fi
    done
    """
}