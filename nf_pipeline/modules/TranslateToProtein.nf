#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TranslateToProtein {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(cds_files), path(aligned_cds_files)

    output:
    tuple val(sample_id), path("samples/${sample_id}/proteins/*_PROT.fasta"), emit: results
    tuple val(sample_id), path("${sample_id}_*_PROT_aligned.fasta"), emit: aligned
    tuple val(sample_id), path("TPerrors.log"), optional: true, emit: errors


    script:
    
    """
    mkdir -p "samples/${sample_id}/proteins"
    for cds_fasta in *_CDS.fasta; do
        prot_fasta="\$(basename "\${cds_fasta}" _CDS.fasta)_PROT.fasta"
        
        if [[ -f "\${cds_fasta}" ]]; then
            # Use seqkit translate to convert the CDS fasta to protein fasta
            seqkit translate "\${cds_fasta}" > "samples/${sample_id}/proteins/\${prot_fasta}"
        else
            echo "TranslateToProtein: CDS FASTA file \${cds_fasta} not found for sample ${sample_id}, skipping translation for this file." >> "TPerrors.log"
        fi
    done
    # Translate the aligned CDS files sequentially right after
    for aligned_cds in ${aligned_cds_files}; do
        if [[ -f "\${aligned_cds}" ]]; then
            protname=\$(echo "\${aligned_cds}" | cut -d'_' -f2)
            output_name="${sample_id}_\${protname}_PROT_aligned.fasta"
            
            seqkit translate "\${aligned_cds}" > "\${output_name}"
        else
            echo "TranslateToProtein: Aligned CDS file \${aligned_cds} not found for sample ${sample_id}, skipping." >> "TPerrors.log"
        fi
    done
    """
}
    

