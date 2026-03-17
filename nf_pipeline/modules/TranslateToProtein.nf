#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TranslateToProtein {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(cds_files)

    output:
    tuple val(sample_id), path("samples/${sample_id}/proteins/*_PROT.fasta"), emit: results
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
    """
}
    

