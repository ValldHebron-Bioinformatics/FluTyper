#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TranslateToProtein {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(cds_files), path(sample_dir)

    output:
    tuple val(sample_id), path("samples/${sample_dir}/proteins/*_PROT.fasta")


    script:
    """
    mkdir -p "samples/${sample_dir}/proteins"
    for cds_fasta in *_CDS.fasta; do
        prot_fasta="\$(basename "\${cds_fasta}" _CDS.fasta)_PROT.fasta"
        
        if [[ -f "\${cds_fasta}" ]]; then
            # Use seqkit translate to convert the CDS fasta to protein fasta
            seqkit translate "\${cds_fasta}" > "samples/${sample_dir}/proteins/\${prot_fasta}"
        fi
    done
    """
}
    

