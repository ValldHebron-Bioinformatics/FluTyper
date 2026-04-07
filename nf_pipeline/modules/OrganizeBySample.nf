#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process OrganizeBySample {
    errorStrategy 'ignore' 
    
    input:
    val(sample_id)

    output:
    tuple val(sample_id), path("samples/${sample_id}"), emit: results
    tuple val(sample_id), path("OSerrors.log"), optional: true, emit: errors

    script:
    // Casting the parameter to a file() object forces Nextflow to stage it into the work directory
    def staged_fasta = file(params.inputFasta)
    """
    mkdir -p "samples/${sample_id}/segments"

    # seqkit -r regex and -p pattern to extract all records for the sample.
    seqkit grep -r -p ${sample_id} "${staged_fasta}" > "samples/${sample_id}/${sample_id}.fasta"

    # Iterate through segments defined in params
    for seg in ${params.segments.join(' ')}; do
        SEG_FILE="${sample_id}_\${seg}.fasta"
        
        
        # Search for the specific sample and segment combination
        seqkit grep -r -p "${sample_id}[_|]\${seg}[_|]" "${staged_fasta}" > "\${SEG_FILE}"
        
        # Check if the generated file is empty (size 0)
        if [ ! -s "\${SEG_FILE}" ]; then
            echo "OrganizeBySample: No records found for sample ${sample_id} segment \${seg}, skipping." >> "OSerrors.log"
        else
            mv "\${SEG_FILE}" "samples/${sample_id}/segments/\${SEG_FILE}"
        fi
    done
    """
}