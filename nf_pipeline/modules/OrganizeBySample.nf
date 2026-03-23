#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process OrganizeBySample {
    // errorStrategy 'ignore' // Ignora errors i continua
    
    input:
    val(sample_id)
    output:
    tuple val(sample_id), path("samples/${sample_id}"), emit: results
    tuple val(sample_id), path("OSerrors.log"), optional: true, emit: errors

    script:
    // Another option is to put it as output... ASK ALEJANDRA
    
    """
    mkdir -p "samples/${sample_id}"
    input_fasta="${params.inputFasta}"

    # seqkit -r regex and -p pattern to extract all records for the sample.
    seqkit grep -r -p ${sample_id} "\${input_fasta}" > "samples/${sample_id}/${sample_id}.fasta"

    # Iterate through segments defined in params
    for seg in ${params.segments.join(' ')}; do
        SEG_DIR="samples/${sample_id}/segments/\${seg}"
        SEG_FILE="\${SEG_DIR}/${sample_id}_\${seg}.fasta"
        
        mkdir -p "\${SEG_DIR}"
        
        # Search for the specific sample and segment combination
        seqkit grep -r -p "${sample_id}[_|]\${seg}[_|]" "\${input_fasta}" > "\${SEG_FILE}"
        
        # Check if the generated file is empty (size 0)
        if [ ! -s "\${SEG_FILE}" ]; then
            echo "OrganizeBySample: No records found for sample ${sample_id} segment \${seg}, skipping." >> "OSerrors.log"
        fi
    done
    """
    }