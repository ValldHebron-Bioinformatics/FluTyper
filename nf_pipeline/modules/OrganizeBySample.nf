#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process OrganizeBySample {
    errorStrategy 'ignore' // Ignora errors i continua
    
    input:
    tuple val(sample_id), path(input_fasta)

    output:
    tuple val(sample_id), path("samples/${sample_id}")

    script:
    // Convert to absolute path so the worker can find it regardless of the work directory
    def logDir = file(params.outDir).toAbsolutePath()
    """
    mkdir -p "samples/${sample_id}"

    # seqkit -r regex and -p pattern to extract all records for the sample.
    seqkit grep -r -p ${sample_id} ${input_fasta} > "samples/${sample_id}/${sample_id}.fasta"

    # Iterate through segments defined in params
    for seg in ${params.segments.join(' ')}; do
        SEG_DIR="samples/${sample_id}/segments/\${seg}"
        SEG_FILE="\${SEG_DIR}/${sample_id}_\${seg}.fasta"
        
        mkdir -p "\${SEG_DIR}"
        
        # Search for the specific sample and segment combination
        seqkit grep -r -p "${sample_id}[_|]\${seg}[_|]" "${input_fasta}" > "\${SEG_FILE}"
        
        # Check if the generated file is empty (size 0)
        if [ ! -s "\${SEG_FILE}" ]; then
            rm "\${SEG_FILE}"
            if [ "\${seg}" = "HA" ]; then
                echo "No HA segment for sample ${sample_id}. Cannot determine subtype or clade. Stopping processing for this sample." >> "${logDir}/errors.log"
                exit 1
            else
                echo "No records found for sample ${sample_id} segment \${seg}, skipping." >> "${logDir}/errors.log"
            fi
        fi
    done
    """
    }