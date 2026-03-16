#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MutationsFinder {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(prot_files), val(h_tag)

    output:
    tuple val(sample_id), path("samples/${sample_id}/NUMERATION/${sample_id}_*_numeration.csv")

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

    FINAL_MARKERS="translated_markers_\${TARGET_H}.xlsx"

    python3 "${params.programs.MutationsDictionary}" \\
        --subtype \${TARGET_H} \\
        --markers "\$MARKERS_dir/HA.csv" \\
        --dictionary "\$DICTIONARY" \\
        --output \$FINAL_MARKERS

    mkdir -p "samples/${sample_id}/NUMERATION"

    for prot_file in ${prot_files}; do
        prot_name=\$(basename "\$prot_file" .fasta)
        temp_csv="samples/${sample_id}/NUMERATION/\${prot_name}_numeration.csv"
            sequence=\$(cat "\$prot_file" | grep -v ">" | tr -d '\n')
            echo "position,aa" > "\$temp_csv"
            num=1
            for ((i=0; i<\${#sequence}; i++)); do
                aa=\${sequence:\$i:1}
                echo "\$num,\$aa" >> "\$temp_csv"
                num=\$((\$num + 1))
            done
        
    done
    """
}