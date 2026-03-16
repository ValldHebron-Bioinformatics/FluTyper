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
        DICTIONARY="${params.protocols[params.protocol].resources}/AA_dictionary_proposal/AA_Sites.csv"
        MARKERS="${params.protocols[params.protocol].resources}/AA_dictionary_proposal/MARKERS.xlsx"

        if [[ "${h_tag}" =~ ^(H1|H3|H5|H7|H9)\$ ]]; then
            TARGET_H="${h_tag}"
        else
            TARGET_H="H5"
            echo "Subtype ${h_tag} not found in dictionary for sample ${sample_id}. Defaulting to H5 numbering." >> "${logDir}/errors.log"
        fi

        FINAL_MARKERS="translated_markers_\${TARGET_H}.xlsx"

        python3 "${params.programs.MutationsDictionary}" \\
            --subtype \${TARGET_H} \\
            --markers "\$MARKERS" \\
            --dictionary "\$DICTIONARY" \\
            --output \$FINAL_MARKERS

        mkdir -p "samples/${sample_id}/NUMERATION"

        for file in *.fasta; do
            if [ ! -e "\$file" ]; then
                echo "No protein FASTA files found for sample ${sample_id}, skipping mutation finding." >> "${logDir}/errors.log"
                break
            fi

            protein=\$(basename "\$file" | sed 's/${sample_id}_//' | sed 's/_PROT\\.fasta//')
            output_file="samples/${sample_id}/NUMERATION/${sample_id}_\${protein}_numeration.csv"

            sequence="\$(grep -v '^>' "\$file" | tr -d '\\n')"

            if [[ -n "\$sequence" ]]; then
                echo -e "AminoAcid\\tPosition" > "\$output_file"

                for ((i=0; i<\${#sequence}; i++)); do
                    aa="\${sequence:\$i:1}"
                    echo -e "\${aa}\\t\$((i+1))" >> "\$output_file"
                done
            fi
            
        done
    """
}