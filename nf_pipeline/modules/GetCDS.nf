#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GetCDS {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(sample_dir), path(inferred_subtypes)
    

    output:
    tuple val(sample_id), path("samples/${sample_dir}/segments/CDS/${sample_id}_CDS_HA.fasta")

    script:
    """
    # Definim la referència i l'entrada de l'HA
    export REFERENCES="${params.protocols[params.protocol].resources}/CDS_references.fasta"

    # Identify HA subtype
    h_subtype=\$(grep -P "^${sample_id}\t" "${inferred_subtypes}" | head -n 1 | cut -f2 | grep -oE 'H[0-9]+' | head -n 1 || true)
    pathotype=""
    if [[ "\$h_subtype" == "H5" || "\$h_subtype" == "H7" ]]; then
        pathotype=\$(grep ">" "\${REFERENCES}" | grep -E "^>\${h_subtype}_" | grep -oE "HPAI|LPAI" | head -n 1 || true)
    fi
    ref_pattern="\$h_subtype"
    if [[ -n "\$pathotype" ]]; then
        ref_pattern="\${h_subtype}.*\${pathotype}"
        echo "DEBUG: sample_dir=${sample_dir}"
        echo "DEBUG: inferred_subtypes=${inferred_subtypes}"
    fi
    mkdir -p "samples/${sample_dir}/segments/CDS"
    seqkit grep -r -p "^\$ref_pattern" "\${REFERENCES}" > HA_CDS_ref.fasta || true
    HA_INPUT_FILE="samples/${sample_dir}/segments/HA/${sample_id}_HA.fasta"
    OUTFILE="samples/${sample_dir}/segments/CDS/${sample_id}_CDS_HA.fasta"
    if [[ -s HA_CDS_ref.fasta && -s "\$HA_INPUT_FILE" ]]; then
        cat HA_CDS_ref.fasta "\$HA_INPUT_FILE" > mafft_input_HA.fasta
        mafft --auto mafft_input_HA.fasta > "\$OUTFILE"
    else
        echo "DEBUG: checking HA_INPUT_FILE: \$HA_INPUT_FILE"
        echo "DEBUG: checking OUTFILE: \$OUTFILE"
        if [ -f "\$HA_INPUT_FILE" ]; then
            echo "DEBUG: \$HA_INPUT_FILE exists"
            ls -l "\$HA_INPUT_FILE"
        else
            echo "DEBUG: \$HA_INPUT_FILE does NOT exist"
        fi
        if [[ ! -s "\$HA_INPUT_FILE" ]]; then
            echo "DEBUG: \$HA_INPUT_FILE is missing or empty" >&2
        fi
        echo ">NO_HA_CDS_FOUND" > "\$OUTFILE"
    fi

    """
}