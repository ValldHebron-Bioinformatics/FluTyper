#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GetCDS {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(sample_dir), path(inferred_subtypes)
    

    output:
    tuple val(sample_id), path("samples/${sample_dir}/segments/CDS/${sample_id}_CDS_HA.fasta"), emit: cds_ha

    script:
    """
    # Definim la referència i l'entrada de l'HA
    export REFERENCES="${params.protocols[params.protocol].resources}/CDS_references.fasta"
    # Define output directory and file
    OUTDIR="samples/${sample_dir}/segments/CDS"
    OUTFILE="\$OUTDIR/${sample_id}_CDS_HA.fasta"
    mkdir -p "\$OUTDIR"

    # Identify HA subtype
    h_subtype=\$(grep -m 1 -E "^${sample_id}\\b" "${inferred_subtypes}" | cut -f2 | grep -oE 'H[0-9]+' | head -n 1 || true)
    pathotype=""
    if [[ "\$h_subtype" == "H5" || "\$h_subtype" == "H7" ]]; then
        pathotype=\$(grep ">" "\${REFERENCES}" | grep -E "_\${h_subtype}_" | grep -oE "HPAI|LPAI" | head -n 1 || true)
    fi
    ref_pattern="\$h_subtype"
    if [[ -n "\$pathotype" ]]; then
        ref_pattern="\${h_subtype}.*\${pathotype}"
    fi
    seqkit grep -r -p "^\$ref_pattern" "\${REFERENCES}" > HA_CDS_ref.fasta || true
    HA_INPUT_FILE="samples/${sample_dir}/segments/HA/${sample_id}_HA.fasta"
    if [[ -s HA_CDS_ref.fasta && -s "$\HA_INPUT_FILE" ]]; then
        cat HA_CDS_ref.fasta "\$HA_INPUT_FILE" > mafft_input_HA.fasta
        mafft --auto mafft_input_HA.fasta > "\$OUTFILE"
    else
        echo ">NO_HA_CDS_FOUND" > "\$OUTFILE"
    fi
    """
}