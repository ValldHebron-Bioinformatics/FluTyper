#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GetDatasets {
    errorStrategy 'ignore'
    debug true

    input:
    path(inferred_subtypes)


    output:
    path("*/nextclade_*_dataset", emit: datasets, optional: true)


    script:
    """
    # Extract unique H tags from the inferred subtypes CSV
    h_tags=\$(grep -oE 'H[0-9]+' "${inferred_subtypes}" || true) 
    
    if [ -z "\$h_tags" ]; then
        echo "GetDatasets: No H tags found in ${inferred_subtypes}."
        exit 0
    fi

    if [ "${params.protocol}" = "HUMAN" ]; then
        echo "\$h_tags" | sort -u | while read -r tag; do
            case "\$tag" in
                H1)
                    DATASET_NAME='flu_h1n1pdm_ha'
                    mkdir -p H1/nextclade_H1_dataset
                    nextclade dataset get --name "\${DATASET_NAME}" --output-dir H1/nextclade_H1_dataset
                    ;;
                H3)
                    DATASET_NAME='flu_h3n2_ha'
                    mkdir -p H3/nextclade_H3_dataset
                    nextclade dataset get --name "\${DATASET_NAME}" --output-dir H3/nextclade_H3_dataset
                    ;;
                *)
                    echo "GetDatasets: No valid H tag found for HUMAN protocol dataset retrieval: \$tag" 
                    ;;
            esac
        done
    else
        echo "\$h_tags" | sort -u | while read -r tag; do
        case "\$tag" in
            H5)
                DATASET_NAME='community/moncla-lab/iav-h5/ha/2.3.4.4'
                mkdir -p H5
                nextclade dataset get --name "\${DATASET_NAME}" --output-dir H5/nextclade_H5_dataset
                ;;
            H7)
                DATASET_NAME='TO_BE_DECIDED_H7' # No H7 dataset yet
                mkdir -p H7/nextclade_H7_dataset
                echo "GetDatasets: H7 dataset retrieval is currently under development."
                #nextclade dataset get --name "\${DATASET_NAME}" --output-dir H7/nextclade_H7_dataset
                ;;
            H9)
                DATASET_NAME='TO_BE_DECIDED_H9' # No H9 dataset yet
                mkdir -p H9/nextclade_H9_dataset
                echo "GetDatasets: H9 dataset retrieval is currently under development."
                #nextclade dataset get --name "\${DATASET_NAME}" --output-dir H9/nextclade_H9_dataset
                ;;
            *)
                echo "GetDatasets: No valid H tag found for dataset retrieval: \$tag" 
                ;;
        esac
    done
    fi
    """
}