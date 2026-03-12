#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GetDatasets {
    errorStrategy 'ignore'

    input:
    path(inferred_subtypes)


    output:
    path("*/nextclade_*_dataset")
    

    script:
    """
    # Extract unique H tags from the inferred subtypes CSV, 
    # This way to have only 1 process instead of 1 per sample giving the sample id and h_tag tuple
    h_tags=\$(grep -oE 'H[0-9]+' "${inferred_subtypes}") 
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
                
                #nextclade dataset get --name "\${DATASET_NAME}" --output-dir H7/nextclade_H7_dataset
                ;;
            H9)
                DATASET_NAME='TO_BE_DECIDED_H9' # No H9 dataset yet
                mkdir -p H9/nextclade_H9_dataset
                #nextclade dataset get --name "\${DATASET_NAME}" --output-dir H9/nextclade_H9_dataset
                ;;
            *)
                echo "Skipping undefined subtype: \$tag"
                ;;
        esac
    done

    """
}