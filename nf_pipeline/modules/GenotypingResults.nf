#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GenotypingResults {
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path("nextclade_results_${sample_id}.csv"), val(h_tag), val(n_tag)
    path("*/nextclade_*_dataset")

    output:
    path("final_genotyping_results.csv")

    script:
    """
    if [[ ${h_tag} == "H5" ]]; then
        dataset_dir = \$(find . -type d -name "nextclade_H5_dataset" | head -n 1)
        version = \$(grep '^## ' "\${dataset_dir}/CHANGELOG.md" | head -n 1 | cut -d ' ' -f 2)
    elif [[ ${h_tag} == "H7" ]]; then
        dataset_dir = \$(find . -type d -name "nextclade_H7_dataset" | head -n 1)
        version = \$(grep '^## ' "\${dataset_dir}/CHANGELOG.md" | head -n 1 | cut -d ' ' -f 2)
    elif [[ ${h_tag} == "H9" ]]; then
        dataset_dir = \$(find . -type d -name "nextclade_H9_dataset" | head -n 1)
        version = \$(grep '^## ' "\${dataset_dir}/CHANGELOG.md" | head -n 1 | cut -d ' ' -f 2)
    else
        dataset_dir = ""
        version =""
    fi
    subtype = "${h_tag}${n_tag}"
    clade = \$(cat "nextclade_results_${sample_id}.csv" | tail -n +2 | cut -d ',' -f5 | head -n 1)
    echo "SampleID,Subtype,Dataset,Version,Clade,qc.status,qc.score" > final_genotyping_results.csv
    echo "\${sample_id},\${subtype},\${dataset_dir},\${version},\${clade}" >> final_genotyping_results.csv



    """
}