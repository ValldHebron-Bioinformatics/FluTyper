#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GenotypingNextclade {
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(ha_fasta), val(h_tag), val(n_tag), val(pathotype), val(dataset_dir)
    output:
    path("nextclade_results_${sample_id}.csv")
    script:
    """
    if [ ! -d "${dataset_dir}" ]; then
        echo "[WARNING] No valid dataset_dir for ${sample_id}, skipping."
        touch nextclade_results_${sample_id}.csv
        exit 0
    fi
    if [[ ${h_tag} == "H5" ]]; then
        dataset_dir="${dataset_dir}/H5/nextclade_H5_dataset"
    elif [[ ${h_tag} == "H7" ]]; then
        ##dataset_dir="${dataset_dir}/H7/nextclade_H7_dataset"
        # This is just to keep the same header, once we have a real H7 dataset we can update this part
        echo "index;seqName;clade;polybasic_cleavage_site;cleavage_site;glycosylation;qc.overallScore;qc.overallStatus;totalSubstitutions;totalDeletions;totalInsertions;totalFrameShifts;totalMissing;totalNonACGTNs;totalAminoacidSubstitutions;totalAminoacidDeletions;totalAminoacidInsertions;totalUnknownAa;alignmentScore;alignmentStart;alignmentEnd;coverage;cdsCoverage;isReverseComplement;substitutions;deletions;insertions;frameShifts;aaSubstitutions;aaDeletions;aaInsertions;privateNucMutations.reversionSubstitutions;privateNucMutations.labeledSubstitutions;privateNucMutations.unlabeledSubstitutions;privateNucMutations.totalReversionSubstitutions;privateNucMutations.totalLabeledSubstitutions;privateNucMutations.totalUnlabeledSubstitutions;privateNucMutations.totalPrivateSubstitutions;privateAaMutations.reversionSubstitutions;privateAaMutations.labeledSubstitutions;privateAaMutations.unlabeledSubstitutions;privateAaMutations.totalReversionSubstitutions;privateAaMutations.totalLabeledSubstitutions;privateAaMutations.totalUnlabeledSubstitutions;privateAaMutations.totalPrivateSubstitutions;missing;founderMuts['clade'].nodeName;founderMuts['clade'].substitutions;founderMuts['clade'].deletions;founderMuts['clade'].aaSubstitutions;founderMuts['clade'].aaDeletions;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].nodeName;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].substitutions;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].deletions;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].aaSubstitutions;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].aaDeletions;unknownAaRanges;nonACGTNs;qc.missingData.missingDataThreshold;qc.missingData.score;qc.missingData.status;qc.missingData.totalMissing;qc.mixedSites.mixedSitesThreshold;qc.mixedSites.score;qc.mixedSites.status;qc.mixedSites.totalMixedSites;qc.privateMutations.cutoff;qc.privateMutations.excess;qc.privateMutations.score;qc.privateMutations.status;qc.privateMutations.total;qc.snpClusters.clusteredSNPs;qc.snpClusters.score;qc.snpClusters.status;qc.snpClusters.totalSNPs;qc.frameShifts.frameShifts;qc.frameShifts.totalFrameShifts;qc.frameShifts.frameShiftsIgnored;qc.frameShifts.totalFrameShiftsIgnored;qc.frameShifts.score;qc.frameShifts.status;qc.stopCodons.stopCodons;qc.stopCodons.totalStopCodons;qc.stopCodons.score;qc.stopCodons.status;totalPcrPrimerChanges;pcrPrimerChanges;failedCdses;warnings;errors" > nextclade_results_${sample_id}.csv
        exit 0
    elif [[ ${h_tag} == "H9" ]]; then
        ##dataset_dir="${dataset_dir}/H9/nextclade_H9_dataset"
         # Same here, just to keep the same header until we have a real H9 dataset
        echo "index;seqName;clade;polybasic_cleavage_site;cleavage_site;glycosylation;qc.overallScore;qc.overallStatus;totalSubstitutions;totalDeletions;totalInsertions;totalFrameShifts;totalMissing;totalNonACGTNs;totalAminoacidSubstitutions;totalAminoacidDeletions;totalAminoacidInsertions;totalUnknownAa;alignmentScore;alignmentStart;alignmentEnd;coverage;cdsCoverage;isReverseComplement;substitutions;deletions;insertions;frameShifts;aaSubstitutions;aaDeletions;aaInsertions;privateNucMutations.reversionSubstitutions;privateNucMutations.labeledSubstitutions;privateNucMutations.unlabeledSubstitutions;privateNucMutations.totalReversionSubstitutions;privateNucMutations.totalLabeledSubstitutions;privateNucMutations.totalUnlabeledSubstitutions;privateNucMutations.totalPrivateSubstitutions;privateAaMutations.reversionSubstitutions;privateAaMutations.labeledSubstitutions;privateAaMutations.unlabeledSubstitutions;privateAaMutations.totalReversionSubstitutions;privateAaMutations.totalLabeledSubstitutions;privateAaMutations.totalUnlabeledSubstitutions;privateAaMutations.totalPrivateSubstitutions;missing;founderMuts['clade'].nodeName;founderMuts['clade'].substitutions;founderMuts['clade'].deletions;founderMuts['clade'].aaSubstitutions;founderMuts['clade'].aaDeletions;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].nodeName;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].substitutions;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].deletions;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].aaSubstitutions;relativeMutations['A/American_wigeon/North_Carolina/AH0182517/2022'].aaDeletions;unknownAaRanges;nonACGTNs;qc.missingData.missingDataThreshold;qc.missingData.score;qc.missingData.status;qc.missingData.totalMissing;qc.mixedSites.mixedSitesThreshold;qc.mixedSites.score;qc.mixedSites.status;qc.mixedSites.totalMixedSites;qc.privateMutations.cutoff;qc.privateMutations.excess;qc.privateMutations.score;qc.privateMutations.status;qc.privateMutations.total;qc.snpClusters.clusteredSNPs;qc.snpClusters.score;qc.snpClusters.status;qc.snpClusters.totalSNPs;qc.frameShifts.frameShifts;qc.frameShifts.totalFrameShifts;qc.frameShifts.frameShiftsIgnored;qc.frameShifts.totalFrameShiftsIgnored;qc.frameShifts.score;qc.frameShifts.status;qc.stopCodons.stopCodons;qc.stopCodons.totalStopCodons;qc.stopCodons.score;qc.stopCodons.status;totalPcrPrimerChanges;pcrPrimerChanges;failedCdses;warnings;errors" > nextclade_results_${sample_id}.csv
        exit 0
    else
        echo "No valid H subtype found for genotyping: ${h_tag}"
        touch nextclade_results_${sample_id}.csv
        exit 0
    fi
    nextclade run \
        --input-dataset "${dataset_dir}" \
        --output-csv nextclade_results_${sample_id}.csv \
        "${ha_fasta}"


    """
}