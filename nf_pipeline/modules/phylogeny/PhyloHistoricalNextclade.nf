#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PhyloHistoricalNextclade {
    // --append + --phylogenetics only: a historical candidate sample published by a run before Nextclade results
    // were persisted per-sample (samples/<id>/nextclade_results.csv) has no such file. Rerun Nextclade on its HA
    // segment against the current subtype's dataset instead of leaving it with subclade NA (mirrors the Nextclade
    // call of GenotypingNextclade, minus genin2, which PHYLOGENETICS does not need).
    errorStrategy 'ignore'
    debug true

    input:
    tuple val(sample_id), path(ha_fasta), path(dataset_dir)

    output:
    tuple val(sample_id), path("nextclade_results_${sample_id}.csv"), emit: results
    tuple val(sample_id), path("PHNerrors.log"), optional: true, emit: errors

    script:
    """
    # The declared nextclade_results_${sample_id}.csv output is not optional: if nextclade fails, still create it
    # (empty, like GenotypingNextclade's no-valid-H-subtype branch) so the task's own output binding still succeeds
    # and PHNerrors.log (errorStrategy 'ignore' would otherwise silently discard both) actually reaches the logs.
    nextclade run \
        --input-dataset "${dataset_dir}" \
        --output-csv nextclade_results_${sample_id}.csv \
        "${ha_fasta}" \
        || { echo "PhyloHistoricalNextclade: nextclade failed for historical sample ${sample_id}" >> PHNerrors.log; touch nextclade_results_${sample_id}.csv; }
    """
}

process PhyloHistoricalNextcladeReuse {
    // --append + --phylogenetics only: a historical candidate sample DOES have a persisted samples/<id>/nextclade_results.csv
    // (from GenotypingNextclade), but on disk it is literally named "nextclade_results.csv" (no per-sample suffix),
    // unlike every other CSV PhyloSubclades collects (nextclade_results_<sample_id>.csv - the name the process's own
    // "nextclade_results_*.csv" glob requires to recover the sample_id). Staging the raw file as-is would silently
    // fail that glob and leave the tip's subclade NA even though its real Nextclade result exists. Just rename it.
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(csv, name: "reused.csv")

    output:
    tuple val(sample_id), path("nextclade_results_${sample_id}.csv"), emit: results

    script:
    """
    cp reused.csv nextclade_results_${sample_id}.csv
    """
}
