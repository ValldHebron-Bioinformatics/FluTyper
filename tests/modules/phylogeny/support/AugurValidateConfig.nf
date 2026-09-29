#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AugurValidateConfig {
    // TEST SUPPORT ONLY (not part of the pipeline): runs `augur validate auspice-config-v2` (augur 34.1.4, schema
    // augur/data/schema-auspice-config-v2.json) on an Auspice config inside the same Nextstrain env the pipeline's
    // augur processes use, so a config edit that breaks the schema fails the test suite (plan Phase 2 acceptance).
    // moduleDir is tests/modules/phylogeny/support, hence four levels up to the repo root (D14: moduleDir, not projectDir).
    conda "${moduleDir}/../../../../envs/nextstrain.yaml"

    input:
    path(config)

    output:
    path("validate.log"), emit: log

    script:
    """
    augur validate auspice-config-v2 "${config}" > validate.log 2>&1 || { cat validate.log >&2; exit 1; }
    """
}
