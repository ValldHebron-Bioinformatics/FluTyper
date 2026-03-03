#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SubtypeDetection {
    errorStrategy 'ignore'

    input:
    tuple val(params.sample), path(params.dirSample)

    output:
    path("minimizers_results.tsv"), emit: minimizers
    path("inferred_subtypes.tsv"), emit: inferred

    script:
    """
    input_fasta="${params.dirSample}/${params.sample}"

    if [[ "${params.protocol}" == "AVIAN" ]]; then
        echo "Subtype detection for AVIAN protocol"
        minimizer_index="${params.protocols.AVIAN.resources}/Avian_minimizers.json"
    elif [[ "${params.protocol}" == "SWINE" ]]; then
        echo "Subtype detection for SWINE protocol"
        minimizer_index="${params.protocols.SWINE.resources}/Swine_minimizers.json"
    else
        echo "No valid protocol specified for subtype detection: ${params.protocol}"
        : > minimizers_results.tsv
        printf 'species\tha_match\tha_score\tna_match\tna_score\tinferred_subtype\n' > inferred_subtypes.tsv
        exit 0
    fi

    nextclade sort -m "\${minimizer_index}" -r minimizers_results.tsv "\${input_fasta}"

    python3 - <<'PY'
import csv
import re

min_file = 'minimizers_results.tsv'
out_file = 'inferred_subtypes.tsv'

subtype_re = re.compile(r'^(H[0-9]+N[0-9]+)_')
h_re = re.compile(r'(H[0-9]+)')
n_re = re.compile(r'(N[0-9]+)')

best = {}

def parse_seq_name(seq_name: str):
    parts = seq_name.split('|')
    if len(parts) >= 2:
        group_id = parts[0].strip()
        segment = parts[1].strip().upper()
        if segment in {'HA', 'NA'}:
            return group_id, segment

    m = re.match(r'^(H[0-9]+N[0-9]+)_(HA|NA)_(.+)', seq_name)
    if m:
        segment = m.group(2)
        group_id = m.group(3).strip()
        return group_id, segment

    return None, None

with open(min_file, newline='', encoding='utf-8') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        seq_name = (row.get('seqName') or '').strip()
        dataset = (row.get('dataset') or '').strip()
        score_raw = (row.get('score') or '').strip()
        if not seq_name or not dataset:
            continue

        species, segment = parse_seq_name(seq_name)
        if not species or not segment:
            continue
        if segment not in {'HA', 'NA'}:
            continue

        try:
            score = float(score_raw)
        except ValueError:
            score = float('-inf')

        key = (species, segment)
        current = best.get(key)
        if current is None or score > current['score']:
            best[key] = {
                'dataset': dataset,
                'score': score,
            }

species_ids = sorted({s for s, _ in best.keys()})

with open(out_file, 'w', newline='', encoding='utf-8') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['species', 'ha_match', 'ha_score', 'na_match', 'na_score', 'inferred_subtype'])

    for species in species_ids:
        ha = best.get((species, 'HA'))
        na = best.get((species, 'NA'))

        ha_match = ha['dataset'] if ha else ''
        na_match = na['dataset'] if na else ''
        ha_score = '' if not ha or ha['score'] == float('-inf') else f"{ha['score']:.6f}"
        na_score = '' if not na or na['score'] == float('-inf') else f"{na['score']:.6f}"

        ha_sub = subtype_re.match(ha_match).group(1) if ha_match and subtype_re.match(ha_match) else ''
        na_sub = subtype_re.match(na_match).group(1) if na_match and subtype_re.match(na_match) else ''

        h_part = h_re.search(ha_sub).group(1) if ha_sub and h_re.search(ha_sub) else ''
        n_part = n_re.search(na_sub).group(1) if na_sub and n_re.search(na_sub) else ''

        inferred = f"{h_part}{n_part}" if (h_part and n_part) else ''
        writer.writerow([species, ha_match, ha_score, na_match, na_score, inferred])
PY

    """
}