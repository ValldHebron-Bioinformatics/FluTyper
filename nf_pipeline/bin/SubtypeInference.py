#! /usr/bin/env python3

import csv
import re
import sys

# 1. Quedar-se amb el millor score
# Aquí guardem la millor peça trobada fins ara per a cada (mostra, segment)
best_match = {}
missing_score_ha_na = []

with open('minimizers_results.tsv', 'r') as f:
    lector = csv.DictReader(f, delimiter='\t')
    for line in lector:
        seqName = (line.get('seqName') or '').strip()
        score_raw = (line.get('score') or '').strip()
        match = (line.get('dataset') or '').strip()

        if not seqName:
            continue

        # Extraiem mostra i segment poden estar separades per '|' o '_'
        parts = re.split(r'[|_]', seqName)
        if len(parts) >= 2:
            id = parts[0]
            segment = parts[1].upper()

            if segment in ["HA", "NA"] and not score_raw:
                missing_score_ha_na.append(seqName)
                continue

            try:
                score = float(score_raw)
            except ValueError:
                if segment in ["HA", "NA"]:
                    missing_score_ha_na.append(seqName)
                continue

            if not match:
                continue

            if segment in ["HA", "NA"]:
                key = (id, segment)
                # Només guardem si el score és més alt que l'anterior
                if key not in best_match or score > best_match[key]['score']:
                    best_match[key] = {'dataset': match, 'score': score}

if missing_score_ha_na:
    unique_missing = sorted(set(missing_score_ha_na))
    print(f"[SubtypeInference] WARNING: Missing/invalid score for {len(unique_missing)} HA/NA sequences", file=sys.stderr)
    for seq_name in unique_missing:
        print(f"[SubtypeInference] Missing score: {seq_name}", file=sys.stderr)

# 2. Unió de H + N per a cada mostra

inferred_subtypes = {}

for (id, segment), dades in best_match.items():
    if id not in inferred_subtypes:
        inferred_subtypes[id] = {'ha': '', 'ha_s': 0, 'na': '', 'na_s': 0}
    
    if segment == 'HA':
        inferred_subtypes[id]['ha'] = dades['dataset']
        inferred_subtypes[id]['ha_s'] = dades['score']
    else:
        inferred_subtypes[id]['na'] = dades['dataset']
        inferred_subtypes[id]['na_s'] = dades['score']

# 3. Escriptura del resultat amb el subtipus inferit
with open('inferred_subtypes.tsv', 'w') as f_out:
    f_out.write("seqName\tha_match\tha_score\tna_match\tna_score\tinferred_subtype\n")
    
    for id in sorted(inferred_subtypes.keys()):
        d = inferred_subtypes[id]
        
        # Busquem el número de H i N (ex: H5 i N1)
        h_num = re.search(r'H[0-9]+', d['ha'])
        n_num = re.search(r'N[0-9]+', d['na'])
        
        h_text = h_num.group(0) if h_num else ""
        n_text = n_num.group(0) if n_num else ""
        subtype = h_text + n_text if (h_text and n_text) else "Incomplete"
        
        f_out.write(f"{id}\t{d['ha']}\t{d['ha_s']}\t{d['na']}\t{d['na_s']}\t{subtype}\n")