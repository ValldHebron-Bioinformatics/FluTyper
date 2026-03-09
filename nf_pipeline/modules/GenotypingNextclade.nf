process GenotypingNextclade {
    errorStrategy 'ignore'
    
    input:
        tuple val(sample_id), path(input_fasta)
        path(inferred_subtypes)

    output:
    path("genotyping_output.csv")

    script:
    """
    # # Primer identifiquem el subtipus H de la mostra actual consultant el fitxer de subtipus inferits
    h_subtype=\$(grep -m 1 -E "^${sample_id}\b" "${inferred_subtypes}" | cut -f2 | grep -oE 'H[0-9]+' | head -n 1 || true)

    # Filtrem el fitxer de subtipus per extreure els IDs que coincideixen amb el subtipus global demanat
    # Fem servir una variable personalitzada en awk per evitar conflictes amb paraules reservades del sistema
    awk -F'\t' -v mysub="${params.subtype}" '\$2 ~ "^"mysub {print \$1}' "${inferred_subtypes}" > target_ids.txt
    
    # Preparem el fitxer FASTA temporal per a les seqüències HA que passaran a l'anàlisi
    touch filtered_HA_subtype.fasta
    
    # Si hem trobat IDs coincidents, extraiem les seves seqüències HA del FASTA d'entrada
    # La regex gestiona tant el format d'ID amb [|_]
    if [ -s target_ids.txt ]; then
        for id in \$(cat target_ids.txt); do
            seqkit grep -r -p "^\${id}[|_]HA[|_]" "${input_fasta}" >> filtered_HA_subtype.fasta
        done
    fi

    # Determinem quin dataset de Nextclade s'ha d'utilitzar en funció del subtipus H identificat
    if [[ "${params.subtype}" == "H5" ]]; then
        DATASET_NAME='community/moncla-lab/iav-h5/ha/2.3.4.4'
    elif [[ "${params.subtype}" == "H7" ]]; then
        DATASET_NAME="TO BE DECIDED"
    elif [[ "${params.subtype}" == "H9" ]]; then
        DATASET_NAME="TO BE DECIDED"
    fi

    # Si tenim una seqüència vàlida i un dataset assignat, procedim amb l'anàlisi
    if [[ -n "\${DATASET_NAME}" && "\${DATASET_NAME}" != "COMING SOON" ]]; then

        # Verifiquem la disponibilitat del dataset localment i el copiem a l'entorn de treball
        LOCAL_DATASET="${params.workDir}/../docs/nextclade_dataset"
        if [ ! -f "\${LOCAL_DATASET}/pathogen.json" ]; then
            nextclade dataset get --name "\${DATASET_NAME}" --output-dir "\${LOCAL_DATASET}"
        fi
        cp -r "\${LOCAL_DATASET}" ./nextclade_dataset

        # Executem Nextclade per generar el fitxer CSV de resultats bruts
        nextclade run \
            --input-dataset nextclade_dataset \
            --output-csv nextclade_results_${sample_id}.csv \
            filtered_HA_subtype.fasta

        # Utilitzem Python per processar el CSV, ja que gestiona correctament els delimitadors interns
        # També netegem l'ID de la seqüència per mantenir només la referència principal
        python3 - <<'EOF'
import csv
import re
with open("nextclade_results_${sample_id}.csv", newline='') as infile, open("genotyping_output.csv", 'w', newline='') as outfile:
    reader = csv.DictReader(infile, delimiter=';')  # Llegim el CSV amb delimitador ;
    writer = csv.writer(outfile)  # Preparem el writer per escriure el CSV de sortida
    writer.writerow(["seqID", "predicted_clade", "qc.overallStatus", "qc.overallScore"])  # Capçalera del CSV
    for row in reader:
        seqid = row.get("seqName", "")  # Agafem el nom de la seqüència
        # Extraiem només l'ID principal (abans de _ o |)
        match = re.match(r"([^_|]+)", seqid)
        only_id = match.group(1) if match else seqid
        writer.writerow([
            only_id,  # ID net
            row.get("clade", ""),  # Clade predit
            row.get("qc.overallStatus", ""),  # Estat QC
            row.get("qc.overallScore", "")    # Puntuació QC
        ])
EOF
    else
        echo "Skipping genotyping: No matching or supported subtype found for \${h_subtype}"
    fi
    """
}