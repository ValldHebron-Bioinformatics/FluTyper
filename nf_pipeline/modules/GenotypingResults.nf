#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GenotypingResults {
    // This process compiles the genotyping results from Nextclade and Genin2 into a final CSV report for each sample.
    
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(h_tag), val(n_tag), val(pathotype), path(csv_path), path(genin_path)
    path("datasets/*")

    output:
    tuple val(sample_id), path("final_genotyping_results_${sample_id}.csv"), emit: results
    tuple val(sample_id), path("GRerrors.log"), optional: true, emit: errors

    script:
    """
#!/usr/bin/env python3
import os, csv, subprocess
# Determine dataset name and version based on the H tag
d_path = "datasets/nextclade_${h_tag}_dataset"
d_name = "nextclade_${h_tag}_dataset" if os.path.isdir(d_path) else "-"

subtype_val = "${h_tag}${n_tag}(${pathotype})" if "${pathotype}" != "" else "${h_tag}${n_tag}"
if "${params.protocol}" == "HUMAN":
    if subtype_val == "H1N1":
        subtype_val = "A(H1N1)pdm09"
    else:
        subtype_val = "A(${h_tag}${n_tag})"

# Prepare the data dictionary with default values
data = {
    "SampleID": "${sample_id}",
    "Subtype": subtype_val,
    "Dataset": d_name,
    "Dataset Version": "-", 
    "Clade": "-",
    "Genin Version": "-",
    "Genotype": "-", 
    "Sub-genotype": "-",
    "Notes": "-"
}

# Extract version from CHANGELOG.md if the dataset exists
if d_name != "-":
    changelog_file = f"{d_path}/CHANGELOG.md"
    if os.path.isfile(changelog_file):
        with open(changelog_file, 'r') as f:
            first_line = f.readline()
            if first_line.startswith('##'):
                data["Dataset Version"] = first_line.replace('##', '').strip()
    else:
        with open("GRerrors.log", 'a') as log_f:
            log_f.write(f"GenotypingResults: CHANGELOG.md not found for dataset {d_name}, cannot extract version.\\n")

if os.path.isfile("${csv_path}"):
    best_row = None
    min_score = float('inf')
    
    with open("${csv_path}", 'r') as f:
        for row in csv.DictReader(f, delimiter=';'):
            score_str = row.get('qc.overallScore', '-')
            try:
                score = float(score_str)
                if score < min_score:
                    min_score, best_row = score, row
            except ValueError:
                with open("GRerrors.log", 'a') as log_f:
                    log_f.write(f"GenotypingResults: Invalid qc.overallScore '{score_str}' for sample ${sample_id}, skipping this row.\\n")
                pass
    
    if best_row:
        data["Clade"] = best_row.get('clade', 'unclassified')

# Genin data extaction (if available)
genin_str = "${genin_path}".replace("[", "").replace("]", "").strip()
if genin_str and os.path.isfile(genin_str):
    with open(genin_str, 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            data["Genotype"] = row.get("Genotype", "-")
            data["Sub-genotype"] = row.get("Sub-genotype", "-")
            data["Notes"] = row.get("Notes", "-")
            break # Amb la primera fila de genin ja en tenim prou
    
    # Only set Genin version if genotype was actually found
    if data["Genotype"] != "-":
        try:
            res = subprocess.run(['genin2', '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if res.returncode == 0:
                data["Genin Version"] = res.stdout.split(',')[1].replace('version', '').strip()
        except Exception as e:
            with open("GRerrors.log", 'a') as log_f:
                log_f.write(f"GenotypingResults: Could not determine Genin2 version. Error: {e}\\n")

with open("final_genotyping_results_${sample_id}.csv", 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=data.keys())
    writer.writeheader() 
    writer.writerow(data)
    """
}
