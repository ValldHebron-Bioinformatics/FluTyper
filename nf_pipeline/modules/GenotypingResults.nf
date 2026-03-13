#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GenotypingResults {
    input:
    tuple val(sample_id), val(h_tag), val(n_tag), val(pathotype), path(csv_path)
    path("datasets/*") // datasets/ directory containing the downloaded datasets in the work dir.

    output:
    path "final_genotyping_results_${sample_id}.csv"

    script:
    """
    #!/usr/bin/env python3
    import os, csv

    # Determine dataset name and version based on the H tag
    d_path = "datasets/nextclade_${h_tag}_dataset"
    d_name = "nextclade_${h_tag}_dataset" if os.path.isdir(d_path) else "-"
    
    subtype_val = "${h_tag}${n_tag}(${pathotype})" if "${h_tag}" in ["H5", "H7", "H9"] and "${pathotype}" else "${h_tag}${n_tag}"

    # Prepare the data dictionary with default values
    data = {
        "SampleID": "${sample_id}", 
        "Subtype": subtype_val, 
        "Dataset": d_name, 
        "Version": "-", 
        "Clade": "-", 
        "qc.status": "-", 
        "qc.score": "-"
    }
    # Extract version from CHANGELOG.md if the dataset exists
    if d_name != "-":
        changelog_file = f"{d_path}/CHANGELOG.md"
        if os.path.isfile(changelog_file):
            with open(changelog_file, 'r') as f:
                first_line = f.readline()
                if first_line.startswith('##'):
                    data["Version"] = first_line.replace('##', '').strip()

        if os.path.isfile("${csv_path}"):
            best_row = None
            min_score = float('inf') # We want to find the best (lowest) score among the rows in the CSV
            
            with open("${csv_path}", 'r') as f:
                for row in csv.DictReader(f, delimiter=';'):
                    score_str = row.get('qc.overallScore', '-')
                    try:
                        score = float(score_str)
                        if score < min_score:
                            min_score, best_row = score, row
                    except ValueError:
                        pass
            
            if best_row:
                data["Clade"] = best_row.get('clade', 'unclassified')
                data["qc.status"] = best_row.get('qc.overallStatus', '-')
                data["qc.score"] = best_row.get('qc.overallScore', '-')

    # Set up the CSV writer using the dictionary keys as column headers
    # This ensures that data is always mapped to the correct column name
    with open("final_genotyping_results_${sample_id}.csv", 'w', newline='') as f:
        # DictWriter uses 'fieldnames' to know which keys to look for in the data dictionary https://docs.python.org/3/library/csv.html
        writer = csv.DictWriter(f, fieldnames=data.keys())
        # Write the header row using the keys provided in fieldnames
        writer.writeheader() 
        # Write the data row; the writer matches the dictionary values to the correct headers
        writer.writerow(data)
    """
}