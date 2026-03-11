#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GenotypingResults {
    input:
    tuple val(sample_id), val(h_tag), val(n_tag), val(pathotype), path(csv_path)
    path("datasets/*") 

    output:
    path "final_genotyping_results_${sample_id}.csv"

    script:
    """
    #!/usr/bin/env python3
    import os, csv

    csv_path = "${csv_path}"
    # Creem la ruta esperada de la base de dades, dins el work dir, per comprovar si existeix
    d_path = f"datasets/nextclade_{'${h_tag}'}_dataset"
    
    # Comprovem si la carpeta existeix realment
    d_name = f"nextclade_{'${h_tag}'}_dataset" if os.path.isdir(d_path) else "-"
    
    # Definim els valors per defecte al principi per si falla alguna cosa
    data = {
        "SampleID": "${sample_id}", 
        "Subtype": "${h_tag}${n_tag}", 
        "Dataset": d_name, 
        "Version": "-", 
        "Clade": "-", 
        "qc.status": "-", 
        "qc.score": "-"
    }

    # Actualitzem els valors només si la base de dades existeix
    if d_name != "-":
        # Apuntem directament al fitxer CHANGELOG.md per extreure la versió
        changelog_file = f"{d_path}/CHANGELOG.md"
        if os.path.isfile(changelog_file):
            # Obrim el fitxer i busquem la versió
            with open(changelog_file, 'r') as f:
                for line in f:
                    if line.startswith('##'):
                        # Extreiem la versió de la línia, traiem els '##' i els espais sobrants
                        data["Version"] = line.replace('##', '').strip()
                        break

        if os.path.isfile(csv_path):
            best_row = None
            min_score = float('inf') # ho posem a infinit perquè qualsevol score real serà més petit
            # Llegim el fitxer CSV línia per línia, utilitzem csv per accedir a les columnes pel nom
            for row in csv.DictReader(open(csv_path), delimiter=';'):
                score_str = row.get('qc.overallScore', row.get('qc.score', '-'))
                try:
                    score = float(score_str)
                except ValueError:
                    score = float('inf') # Error handling
                # Ens quedem només amb la fila que tingui el score més baix
                if score < min_score:
                    min_score = score
                    best_row = row
            # Actualitzem les dades de clade i QC només si hem trobat una fila vàlida amb un score numèric
            if best_row:
                data["Clade"] = best_row.get('clade', 'unclassified')
                data["qc.status"] = best_row.get('qc.overallStatus', best_row.get('qc.status', '-'))
                data["qc.score"] = best_row.get('qc.overallScore', best_row.get('qc.score', '-'))
    # Escrivim el header i la fila amb el millor resultat en el nou CSV
    with open("final_genotyping_results_${sample_id}.csv", 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=data.keys()) # Escrivim el header només una vegada
        writer.writeheader() 
        writer.writerow(data)
    """
}