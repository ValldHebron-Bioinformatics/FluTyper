#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process OrganizeBySpecies {
    errorStrategy 'ignore' // Ignora errors i continua

    input:
    tuple val(params.sample), path(params.dirSample)

    output:
    path("sequences")

    script:
    """
    # Crea la carpeta "pare"
    mkdir -p sequences

    awk '
      /^>/ {
        # Extreu espècie i proteïna del header
        split(substr(\$0, 2), a, "_")
        species = a[1]
        protein = a[2]

        # Rutes de sortida
        dir = "sequences/" species
        prot_dir = dir "/proteins"
        system("mkdir -p " prot_dir)

        cons_file = dir "/consensus.fasta"
        prot_file = prot_dir "/" protein ".fasta"
      }
      {
        # Escriu a consensus i per proteïna
        if (cons_file != "") {
          print \$0 >> cons_file
          print \$0 >> prot_file
        }
      }
    ' ${params.dirSample}/${params.sample}
    """
}