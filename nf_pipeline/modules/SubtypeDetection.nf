#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SubtypeDetection {
    errorStrategy 'ignore'

    input:
    tuple val(params.sample), path(params.dirSample)

    output:
    path("inferred_subtypes.tsv")

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
        printf 'id\tha_match\tha_score\tna_match\tna_score\tinferred_subtype\n' > inferred_subtypes.tsv
        exit 0
    fi

    nextclade sort -m "\${minimizer_index}" -r minimizers_results.tsv "\${input_fasta}"

        # 1) Tria millor match (score més alt) per espècie i segment HA/NA
        awk -F '\t' '
            BEGIN { OFS = "\t" }
            NR == 1 {next}
            {
                seq_name = \$2
                dataset  = \$3
                score    = \$4 + 0
                species = ""
                segment = ""

                if (index(seq_name, "|") > 0) {
                    split(seq_name, a, "[|]")
                    species = a[1]
                    segment = toupper(a[2])
                } else if (seq_name ~ /^H[0-9]+N[0-9]+_(HA|NA)_/) {
                    n = split(seq_name, a, "_")
                    segment = toupper(a[2])
                    species = substr(seq_name, length(a[1]) + length(a[2]) + 3)
                }

                if (segment != "HA" && segment != "NA") next
                key = species SUBSEP segment

                if (!(key in best_score) || score > best_score[key]) {
                    best_score[key] = score
                    best_dataset[key] = dataset
                    best_species[key] = species
                    best_segment[key] = segment
                }
            }
            END {
                for (k in best_score) {
                    print best_species[k], best_segment[k], best_dataset[k], best_score[k]
                }
            }
        ' minimizers_results.tsv | sort -k1,1 -k2,2 > best_segments.tsv

        # 2) Combina H (de HA) + N (de NA) i escriu subtipus inferit
        awk -F '\t' 'BEGIN {
            OFS = "\t"
            print "id", "ha_match", "ha_score", "na_match", "na_score", "inferred_subtype"
        }
        {
            species = \$1
            segment = \$2
            dataset = \$3
            score   = \$4

            if (!(species in seen)) {
                seen[species] = 1
                order[++n] = species
            }

            if (segment == "HA") {
                ha_match[species] = dataset
                ha_score[species] = score
            } else if (segment == "NA") {
                na_match[species] = dataset
                na_score[species] = score
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                species = order[i]
                hm = ha_match[species]
                hs = ha_score[species]
                nm = na_match[species]
                ns = na_score[species]

                h = ""
                nseg = ""
                inferred = ""

                if (match(hm, /H[0-9]+/, m1)) h = m1[0]
                if (match(nm, /N[0-9]+/, m2)) nseg = m2[0]
                if (h != "" && nseg != "") inferred = h nseg

                print species, hm, hs, nm, ns, inferred
            }
        }' best_segments.tsv > inferred_subtypes.tsv

    """
}