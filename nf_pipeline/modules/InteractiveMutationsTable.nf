#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process InteractiveMutationsTable {

    input:
    path excel_file

    output:
    path "Interactive_Mutations_Table.html", emit: table

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import json
    import re

    input_excel = "${excel_file}"
    output_html = "Interactive_Mutations_Table.html"

    # Llegim la primera fulla de l'Excel, evitant valors NaN
    df = pd.read_excel(input_excel, sheet_name=0, keep_default_na=False)
    df.columns = df.columns.str.strip()

    results = {}

    def get_color_class(protein, inferred_subtype, real_subtype, found_in):
        combined_sub = f"{inferred_subtype} {real_subtype}"
        inf_h = set(re.findall(r'H(\\d+)', combined_sub))
        inf_n = set(re.findall(r'N(\\d+)', combined_sub))
        
        found_h = set(re.findall(r'H(\\d+)', str(found_in)))
        found_n = set(re.findall(r'N(\\d+)', str(found_in)))
        
        if protein in ['HA1', 'HA2']:
            if inf_h and (inf_h & found_h):
                return 'green'
            else:
                return 'red'
        elif protein == 'NA':
            if inf_n and (inf_n & found_n):
                return 'green'
            else:
                return 'red'
        else: 
            h_match = bool(inf_h and (inf_h & found_h))
            n_match = bool(inf_n and (inf_n & found_n))
            
            if h_match and n_match:
                return 'green'
            elif h_match or n_match:
                return 'yellow'
            else:
                return 'red'

    # Inicialitzem l'estructura de dades per cada mostra
    for idx, row in df.drop_duplicates(subset=['SAMPLE_ID']).iterrows():
        sample_id = str(row['SAMPLE_ID'])
        results[sample_id] = {
            "id_cresa": sample_id,
            "subtype": str(row.get('SUBTYPE', '')),
            "pb2_mutations": [],
            "ha1_mutations": [],
            "ha2_mutations": [],
            "na_mutations": []
        }

    # Filtrem només les mutacions que són marcadors
    df_markers = df[df.get('MARKER', '') == 'Yes']

    for index, row in df_markers.iterrows():
        sample_id = str(row['SAMPLE_ID'])
        protein = str(row.get('PROTEIN', ''))
        
        inferred_subtype = str(row.get('SUBTYPE', ''))
        real_subtype = str(row.get('REAL SUBTYPE', ''))
        found_in = str(row.get('FOUND_IN', ''))
        
        ref_info = str(row.get('POSITION_REF', ''))
        
        color_class = get_color_class(protein, inferred_subtype, real_subtype, found_in)
        
        mut_obj = {
            "mutation": str(row.get('AA_MUTATION', '')),
            "effect": str(row.get('EFFECT', '')),
            "found_in": found_in,
            "reference": str(row.get('REFERENCE', '')),
            "ref_pos": ref_info,
            "color": color_class
        }
        
        if protein == 'PB2':
            results[sample_id]["pb2_mutations"].append(mut_obj)
        elif protein == 'HA1':
            results[sample_id]["ha1_mutations"].append(mut_obj)
        elif protein == 'HA2':
            results[sample_id]["ha2_mutations"].append(mut_obj)
        elif protein == 'NA':
            results[sample_id]["na_mutations"].append(mut_obj)

    # Convertim a llista per al JSON
    output_data = list(results.values())
    
    # Provem d'ordenar si hi ha una columna d'ordre numèrica, si no, per ID
    try:
        output_data.sort(key=lambda x: int(x["id_cresa"]) if x["id_cresa"].isdigit() else x["id_cresa"])
    except:
        pass

    json_data_string = json.dumps(output_data)

    html_template = f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Interactive Mutations Table</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 40px;
                background-color: #f9f9f9;
            }}
            h2 {{
                margin-bottom: 5px;
            }}
            .subtitle {{
                margin-top: 0;
                margin-bottom: 20px;
                color: #555;
                font-size: 14px;
            }}
            table {{
                border-collapse: collapse;
                width: 100%;
                background-color: white;
                box-shadow: 0 1px 3px rgba(0,0,0,0.1);
                margin-bottom: 30px;
            }}
            th, td {{
                border: 1px solid #e0e0e0;
                padding: 12px;
                text-align: left;
            }}
            th {{
                background-color: #f2f2f2;
                font-weight: bold;
            }}
            .mutation-tag {{
                display: inline-block;
                padding: 4px 8px;
                border-radius: 4px;
                margin: 2px;
                cursor: help;
                font-weight: bold;
                border: 1px solid #bbdefb;
                background: #e3f2fd;
                color: #0d47a1;
            }}
            .mutation-tag.green {{
                background-color: #d4edda;
                color: #155724;
                border-color: #c3e6cb;
            }}
            .mutation-tag.yellow {{
                background-color: #fff3cd;
                color: #856404;
                border-color: #ffeeba;
            }}
            .mutation-tag.red {{
                background-color: #f8d7da;
                color: #721c24;
                border-color: #f5c6cb;
            }}
            .tooltip {{
                position: relative;
                display: inline-block;
            }}
            .tooltip .tooltiptext {{
                visibility: hidden;
                width: 320px;
                background-color: #333;
                color: #fff;
                text-align: left;
                border-radius: 6px;
                padding: 12px;
                position: absolute;
                z-index: 1;
                bottom: 125%;
                top: auto;
                left: 50%;
                margin-left: -160px;
                opacity: 0;
                transition: opacity 0.2s;
                font-size: 0.9em;
                line-height: 1.5;
                box-shadow: 0 4px 6px rgba(0,0,0,0.3);
                font-weight: normal;
            }}
            .tooltip .tooltiptext::after {{
                content: "";
                position: absolute;
                top: 100%;
                bottom: auto;
                left: 50%;
                margin-left: -5px;
                border-width: 5px;
                border-style: solid;
                border-color: #333 transparent transparent transparent;
            }}
            tbody tr:nth-child(-n+3) .tooltip .tooltiptext {{
                top: 130%;
                bottom: auto;
            }}
            tbody tr:nth-child(-n+3) .tooltip .tooltiptext::after {{
                bottom: 100%;
                top: auto;
                border-color: transparent transparent #333 transparent;
            }}
            .tooltip:hover .tooltiptext {{
                visibility: visible;
                opacity: 1;
            }}
            .legend-container {{
                background-color: #ffffff;
                border: 1px solid #e0e0e0;
                border-left: 4px solid #0d47a1;
                padding: 15px 20px;
                border-radius: 4px;
                margin-top: 20px;
                box-shadow: 0 1px 3px rgba(0,0,0,0.05);
                display: flex;
                flex-direction: column;
                gap: 8px;
            }}
            .legend-title {{
                font-weight: bold;
                color: #333;
                margin-bottom: 5px;
                font-size: 16px;
            }}
            .legend-item {{
                line-height: 1.4;
                font-size: 14px;
                color: #444;
            }}
        </style>
    </head>
    <body>
        <h2>Influenza Mutations Interactive Table</h2>
        <p class="subtitle">
            Information about the markers used can be consulted at FluMut: 
            <a href="https://izsvenezie-virology.github.io/FluMut/docs/markers" target="_blank">https://izsvenezie-virology.github.io/FluMut/docs/markers</a>
        </p>
        <table>
            <thead>
                <tr>
                    <th>Sample ID</th>
                    <th>Subtype</th>
                    <th>PB2 Mutations</th>
                    <th>HA1 mutations</th>
                    <th>HA2 mutations</th>
                    <th>NA mutations</th>
                </tr>
            </thead>
            <tbody id="table-body">
            </tbody>
        </table>

        <div class="legend-container">
            <div class="legend-title">Hover Definitions</div>
            <div class="legend-item"><strong>Reference Position (H5N1 numbering):</strong> Position translated via a standardized numbering schema.</div>
            <div class="legend-item"><strong>Effect:</strong> Effects found for that particular mutation.</div>
            <div class="legend-item"><strong>FOUND IN:</strong> In which subtypes those effects were found.</div>
            <div class="legend-item"><strong>REFERENCE:</strong> Articles where those effects are mentioned.</div>
            
            <div class="legend-title" style="margin-top: 15px;">Mutation Colors (FOUND IN vs Sample Subtype)</div>
            <div class="legend-item"><span class="mutation-tag green" style="cursor:default; padding: 2px 6px;">Green</span> Full match: Both H and N match (PB2), only H matches (HA1/HA2), or only N matches (NA).</div>
            <div class="legend-item"><span class="mutation-tag yellow" style="cursor:default; padding: 2px 6px;">Yellow</span> Partial match: Only H or only N matches (applicable to internal proteins like PB2).</div>
            <div class="legend-item"><span class="mutation-tag red" style="cursor:default; padding: 2px 6px;">Red</span> No match: The mutation was originally found in a completely different subtype.</div>
        </div>

        <script>
            const data = {json_data_string};

            function createMutationTags(mutations) {{
                if (mutations.length === 0) return '-';
                return mutations.map(function(m) {{
                    return '<div class="tooltip mutation-tag ' + m.color + '">' + m.mutation + 
                           '<div class="tooltiptext">' + 
                           '<strong>Reference Position (H5N1 numbering):</strong> ' + m.ref_pos + '<br><br>' + 
                           '<strong>Effect:</strong> ' + m.effect + '<br><br>' + 
                           '<strong>FOUND IN:</strong> ' + m.found_in + '<br><br>' + 
                           '<strong>REFERENCE:</strong> ' + m.reference + 
                           '</div></div>';
                }}).join(' ');
            }}

            const tbody = document.getElementById('table-body');
            data.forEach(function(row) {{
                const tr = document.createElement('tr');
                
                const tdId = document.createElement('td');
                tdId.textContent = row.id_cresa;

                const tdSubtype = document.createElement('td');
                tdSubtype.textContent = row.subtype;
                
                const tdPB2 = document.createElement('td');
                tdPB2.innerHTML = createMutationTags(row.pb2_mutations);
                
                const tdHA1 = document.createElement('td');
                tdHA1.innerHTML = createMutationTags(row.ha1_mutations);
                
                const tdHA2 = document.createElement('td');
                tdHA2.innerHTML = createMutationTags(row.ha2_mutations);
                
                const tdNA = document.createElement('td');
                tdNA.innerHTML = createMutationTags(row.na_mutations);
                
                tr.appendChild(tdId);
                tr.appendChild(tdSubtype);
                tr.appendChild(tdPB2);
                tr.appendChild(tdHA1);
                tr.appendChild(tdHA2);
                tr.appendChild(tdNA);
                
                tbody.appendChild(tr);
            }});
        </script>
    </body>
    </html>
    '''

    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html_template)
    """
}