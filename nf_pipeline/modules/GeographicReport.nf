#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GeographicReport {
    // This process generates an interactive geographic report using Folium and Plotly based on genotyping and metadata files.
    // It is optionally activated via a column in the metadata file (LOCATION), and it can handle both human and avian influenza protocols.
    // It uses coordinates from a provided coordinates file to map towns and provinces, and it creates a dynamic legend and controls for filtering by season, geographic level, and classification.
    errorStrategy 'ignore'
    debug true

    input:
    path(genotyping_file)
    path(metadata_file) 
    path(coordinates_file)

    output:
    path("GeographicReport.html"), emit: geo_report

    script:
    def meta_str = metadata_file ? metadata_file.toString() : ""
    """
    #!/usr/bin/env python3
    import pandas as pd
    import folium
    import json
    import unicodedata
    import re
    import colorsys
    import random
    from plotly.colors import qualitative

    # Color palette for subtypes including human specific labels
    subtype_base_colors = {
        'H1': '#1F77B4', 'A(H1)pdm09': '#1F77B4',
        'H2': '#D62728',
        'H3': '#FF7F0E', 'A(H3)': '#FF7F0E', 
        'H4': '#2CA02C','H5': '#9467BD','H6': '#8C564B', 
        'H7': '#E377C2','H8': '#7F7F7F','H9': '#BCBD22','H10': '#17BECF','H11': '#393B79','H12': '#637939',
        'H13': '#8C6D31','H14': '#843C39','H15': '#7B4173','H16': '#5254A3','H17': '#8CA252','H18': '#BD9E39' 
    }

    def generate_shades(base_hex, n):
        '''
        Generate a monotone sequence of n shades from a base color
        '''
        if n <= 0: return []
        if n == 1: return [base_hex]
        
        clean_hex = str(base_hex).lstrip('#')
        
        try:
            # Convert the hex string into basic Red, Green, and Blue values
            r, g, b = [int(clean_hex[i:i+2], 16) / 255.0 for i in (0, 2, 4)]
        except ValueError:
            return ['#888888'] * n
            
        # Extract the Hue (color identity), Lightness, and Saturation (intensity)
        hue, lightness, saturation = colorsys.rgb_to_hls(r, g, b)
        
        # Create an evenly spaced sequence of light to dark values between 0.1 and 0.9
        step = 0.8 / (n - 1)
        brightness_levels = [0.1 + (step * i) for i in range(n)]
        
        # Shuffle the brightness levels using the base color as a predictable seed
        random.Random(base_hex).shuffle(brightness_levels)
        
        shades = []
        for new_lightness in brightness_levels:
            # Rebuild the color using the new lightness, and convert it back to a hex string
            new_r, new_g, new_b = colorsys.hls_to_rgb(hue, new_lightness, saturation)
            shades.append(f"#{int(new_r * 255):02x}{int(new_g * 255):02x}{int(new_b * 255):02x}")
            
        return shades

    def normalize_str(s):
        '''
        Normalize a string for consistent matching
        '''
        if pd.isna(s): return ""
        s = str(s).strip().lower()
        s = ''.join(c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn') # Delete "accents"
        s = s.replace("l'", " ").replace("d'", " ") # Handle contractions like "L'Escala" -> "Escala"
        s = re.sub(r'\\b(?:el|la|els|les)\\b', ' ', s) # Remove common articles
        s = re.sub(r'[^a-z0-9]', ' ', s) # Replace non-alphanumeric with space
        return " ".join(s.split()) # Normalize whitespace

    # Data preprocessing
    df_loc = pd.read_csv("${coordinates_file}", sep="\\t")
    df_loc['Norm_Pob'] = df_loc['Población'].apply(normalize_str)

    # Create mapping dictionaries
    town_to_prov_raw = df_loc.set_index('Norm_Pob')['Provincia'].to_dict()
    town_to_prov = {k: town_to_prov_raw[k] for k in sorted(town_to_prov_raw.keys(), key=len, reverse=True) if k}

    town_mapping_raw = df_loc.set_index('Norm_Pob')['Población'].to_dict()
    town_mapping = {k: town_mapping_raw[k] for k in sorted(town_mapping_raw.keys(), key=len, reverse=True) if k}

    # Get distinct coordinate dictionaries for Provinces and Towns
    capitals = df_loc[df_loc['Población'] == df_loc['Provincia']]
    prov_coord_dict = capitals.set_index('Provincia')[['Latitud', 'Longitud']].to_dict('index')
    town_coord_dict = df_loc.set_index('Población')[['Latitud', 'Longitud']].to_dict('index')

    df_geno = pd.read_csv("${genotyping_file}")
    df_meta = pd.read_csv("${meta_str}", skipinitialspace=True) if "${meta_str}" else pd.DataFrame()

    if not df_meta.empty:
        df_meta.columns = [str(c).strip().upper() for c in df_meta.columns]
        if 'LOCATION' in df_meta.columns:
            def extract_location(loc):
                '''
                Extracts the town and province from a given location string.
                It first normalizes the string, then checks for specific town matches, and finally falls back to broad province matches if no town is found.
                '''
                if pd.isna(loc): return pd.Series(['Sense dades', 'Sense dades'])
                norm_spaced = f" {normalize_str(loc)} "

                # Match specific town. Assigns both Town and Province.
                for town_norm, town_real in town_mapping.items():
                    if f" {town_norm} " in norm_spaced: 
                        return pd.Series([town_real, town_to_prov[town_norm]])

                # If no town is found, check for broad province match
                if ' barcelona ' in norm_spaced: return pd.Series(['Sense dades', 'Barcelona'])
                if ' girona ' in norm_spaced or ' gerona ' in norm_spaced: return pd.Series(['Sense dades', 'Girona'])
                if ' tarragona ' in norm_spaced: return pd.Series(['Sense dades', 'Tarragona'])
                if ' lleida ' in norm_spaced or ' lerida ' in norm_spaced: return pd.Series(['Sense dades', 'Lleida'])

                return pd.Series(['Sense dades', 'Sense dades'])

            df_meta[['TOWN_GROUP', 'PROV_GROUP']] = df_meta['LOCATION'].apply(extract_location)

    # Normalize missing values and ensure string type for key columns
    for col in ['Clade', 'Genotype', 'Sub-genotype']:
        df_geno[col] = df_geno[col].fillna("Unassigned").astype(str).str.strip() if col in df_geno.columns else "Unassigned"

    # Logic to extract H subtype and define clade grouping based on protocol
    df_geno['H_Subtype'] = df_geno['Subtype'].astype(str).str.extract(r'(H[0-9]+)', expand=False).fillna('Unknown')
    
    if "${params.protocol}".upper() == "HUMAN":
        df_geno['H_Subtype'] = df_geno['H_Subtype'].replace({
            'H1': 'A(H1)pdm09',
            'H3': 'A(H3)'
        })
        df_geno['Root_Clade'] = df_geno['Clade'].apply(lambda c: ".".join(str(c).split('.')[:3]) + "-like" if str(c).count('.') > 2 else c)
    else:
        df_geno['Root_Clade'] = df_geno['Clade']

    df = pd.merge(df_geno, df_meta, left_on='SampleID', right_on='ID') if not df_meta.empty else df_geno.copy()
    
    # Establish columns for both resolution layers
    df['Poblacion'] = df['TOWN_GROUP'] if 'TOWN_GROUP' in df.columns else 'Sense dades'
    df['Provincia'] = df['PROV_GROUP'] if 'PROV_GROUP' in df.columns else 'Sense dades'

    # Define Season based on DATE column if it exists, otherwise assign a default
    if 'DATE' in df.columns:
        df['DATE'] = pd.to_datetime(df['DATE'], errors='coerce')
        iso = df['DATE'].dt.isocalendar()
        s_year = iso.year.where(iso.week >= 40, iso.year - 1) # Assign season based on ISO week (season starts in week 40)
        df['Season'] = s_year.astype(str) + "-" + (s_year + 1).astype(str) # Season format "2020-2021"
    else:
        df['Season'] = "Unknown Season"
    
    # Sort seasons in reverse order to automatically put the newest season at the top
    seasons = sorted([s for s in df['Season'].unique() if pd.notna(s) and s != "Unknown Season"], reverse=True)
    geo_levels = ['Province', 'Town']
    
    # Define valid views based on available data and protocol
    base_label = 'Subtype' if "${params.protocol}".upper() == "HUMAN" else 'Subtype (H)'
    valid_classifications = [{'label': base_label, 'id': 'subtypes', 'col': 'H_Subtype', 'filter': None}]
    
    # Add clade views for each H subtype
    for h in sorted([h for h in df['H_Subtype'].unique() if h != 'Unknown']):
        if not all(c == "-" for c in df[df['H_Subtype'] == h]['Clade'].unique()):
            valid_classifications.append({'label': f'Clades ({h})', 'id': f'clades_{h}', 'col': 'Root_Clade', 'filter': h})
            
    # Genotypes view is exclusive to AVIAN protocol
    if "${params.protocol}".upper() == "AVIAN":
        if df[(df['Clade'] == '2.3.4.4b') & (df['Genotype'] != '-')].shape[0] > 0:
            valid_classifications.append({'label': 'Genotypes (2.3.4.4b)', 'id': 'genotypes', 'col': 'Genotype', 'filter': '2.3.4.4b_clade'})

    # Consistent color mapping accross all seasons
    global_color_map = {}
    for classification in valid_classifications:
        view_map = {}
        if classification['id'] == 'genotypes':
            g_labels = sorted([g for g in df[df['Clade'] == '2.3.4.4b']['Genotype'].unique() if g not in ['-', 'Unassigned']])
            for i, g in enumerate(g_labels): view_map[g] = qualitative.Vivid[i % len(qualitative.Vivid)]
        elif classification['id'] == 'subtypes':
            s_labels = sorted([s for s in df['H_Subtype'].unique() if s != '-'])
            for s in s_labels:
                if s in ['Unassigned', 'Unknown', 'unassigned']:
                    view_map[s] = '#888888'
                else:
                    view_map[s] = subtype_base_colors.get(s, '#888888')
        else:
            base_hex = subtype_base_colors.get(classification['filter'], '#888888')
            c_labels = sorted([c for c in df[df['H_Subtype'] == classification['filter']]['Root_Clade'].unique() if c != '-'])
            # Generate shades only for clades that are not Unassigned or Unknown
            shades = generate_shades(base_hex, len([c for c in c_labels if c not in ['Unassigned', 'Unknown', 'unassigned']]))
            view_map = {}
            shade_iter = iter(shades)
            for c in c_labels:
                if c in ['Unassigned', 'Unknown', 'unassigned']:
                    view_map[c] = '#888888'
                else:
                    view_map[c] = next(shade_iter)
        global_color_map[classification['id']] = view_map

    # MAP CREATION
    m = folium.Map(location=[41.7, 1.8], zoom_start=8, tiles='Cartodb Positron')
    m.get_root().html.add_child(folium.Element("<h3 align='center' style='font-family: Arial; font-weight: bold; margin-top: 15px; color: #333;'>Influenza Geographic Distribution</h3>"))

    # Registry to keep track of layer names and colors
    trace_registry = {}
    all_view_colors = {}

    # Generate layers for each combination of season, geographic level, and classification
    for season in seasons:
        df_season = df[df['Season'] == season]
        
        for level in geo_levels:
            for classification in valid_classifications:
                # 3-Dimensional ID: Season | GeoLevel | Classification
                layer_id = f"{season}|{level}|{classification['id']}"
                fg = folium.FeatureGroup(name=layer_id, show=False)
                
                # Filter data based on current view
                df_view = df_season[df_season[classification['col']] != "-"].copy()
                if classification['filter'] == '2.3.4.4b_clade':
                    df_view = df_view[df_view['Clade'] == '2.3.4.4b']
                elif classification['filter']:
                    df_view = df_view[df_view['H_Subtype'] == classification['filter']]
                    
                if df_view.empty: continue

                color_map = global_color_map[classification['id']]
                current_labels = df_view[classification['col']].unique()
                all_view_colors[layer_id] = {label: color_map.get(label, '#999999') for label in current_labels if label != '-'}

                # Switch targeting logic based on geographic level
                coords_dict_to_use = prov_coord_dict if level == 'Province' else town_coord_dict
                loc_col = 'Provincia' if level == 'Province' else 'Poblacion'

                # Generate visual markers
                for loc_name, coords in coords_dict_to_use.items():
                    loc_df = df_view[df_view[loc_col] == loc_name]
                    if loc_df.empty: continue
                    
                    counts = loc_df[classification['col']].value_counts()
                    total = int(counts.sum())
                    
                    pie_colors = []
                    hover_details = []
                    current_pct = 0
                    
                    for label, count in counts.items():
                        pct = (count / total) * 100
                        color = color_map.get(label, '#999999')
                        
                        pie_colors.append(f"{color} {current_pct:.2f}% {current_pct + pct:.2f}%")
                        hover_line = f"<span style='color:{color}'>&#9608;</span> <b>{label}</b>: {int(count)} / {total} ({pct:.1f}%)"
                        
                        breakdown = []
                        
                        if classification['col'] == 'Root_Clade' and str(label).endswith("-like") and "${params.protocol}".upper() == "HUMAN":
                            sub_counts = loc_df[loc_df['Root_Clade'] == label]['Clade'].value_counts()
                            for sub_label, sub_count in sub_counts.items():
                                if str(sub_label) not in ["Unassigned", "-", "No dataset available", "nan"]:
                                    sub_pct = (sub_count / total) * 100
                                    breakdown.append(f"&nbsp;&nbsp;&nbsp;&nbsp;- <b>{sub_label}:</b> {int(sub_count)} / {total} ({sub_pct:.2f}%)")
                                    
                        elif classification['col'] == 'Genotype':
                            sub_counts = loc_df[loc_df['Genotype'] == label]['Sub-genotype'].value_counts()
                            for sub_label, sub_count in sub_counts.items():
                                if str(sub_label) not in ["Unassigned", "-", "None", "", "nan"]:
                                    sub_pct = (sub_count / total) * 100
                                    breakdown.append(f"&nbsp;&nbsp;&nbsp;&nbsp;- <b>{sub_label}:</b> {int(sub_count)} / {total} ({sub_pct:.2f}%)")
                        
                        if breakdown:
                            hover_line += "<br>" + "<br>".join(breakdown)
                            
                        hover_details.append(hover_line + "<br>")
                        current_pct += pct
                    
                    icon_size = int(55 + min(45, total * 2.5))
                    pie_html = f'''<div style="width:{icon_size}px; height:{icon_size}px; border-radius:50%; 
                                   background:conic-gradient({", ".join(pie_colors)}); 
                                   border:2px solid white; box-shadow:0 0 5px rgba(0,0,0,0.3);"></div>'''
                                   
                    hover_box_html = f"<div style='font-family:Arial; min-width:150px;'><b>{loc_name}</b><hr style='margin: 4px 0;'><b>Occurrences:</b> {total}<br><br>{''.join(hover_details)}</div>"                
                    folium.Marker(
                        location=[float(coords['Latitud']), float(coords['Longitud'])],
                        icon=folium.DivIcon(html=pie_html, icon_anchor=(icon_size/2, icon_size/2)),
                        tooltip=folium.Tooltip(hover_box_html)
                    ).add_to(fg)
                
                fg.add_to(m)
                trace_registry[layer_id] = fg.get_name()

    # HTML controls for filtering by season, geographic level, and classification
    # (Plotly's built-in controls are not compatible with multiple layers in folium)
    season_options = "".join([f'<option value="{s}">Season {s}</option>' for s in seasons])
    level_options = '<option value="Province">Province</option><option value="Town">City/Town</option>'
    class_options = "".join([f'<option value="{v["id"]}">{v["label"]}</option>' for v in valid_classifications])

    control_html = f'''
    <div style="position:fixed; top:20px; left:60px; z-index:9999; background:white; padding:15px; border-radius:8px; display:flex; flex-direction:column; gap:10px; box-shadow:0 4px 15px rgba(0,0,0,0.1); font-family:Arial; min-width:360px;">
        
        <div style="display:flex; gap:10px;">
            <div style="flex:1.2; min-width:120px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">SEASON</label><br>
                <select id="seasonSel" style="padding:5px; border-radius:4px; width:100%; box-sizing:border-box;">{season_options}</select>
            </div>
            <div style="flex:1; min-width:120px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">GEOGRAPHIC LEVEL</label><br>
                <select id="levelSel" style="padding:5px; border-radius:4px; width:100%; box-sizing:border-box;">{level_options}</select>
            </div>
            <div style="flex:1; min-width:120px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">CLASSIFICATION</label><br>
                <select id="classSel" style="padding:5px; border-radius:4px; width:100%; box-sizing:border-box;">{class_options}</select>
            </div>
        </div>
        
        <div id="legendContainer" style="border-top: 1px solid #eee; padding-top: 10px; max-height: 250px; overflow-y: auto;">
            <label style="font-size:10px; font-weight:bold; color:#666; text-transform: uppercase;">Legend</label>
            <div id="legendList" style="margin-top: 5px; display: flex; flex-direction: column; gap: 3px;"></div>
        </div>
    </div>

    <script>
        // REGISTRY: Maps our 3D logical keys (Season|Level|Classification) to Folium's internal map layer IDs
        const REGISTRY = {json.dumps(trace_registry)};
        // COLORS: Stores the specific hex codes for each classification label to build the legend
        const COLORS = {json.dumps(all_view_colors)};
        
        function update() {{
            // Generate the unique key based on the current dropdown selections
            const key = document.getElementById('seasonSel').value + "|" + document.getElementById('levelSel').value + "|" + document.getElementById('classSel').value;
            
            // Iterate through all map layers and show only the one matching our key, hiding the rest
            for (const k in REGISTRY) {{
                const layer = window[REGISTRY[k]];
                if (layer) k === key ? layer.addTo(window.map_instance) : window.map_instance.removeLayer(layer);
            }}
            
            // Clear the existing HTML legend
            const legendList = document.getElementById('legendList');
            legendList.innerHTML = "";
            
            // Rebuild the legend items if the current map view has assigned colors
            if (COLORS[key]) {{
                Object.keys(COLORS[key]).sort().forEach(label => {{
                    const color = COLORS[key][label];
                    const item = document.createElement('div');
                    item.style.display = 'flex'; item.style.alignItems = 'center'; item.style.fontSize = '12px';
                    item.innerHTML = `<span style="display:inline-block; width:12px; height:12px; background:\${{color}}; margin-right:8px; border-radius:2px; border:1px solid #ddd;"></span><span>\${{label}}</span>`;
                    legendList.appendChild(item);
                }});
            }}
        }}
        
        // Trigger the map update function whenever any dropdown value changes
        document.getElementById('seasonSel').addEventListener('change', update);
        document.getElementById('levelSel').addEventListener('change', update);
        document.getElementById('classSel').addEventListener('change', update);
        
        // Wait for Folium to initialize the map, store its reference, and apply the initial filter
        window.addEventListener('load', () => {{
            for (let o in window) if (o.startsWith('map_')) {{ window.map_instance = window[o]; update(); break; }}
        }});
    </script>
    '''
    m.get_root().html.add_child(folium.Element(control_html))
    m.save("GeographicReport.html")
    """
}
