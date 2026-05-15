process GeographicReport {
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

    # Color palette for subtypes
    subtype_base_colors = {
        'H1': '#1F77B4', 'H1N1': '#1F77B4','H2': '#D62728','H3': '#FF7F0E', 'H3N2': '#FF7F0E', 'H4': '#2CA02C','H5': '#9467BD','H6': '#8C564B', 
        'H7': '#E377C2','H8': '#7F7F7F','H9': '#BCBD22','H10': '#17BECF','H11': '#393B79','H12': '#637939',
        'H13': '#8C6D31','H14': '#843C39','H15': '#7B4173','H16': '#5254A3','H17': '#8CA252','H18': '#BD9E39' 
    }

    def generate_shades(base_hex, n):
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

    # Create a mapping from normalized town names to provinces, sorted by length to prioritize longer matches
    town_to_prov_raw = df_loc.set_index('Norm_Pob')['Provincia'].to_dict()
    town_to_prov = {k: town_to_prov_raw[k] for k in sorted(town_to_prov_raw.keys(), key=len, reverse=True) if k}

    # Get coordinates of provincial capitals for marker placement
    capitals = df_loc[df_loc['Población'] == df_loc['Provincia']]
    coord_dict = capitals.set_index('Provincia')[['Latitud', 'Longitud']].to_dict('index')

    df_geno = pd.read_csv("${genotyping_file}")
    df_meta = pd.read_csv("${meta_str}", skipinitialspace=True) if "${meta_str}" else pd.DataFrame()

    if not df_meta.empty:
        df_meta.columns = [str(c).strip().upper() for c in df_meta.columns]
        if 'LOCATION' in df_meta.columns:
            def extract_group(loc):
                if pd.isna(loc): return 'Sense dades'
                norm_spaced = f" {normalize_str(loc)} "
                # If location contains a province name, assign to that province group
                if ' barcelona ' in norm_spaced: return 'Barcelona'
                if ' girona ' in norm_spaced or ' gerona ' in norm_spaced: return 'Girona'
                if ' tarragona ' in norm_spaced: return 'Tarragona'
                if ' lleida ' in norm_spaced or ' lerida ' in norm_spaced: return 'Lleida'

                # If location contains a town name, assign to the corresponding province group
                for town_norm, prov in town_to_prov.items():
                    if f" {town_norm} " in norm_spaced: 
                        return prov

                return 'Sense dades'

            df_meta['PROV_GROUP'] = df_meta['LOCATION'].apply(extract_group)


    # Normalize missing values and ensure string type for key columns
    for col in ['Clade', 'Genotype', 'Sub-genotype']:
        df_geno[col] = df_geno[col].fillna("Unassigned").astype(str).str.strip() if col in df_geno.columns else "Unassigned"

    # Logic to extract H subtype and define clade grouping based on protocol
    if "${params.protocol}".upper() == "HUMAN":
        df_geno['H_Subtype'] = df_geno['Subtype'].astype(str).str.extract(r'(H[0-9]+N[0-9]+)', expand=False).fillna('Unknown')
        
        # Filtre estricte: Si no és H1N1 o H3N2, ho marquem com a error/desconegut
        df_geno.loc[~df_geno['H_Subtype'].isin(['H1N1', 'H3N2']), 'H_Subtype'] = 'Unknown'
        
        df_geno['Root_Clade'] = df_geno['Clade'].apply(lambda c: ".".join(c.split('.')[:3]) + "-like" if c.count('.') > 2 else c)
    else:
        df_geno['H_Subtype'] = df_geno['Subtype'].astype(str).str.extract(r'(H[0-9]+)', expand=False).fillna('Unknown')
        df_geno['Root_Clade'] = df_geno['Clade']

    df = pd.merge(df_geno, df_meta, left_on='SampleID', right_on='ID') if not df_meta.empty else df_geno.copy()
    df['Provincia'] = df['PROV_GROUP'] if 'PROV_GROUP' in df.columns else 'Sense dades'

    # Define Season based on DATE column if it exists, otherwise assign "All Time"
    if 'DATE' in df.columns:
        df['DATE'] = pd.to_datetime(df['DATE'], errors='coerce')
        iso = df['DATE'].dt.isocalendar()
        s_year = iso.year.where(iso.week >= 40, iso.year - 1) # Assign season based on ISO week (season starts in week 40)
        df['Season'] = s_year.astype(str) + "-" + (s_year + 1).astype(str) # Season format "2020-2021"
    else:
        df['Season'] = "All Time"
    
    seasons = ["All Time"] + sorted([s for s in df['Season'].unique() if pd.notna(s) and s != "Unknown Season"])
    
    # Define valid views based on available data and protocol
    base_label = 'Subtype' if "${params.protocol}".upper() == "HUMAN" else 'Subtype (H)'
    valid_views = [{'label': base_label, 'id': 'subtypes', 'col': 'H_Subtype', 'filter': None}]
    
    # Add clade views for each H subtype
    for h in sorted([h for h in df['H_Subtype'].unique() if h != 'Unknown']):
        if not all(c == "-" for c in df[df['H_Subtype'] == h]['Clade'].unique()):
            valid_views.append({'label': f'Clades ({h})', 'id': f'clades_{h}', 'col': 'Root_Clade', 'filter': h})
            
    # Genotypes view is exclusive to AVIAN protocol
    if "${params.protocol}".upper() == "AVIAN":
        if df[(df['Clade'] == '2.3.4.4b') & (df['Genotype'] != '-')].shape[0] > 0:
            valid_views.append({'label': 'Genotypes (2.3.4.4b)', 'id': 'genotypes', 'col': 'Genotype', 'filter': '2.3.4.4b_clade'})

    # MAP CREATION
    m = folium.Map(location=[41.7, 1.8], zoom_start=8, tiles='Cartodb Positron')
    m.get_root().html.add_child(folium.Element("<h3 align='center' style='font-family: Arial; font-weight: bold; margin-top: 15px; color: #333;'>Influenza Provincial Distribution</h3>"))

    # Registry to keep track of layer names and colors
    trace_registry = {}
    all_view_colors = {}

    for season in seasons:
        df_season = df if season == "All Time" else df[df['Season'] == season]
        
        for view in valid_views:
            season_view_id = f"{season}|{view['id']}"
            fg = folium.FeatureGroup(name=season_view_id, show=False)
            
            # Filter data based on current view
            df_view = df_season[df_season[view['col']] != "-"].copy()
            if view['filter'] == '2.3.4.4b_clade':
                df_view = df_view[df_view['Clade'] == '2.3.4.4b']
            elif view['filter']:
                df_view = df_view[df_view['H_Subtype'] == view['filter']]
                
            if df_view.empty: continue

            # Build a global color map for legend consistency
            valid_labels = sorted([l for l in df_view[view['col']].unique() if l not in ['Unassigned', 'Unknown']])
            color_map = {}
            # Pre-generate shades if this is a clade view
            if view['id'] not in ['genotypes', 'subtypes'] and len(valid_labels) > 0:
                base_hex = subtype_base_colors.get(view['filter'], '#888888')
                clade_shades = generate_shades(base_hex, len(valid_labels))
            
            for i, label in enumerate(valid_labels):
                if view['id'] == 'genotypes': 
                    color_map[label] = qualitative.Vivid[i % len(qualitative.Vivid)]
                elif view['id'] == 'subtypes': 
                    color_map[label] = subtype_base_colors.get(label, '#888888')
                else: 
                    # Apply the true hex shades generated by the function
                    color_map[label] = clade_shades[i]
            
            all_view_colors[season_view_id] = color_map

            # Generate visual markers per province
            for prov_name, coords in coord_dict.items():
                prov_df = df_view[df_view['Provincia'] == prov_name]
                if prov_df.empty: continue
                
                counts = prov_df[view['col']].value_counts()
                total = int(counts.sum())
                
                # Build the pie chart colors and the text that appears when you hover over the marker
                pie_colors = []
                hover_details = []
                current_pct = 0
                
                for label, count in counts.items():
                    pct = (count / total) * 100
                    color = color_map.get(label, '#999999')
                    
                    # Add color slice to the pie chart
                    pie_colors.append(f"{color} {current_pct:.2f}% {current_pct + pct:.2f}%")
                    # Base string for the hover box
                    hover_line = f"<span style='color:{color}'>&#9608;</span> <b>{label}</b>: {int(count)} / {total} ({pct:.1f}%)"
                    
                    # Dynamic Breakdown Logic
                    breakdown = []
                    
                    # Check if we are in a clade view and the label ends with '-like'
                    if view['col'] == 'Root_Clade' and str(label).endswith("-like") and "${params.protocol}".upper() == "HUMAN":
                        sub_counts = prov_df[prov_df['Root_Clade'] == label]['Clade'].value_counts()
                        for sub_label, sub_count in sub_counts.items():
                            if str(sub_label) not in ["Unassigned", "-", "No dataset available", "nan"]:
                                sub_pct = (sub_count / total) * 100
                                breakdown.append(f"&nbsp;&nbsp;&nbsp;&nbsp;- <b>{sub_label}:</b> {int(sub_count)} / {total} ({sub_pct:.2f}%)")
                                
                    # Check if we are in the genotype view to break down sub-genotypes
                    elif view['col'] == 'Genotype':
                        sub_counts = prov_df[prov_df['Genotype'] == label]['Sub-genotype'].value_counts()
                        for sub_label, sub_count in sub_counts.items():
                            if str(sub_label) not in ["Unassigned", "-", "None", "", "nan"]:
                                sub_pct = (sub_count / total) * 100
                                breakdown.append(f"&nbsp;&nbsp;&nbsp;&nbsp;- <b>{sub_label}:</b> {int(sub_count)} / {total} ({sub_pct:.2f}%)")
                    
                    # If we found valid sub-details, append them to the hover line
                    if breakdown:
                        hover_line += "<br>" + "<br>".join(breakdown)
                        
                    # Add the finished line (and a trailing break) to our hover box list
                    hover_details.append(hover_line + "<br>")
                    current_pct += pct
                
                icon_size = int(55 + min(45, total * 2.5))
                
                # Create a pie chart using CSS conic-gradient for the marker icon
                pie_html = f'''<div style="width:{icon_size}px; height:{icon_size}px; border-radius:50%; 
                               background:conic-gradient({", ".join(pie_colors)}); 
                               border:2px solid white; box-shadow:0 0 5px rgba(0,0,0,0.3);"></div>'''
                               
                # The hover
                hover_box_html = f"<div style='font-family:Arial; min-width:150px;'><b>{prov_name}</b><hr style='margin: 4px 0;'><b>Occurrences:</b> {total}<br><br>{''.join(hover_details)}</div>"                
                folium.Marker(
                    location=[float(coords['Latitud']), float(coords['Longitud'])],
                    icon=folium.DivIcon(html=pie_html, icon_anchor=(icon_size/2, icon_size/2)),
                    tooltip=folium.Tooltip(hover_box_html)
                ).add_to(fg)
            
            fg.add_to(m)
            trace_registry[season_view_id] = fg.get_name()

    # HTML controls for season and view selection and legend management
    season_options = "".join([f'<option value="{season}">{"All Seasons" if season=="All Time" else "Season "+season}</option>' for season in seasons])
    view_options = "".join([f'<option value="{view["id"]}">{view["label"]}</option>' for view in valid_views])

    control_html = f'''
    <div style="position:fixed; top:20px; left:60px; z-index:9999; background:white; padding:15px; border-radius:8px; display:flex; flex-direction:column; gap:10px; box-shadow:0 4px 15px rgba(0,0,0,0.1); font-family:Arial; min-width:220px;">
        <div style="display:flex; gap:20px;">
            <div><label style="font-size:10px; font-weight:bold; color:#666;">SEASON</label><br><select id="seasonSel" style="padding:5px; border-radius:4px; width:100%;">{season_options}</select></div>
            <div><label style="font-size:10px; font-weight:bold; color:#666;">VIEW</label><br><select id="viewSel" style="padding:5px; border-radius:4px; width:100%;">{view_options}</select></div>
        </div>
        <div id="legendContainer" style="border-top: 1px solid #eee; padding-top: 10px; max-height: 250px; overflow-y: auto;">
            <label style="font-size:10px; font-weight:bold; color:#666; text-transform: uppercase;">Legend</label>
            <div id="legendList" style="margin-top: 5px; display: flex; flex-direction: column; gap: 3px;"></div>
        </div>
    </div>
    <script>
        // Registry of layer names and colors to manage visibility and legend based on user selection
        const REGISTRY = {json.dumps(trace_registry)};
        const COLORS = {json.dumps(all_view_colors)};
        // Function to update map layers based on selected season and view
        function update() {{
            const key = document.getElementById('seasonSel').value + "|" + document.getElementById('viewSel').value;
            for (const k in REGISTRY) {{
                const layer = window[REGISTRY[k]];
                if (layer) k === key ? layer.addTo(window.map_instance) : window.map_instance.removeLayer(layer);
            }}
            const legendList = document.getElementById('legendList');
            legendList.innerHTML = "";
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
        // Event listeners for dropdown changes
        document.getElementById('seasonSel').addEventListener('change', update);
        document.getElementById('viewSel').addEventListener('change', update);
        // Initialize map with default selections
        window.addEventListener('load', () => {{
            for (let o in window) if (o.startsWith('map_')) {{ window.map_instance = window[o]; update(); break; }}
        }});
    </script>
    '''
    m.get_root().html.add_child(folium.Element(control_html))
    m.save("GeographicReport.html")
    """
}