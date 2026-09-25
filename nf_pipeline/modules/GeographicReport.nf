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
        s = s.replace("l'", " ").replace("d'", " ").replace("l\\u2019", " ").replace("d\\u2019", " ") # Handle contractions like "L'Escala" -> "Escala" (straight or curly apostrophe)
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

        # Remove overlapping columns from genotyping data to avoid conflicts
        geno_overlap = [
            c for c in df_geno.columns
            if c.strip().upper() in df_meta.columns.tolist() and c.strip().upper() != 'SAMPLEID'
        ]
        if geno_overlap:
            df_geno = df_geno.drop(columns=geno_overlap)

        if 'LOCATION' in df_meta.columns:
            def extract_location(loc):
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
    # Normalize missing values and ensure string type for key columns
    for col in ['Clade', 'Genotype', 'Sub-genotype']:
        df_geno[col] = df_geno[col].fillna("Unassigned").astype(str).str.strip() if col in df_geno.columns else "Unassigned"
        
    df_geno['Clade'] = df_geno['Clade'].replace('-', 'Unassigned')
    df_geno['Genotype'] = df_geno['Genotype'].replace('-', 'Unassigned')
    # Logic to extract H subtype and define clade grouping based on protocol
    df_geno['H_Subtype'] = df_geno['Subtype'].astype(str).str.extract(r'(H[0-9]+)', expand=False).fillna('Unknown')
    
    if "${params.protocol}".upper() == "HUMAN":
        df_geno['H_Subtype'] = df_geno['H_Subtype'].replace({
            'H1': 'A(H1)pdm09',
            'H3': 'A(H3)'
        })
        all_clades = df_geno['Clade'].dropna().unique()
        def get_root_clade(c):
            if pd.isna(c) or c in ["Unassigned", "No dataset available"]: return c
            parts = str(c).split('.')
            if len(parts) > 3: return ".".join(parts[:3]) + "-like"
            elif len(parts) == 3 and any(str(x).startswith(str(c) + ".") for x in all_clades): return str(c) + "-like"
            return c
        df_geno['Root_Clade'] = df_geno['Clade'].apply(get_root_clade)
    else:
        df_geno['Root_Clade'] = df_geno['Clade']

    # One metadata row per sample (as in the other reports) so duplicated IDs are not counted twice on the map
    if not df_meta.empty and 'ID' in df_meta.columns:
        df_meta = df_meta.drop_duplicates(subset=['ID'], keep='first')

    df = pd.merge(df_geno, df_meta, left_on='SampleID', right_on='ID') if not df_meta.empty else df_geno.copy()
    
    # Establish columns for resolution layers and demographics
    df['Poblacion'] = df['TOWN_GROUP'] if 'TOWN_GROUP' in df.columns else 'Sense dades'
    df['Provincia'] = df['PROV_GROUP'] if 'PROV_GROUP' in df.columns else 'Sense dades'
    df['Originating Lab'] = df['ORIGINATING_LAB'].fillna('Sense dades').astype(str).str.strip() if 'ORIGINATING_LAB' in df.columns else 'Sense dades'
    # Normalize laboratory names using the town_mapping dictionary
    df['Originating Lab'] = df['Originating Lab'].apply(
        lambda v: town_mapping.get(normalize_str(v), v)
    )
    
    if 'AGE GROUP' in df.columns:
        df['Age_Group'] = df['AGE GROUP'].fillna('Sense dades').astype(str).str.strip()
    else:
        df['Age_Group'] = 'Sense dades'
        
    if 'SEX' in df.columns:
        df['Sex'] = df['SEX'].fillna('Sense dades').astype(str).str.strip()
    else:
        df['Sex'] = 'Sense dades'

    # Define Season based on DATE column if it exists, otherwise assign a default
    if 'DATE' in df.columns:
        df['DATE'] = pd.to_datetime(df['DATE'], errors='coerce')
        iso = df['DATE'].dt.isocalendar()
        s_year = iso.year.where(iso.week >= 40, iso.year - 1) # Assign season based on ISO week (season starts in week 40)
        df['Season'] = s_year.astype(str) + "-" + (s_year + 1).astype(str) # Season format "2020-2021"
        df.loc[df['DATE'].isna(), 'Season'] = "Unknown Season" # Undated samples only count in "All Time"
    else:
        df['Season'] = "Unknown Season"
    
    # Setup iterative parameters dynamically from the dataset
    seasons = sorted([s for s in df['Season'].unique() if pd.notna(s) and s != "Unknown Season"], reverse=True)
    geo_levels = ['Province', 'Town', 'Originating Lab']
    age_order = {'0-2': 0, '3-4': 1, '5-14': 2, '15-65': 3, '>65': 4}
    age_groups = ['All'] + sorted(
        [a for a in df['Age_Group'].unique() if str(a) != 'nan'],
        key=lambda x: age_order.get(str(x).strip(), 99)
    )
    sexs = ['All'] + sorted([g for g in df['Sex'].unique() if str(g) != 'nan'])
    
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
            s_labels = sorted([s for s in df['H_Subtype'].unique() if s not in ['-', 'Unknown']])
            for s in s_labels: view_map[s] = subtype_base_colors.get(s, '#888888')
        else:
            base_hex = subtype_base_colors.get(classification['filter'], '#888888')
            c_labels = sorted([c for c in df[df['H_Subtype'] == classification['filter']]['Root_Clade'].unique() if c not in ['-', 'Unassigned']])
            shades = generate_shades(base_hex, len(c_labels) + 2)
            for clade, shade in zip(c_labels, shades): view_map[clade] = shade
            
        # Manually force the Unassigned category to light grey
        view_map['Unassigned'] = '#d3d3d3'
        global_color_map[classification['id']] = view_map

    # MAP CREATION
    m = folium.Map(
        location=[41.7, 1.8],
        zoom_start=8,
        tiles='https://basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png?key=${params.carto_api_key}',
        attr='&copy; <a href="https://carto.com/attributions">CARTO</a>'
    )
    m.get_root().html.add_child(folium.Element("<h3 align='center' style='font-family: Arial; font-weight: bold; margin-top: 15px; color: #333;'>Influenza Geographic Distribution</h3>"))

    # Registry to keep track of layer names and colors
    trace_registry = {}
    all_view_colors = {}
    # Per-season marker counts ("level|class|age|sex" -> season -> markers) so the page can sum a season range
    range_data = {}

    for season in ['All Time'] + seasons:
        df_season = df if season == 'All Time' else df[df['Season'] == season]

        for level in geo_levels:
            for classification in valid_classifications:
                for age in age_groups:
                    for sex in sexs:
                        # 5-Dimensional ID ensures uniqueness for dynamic filtering
                        layer_id = f"{season}|{level}|{classification['id']}|{age}|{sex}"
                        
                        df_view = df_season.copy()
                        if age != 'All': df_view = df_view[df_view['Age_Group'] == age]
                        if sex != 'All': df_view = df_view[df_view['Sex'] == sex]
                        
                        # Filter data based on current view
                        if classification['filter'] == '2.3.4.4b_clade':
                            df_view = df_view[df_view['Clade'] == '2.3.4.4b']
                        elif classification['filter']:
                            df_view = df_view[df_view['H_Subtype'] == classification['filter']]
                            
                        if df_view.empty: continue
                        
                        fg = folium.FeatureGroup(name=layer_id, show=False)

                        color_map = global_color_map[classification['id']]
                        current_labels = df_view[classification['col']].unique()
                        all_view_colors[layer_id] = {label: color_map.get(label, '#999999') for label in current_labels if label != '-'}

                        # Switch targeting logic based on geographic level
                        if level == 'Province':
                            coords_dict_to_use = prov_coord_dict
                            loc_col = 'Provincia'
                        elif level == 'Town':
                            coords_dict_to_use = town_coord_dict
                            loc_col = 'Poblacion'
                        else:
                            coords_dict_to_use = town_coord_dict
                            loc_col = 'Originating Lab'

                        # Generate visual markers
                        range_markers = []
                        for loc_name, coords in coords_dict_to_use.items():
                            loc_df = df_view[df_view[loc_col] == loc_name]
                            if loc_df.empty: continue
                            
                            counts = loc_df[classification['col']].value_counts()
                            total = int(counts.sum())
                            
                            pie_colors = []
                            hover_details = []
                            current_pct = 0
                            
                            label_rows = []
                            for label, count in counts.items():
                                pct = (count / total) * 100
                                color = color_map.get(label, '#999999')
                                
                                pie_colors.append(f"{color} {current_pct:.2f}% {current_pct + pct:.2f}%")
                                hover_line = f"<span style='color:{color}'>&#9608;</span> <b>{label}</b>: {int(count)} / {total} ({pct:.1f}%)"
                                
                                breakdown = []
                                sub_rows = []
                                
                                if classification['col'] == 'Root_Clade' and str(label).endswith("-like") and "${params.protocol}".upper() == "HUMAN":
                                    sub_counts = loc_df[loc_df['Root_Clade'] == label]['Clade'].value_counts()
                                    for sub_label, sub_count in sub_counts.items():
                                        if str(sub_label) not in ["Unassigned", "-", "No dataset available", "nan"]:
                                            sub_pct = (sub_count / total) * 100
                                            sub_rows.append([sub_label, int(sub_count)])
                                            breakdown.append(f"&nbsp;&nbsp;&nbsp;&nbsp;- <b>{sub_label}:</b> {int(sub_count)} / {total} ({sub_pct:.2f}%)")
                                            
                                elif classification['col'] == 'Genotype':
                                    sub_counts = loc_df[loc_df['Genotype'] == label]['Sub-genotype'].value_counts()
                                    for sub_label, sub_count in sub_counts.items():
                                        if str(sub_label) not in ["Unassigned", "-", "None", "", "nan"]:
                                            sub_pct = (sub_count / total) * 100
                                            sub_rows.append([sub_label, int(sub_count)])
                                            breakdown.append(f"&nbsp;&nbsp;&nbsp;&nbsp;- <b>{sub_label}:</b> {int(sub_count)} / {total} ({sub_pct:.2f}%)")
                                
                                if breakdown:
                                    hover_line += "<br>" + "<br>".join(breakdown)
                                    
                                hover_details.append(hover_line + "<br>")
                                current_pct += pct
                                label_rows.append([label, int(count), sub_rows])
                            range_markers.append({'n': loc_name, 'c': [float(coords['Latitud']), float(coords['Longitud'])], 'l': label_rows})
                            
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
                        if season != 'All Time':
                            range_data.setdefault(f"{level}|{classification['id']}|{age}|{sex}", {})[season] = range_markers

    # HTML controls
    # Season range: FROM = "All Time" or a season, TO = a season (disabled while FROM is "All Time").
    # Initial view is the latest season (FROM = TO = latest), or "All Time" when no season is dated.
    default_season = seasons[0] if seasons else 'All Time'
    season_from_options = "".join([f'<option value="{s}"{" selected" if s == default_season else ""}>{s if s == "All Time" else "Season " + s}</option>' for s in ['All Time'] + seasons])
    season_to_options = "".join([f'<option value="{s}">Season {s}</option>' for s in seasons])
    season_to_disabled = " disabled" if default_season == 'All Time' else ""
    level_options = '<option value="Province">Province</option><option value="Town">City/Town</option><option value="Originating Lab">Originating Lab</option>'
    class_options = "".join([f'<option value="{v["id"]}">{v["label"]}</option>' for v in valid_classifications])
    age_options = "".join([f'<option value="{a}">{a}</option>' for a in age_groups])
    sex_options = "".join([f'<option value="{g}">{g}</option>' for g in sexs])

    control_html = f'''
    <div style="position:fixed; top:20px; left:60px; z-index:9999; background:white; padding:15px; border-radius:8px; display:flex; flex-direction:column; gap:10px; box-shadow:0 4px 15px rgba(0,0,0,0.1); font-family:Arial; min-width:650px;">
        <div style="display:flex; gap:10px;">
            <div style="flex:1.5; min-width:130px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">SEASON FROM</label><br>
                <select id="seasonSel" style="padding:5px; border-radius:4px; width:100%; box-sizing:border-box;">{season_from_options}</select>
            </div>
            <div style="flex:1.5; min-width:130px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">SEASON TO</label><br>
                <select id="seasonToSel" style="padding:5px; border-radius:4px; width:100%; box-sizing:border-box;"{season_to_disabled}>{season_to_options}</select>
            </div>
            <div style="flex:1; min-width:100px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">GEOGRAPHIC LEVEL</label><br>
                <select id="levelSel" style="padding:5px; border-radius:4px; width:100%; box-sizing:border-box;">{level_options}</select>
            </div>
            <div style="flex:1; min-width:100px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">CLASSIFICATION</label><br>
                <select id="classSel" style="padding:5px; border-radius:4px; width:100%; box-sizing:border-box;">{class_options}</select>
            </div>
            <div style="flex:1; min-width:90px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">AGE GROUP</label><br>
                <select id="ageSel" style="padding:5px; border-radius:4px; width:100%; box-sizing:border-box;">{age_options}</select>
            </div>
            <div style="flex:1; min-width:90px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">SEX</label><br>
                <select id="sexSel" style="padding:5px; border-radius:4px; width:100%; box-sizing:border-box;">{sex_options}</select>
            </div>
        </div>
        <div id="legendContainer" style="border-top: 1px solid #eee; padding-top: 10px; max-height: 250px; overflow-y: auto;">
            <label style="font-size:10px; font-weight:bold; color:#666; text-transform: uppercase;">Legend</label>
            <div id="legendList" style="margin-top: 5px; display: flex; flex-direction: column; gap: 3px;"></div>
        </div>
    </div>
    <script>
        const REGISTRY = {json.dumps(trace_registry)};
        const COLORS = {json.dumps(all_view_colors)};
        const SEASONS = {json.dumps(sorted(seasons))};
        const RANGE_DATA = {json.dumps(range_data)};
        const COLOR_MAPS = {json.dumps(global_color_map)};
        const RANGE_LAYERS = {{}};
        let rangeLayer = null;

        // Number formatting as in Python (exact halves round to even: 6.25 -> "6.2"), unlike toFixed
        function pyFixed(x, d) {{
            const full = x.toFixed(d + 30);
            const cut = full.indexOf('.') + d + 1;
            const kept = full.slice(0, cut);
            if (/^50*\$/.test(full.slice(cut)) && /[02468]\$/.test(kept)) return kept;
            return x.toFixed(d);
        }}

        // Sum the per-season counts of every season from..to (inclusive) into one marker per location
        function buildRangeLayer(from, to, rest) {{
            const perSeason = RANGE_DATA[rest] || {{}};
            const colorMap = COLOR_MAPS[rest.split('|')[1]] || {{}};
            const locs = {{}}, locOrder = [], legend = {{}};
            SEASONS.filter(s => s >= from && s <= to).forEach(s => {{
                (perSeason[s] || []).forEach(mk => {{
                    let loc = locs[mk.n];
                    if (!loc) {{ loc = locs[mk.n] = {{ c: mk.c, labels: {{}}, order: [] }}; locOrder.push(mk.n); }}
                    mk.l.forEach(row => {{
                        let e = loc.labels[row[0]];
                        if (!e) {{ e = loc.labels[row[0]] = {{ n: 0, subs: {{}}, subOrder: [] }}; loc.order.push(row[0]); }}
                        e.n += row[1];
                        row[2].forEach(sub => {{
                            if (!(sub[0] in e.subs)) {{ e.subs[sub[0]] = 0; e.subOrder.push(sub[0]); }}
                            e.subs[sub[0]] += sub[1];
                        }});
                    }});
                }});
            }});
            const group = L.featureGroup();
            locOrder.forEach(name => {{
                const loc = locs[name];
                // Same order as pandas value_counts: highest count first
                const labels = loc.order.slice().sort((a, b) => loc.labels[b].n - loc.labels[a].n);
                const total = labels.reduce((acc, lab) => acc + loc.labels[lab].n, 0);
                const pieColors = [], hover = [];
                let currentPct = 0;
                labels.forEach(lab => {{
                    const e = loc.labels[lab];
                    const pct = (e.n / total) * 100;
                    const color = colorMap[lab] || '#999999';
                    if (lab !== '-') legend[lab] = color;
                    pieColors.push(color + " " + pyFixed(currentPct, 2) + "% " + pyFixed(currentPct + pct, 2) + "%");
                    let line = "<span style='color:" + color + "'>&#9608;</span> <b>" + lab + "</b>: " + e.n + " / " + total + " (" + pyFixed(pct, 1) + "%)";
                    const subs = e.subOrder.slice().sort((a, b) => e.subs[b] - e.subs[a]).map(sl =>
                        "&nbsp;&nbsp;&nbsp;&nbsp;- <b>" + sl + ":</b> " + e.subs[sl] + " / " + total + " (" + pyFixed((e.subs[sl] / total) * 100, 2) + "%)");
                    if (subs.length) line += "<br>" + subs.join("<br>");
                    hover.push(line + "<br>");
                    currentPct += pct;
                }});
                const iconSize = Math.trunc(55 + Math.min(45, total * 2.5));
                const pieHtml = '<div style="width:' + iconSize + 'px; height:' + iconSize + 'px; border-radius:50%; background:conic-gradient(' + pieColors.join(", ") + '); border:2px solid white; box-shadow:0 0 5px rgba(0,0,0,0.3);"></div>';
                const hoverHtml = "<div style='font-family:Arial; min-width:150px;'><b>" + name + "</b><hr style='margin: 4px 0;'><b>Occurrences:</b> " + total + "<br><br>" + hover.join("") + "</div>";
                L.marker(loc.c, {{ icon: L.divIcon({{ html: pieHtml, iconAnchor: [iconSize / 2, iconSize / 2], className: 'empty' }}) }})
                    .bindTooltip("<div>" + hoverHtml + "</div>", {{ sticky: true }})
                    .addTo(group);
            }});
            return {{ layer: group, colors: legend }};
        }}

        function onSeasonFromChange() {{
            const fromSel = document.getElementById('seasonSel');
            const toSel = document.getElementById('seasonToSel');
            toSel.disabled = fromSel.value === 'All Time';
            if (!toSel.disabled && toSel.value < fromSel.value) toSel.value = fromSel.value;
            update();
        }}

        function onSeasonToChange() {{
            const fromSel = document.getElementById('seasonSel');
            const toSel = document.getElementById('seasonToSel');
            if (toSel.value < fromSel.value) toSel.value = fromSel.value;
            update();
        }}
        
        function update() {{
            const from = document.getElementById('seasonSel').value;
            const to = document.getElementById('seasonToSel').value;
            const rest = document.getElementById('levelSel').value + "|" + 
                        document.getElementById('classSel').value + "|" +
                        document.getElementById('ageSel').value + "|" +
                        document.getElementById('sexSel').value;
            // "All Time" or a single season use the precomputed layers; a longer range is summed on the fly
            const isRange = from !== 'All Time' && to && to !== from;
            const key = (from === 'All Time' || !to ? from : to === from ? from : from + ".." + to) + "|" + rest;
                        
            for (const k in REGISTRY) {{
                const layer = window[REGISTRY[k]];
                if (layer) k === key ? layer.addTo(window.map_instance) : window.map_instance.removeLayer(layer);
            }}
            if (rangeLayer) {{ window.map_instance.removeLayer(rangeLayer); rangeLayer = null; }}
            let colors = COLORS[key];
            if (isRange) {{
                if (!RANGE_LAYERS[key]) RANGE_LAYERS[key] = buildRangeLayer(from, to, rest);
                rangeLayer = RANGE_LAYERS[key].layer;
                rangeLayer.addTo(window.map_instance);
                colors = RANGE_LAYERS[key].colors;
            }}
            
            const legendList = document.getElementById('legendList');
            legendList.innerHTML = "";
            if (colors) {{
                Object.keys(colors).sort().forEach(label => {{
                    const color = colors[label];
                    const item = document.createElement('div');
                    item.style.display = 'flex'; item.style.alignItems = 'center'; item.style.fontSize = '12px';
                    item.innerHTML = `<span style="display:inline-block; width:12px; height:12px; background:\${{color}}; margin-right:8px; border-radius:2px; border:1px solid #ddd;"></span><span>\${{label}}</span>`;
                    legendList.appendChild(item);
                }});
            }}
        }}
        
        document.getElementById('seasonSel').addEventListener('change', onSeasonFromChange);
        document.getElementById('seasonToSel').addEventListener('change', onSeasonToChange);
        document.getElementById('levelSel').addEventListener('change', update);
        document.getElementById('classSel').addEventListener('change', update);
        document.getElementById('ageSel').addEventListener('change', update);
        document.getElementById('sexSel').addEventListener('change', update);
        
        window.addEventListener('load', () => {{
            for (let o in window) if (o.startsWith('map_')) {{ window.map_instance = window[o]; break; }}
            update();
        }});
    </script>
    '''
    m.get_root().html.add_child(folium.Element(control_html))
    m.save("GeographicReport.html")
    """
}
