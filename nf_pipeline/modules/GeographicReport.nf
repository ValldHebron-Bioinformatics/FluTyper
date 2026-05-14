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
    import os
    import re
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from plotly import colors as plotly_colors
    import colorsys
    import random
    import unicodedata

    # 1. HELPERS DE COLORS I NORMALITZACIÓ EXTREMA
    is_colorblind = "${params.colorblind}".lower() == "true"
    okabe_ito_colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999', '#000000']
    
    subtype_base_colors = {
        'H1': '#1F77B4','H2': '#D62728','H3': '#FF7F0E','H4': '#2CA02C','H5': '#9467BD','H6': '#8C564B', 
        'H7': '#E377C2','H8': '#7F7F7F','H9': '#BCBD22','H10': '#17BECF','H11': '#393B79','H12': '#637939',
        'H13': '#8C6D31','H14': '#843C39','H15': '#7B4173','H16': '#5254A3','H17': '#8CA252','H18': '#BD9E39' 
    }

    def generate_shades(base_hex, n):
        if n <= 0: return []
        if n == 1: return [base_hex]
        clean_hex = str(base_hex).lstrip('#')
        try:
            r, g, b = [int(clean_hex[i:i+2], 16) / 255.0 for i in (0, 2, 4)]
        except: return ['#888888'] * n
        h, l, s = colorsys.rgb_to_hls(r, g, b)
        levels = [0.1 + (0.8 / (n - 1) * i) for i in range(n)] if n > 1 else [0.5]
        random.Random(base_hex).shuffle(levels)
        return [f"#{int(colorsys.hls_to_rgb(h, lev, s)[0]*255):02x}{int(colorsys.hls_to_rgb(h, lev, s)[1]*255):02x}{int(colorsys.hls_to_rgb(h, lev, s)[2]*255):02x}" for lev in levels]

    def normalize_str(s):
        if pd.isna(s): return ""
        s = str(s).strip().lower()
        s = re.sub(r'\\s*\\([^)]*\\)', '', s).strip()
        for art in ["el ", "la ", "els ", "les ", "l'"]:
            if s.startswith(art): s = s[len(art):]; break
        s = ''.join(c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn')
        s = s.replace('ch', 'c')
        return re.sub(r'[^a-z0-9]', '', s)

    # 2. PROCESSAMENT DE DADES
    df_geno = pd.read_csv("${genotyping_file}")
    df_meta = pd.read_csv("${meta_str}", skipinitialspace=True)
    df_loc = pd.read_csv("${coordinates_file}", sep="\\\\t")
    
    df_meta.columns = [str(c).strip().upper() for c in df_meta.columns]
    df_geno['H_Subtype'] = df_geno['Subtype'].astype(str).str.extract(r'(H[0-9]+)', expand=False).fillna('Unknown')
    
    for col in ['Clade', 'Genotype']:
        df_geno[col] = df_geno[col].fillna("Unassigned").astype(str).str.strip()

    if "${params.protocol}".upper() == "HUMAN":
        all_clades = df_geno['Clade'].dropna().unique()
        def get_root(c):
            if pd.isna(c) or c in ["Unassigned", "No dataset available", "-"]: return "Unassigned"
            parts = str(c).split('.')
            if len(parts) > 3: return ".".join(parts[:3]) + "-like"
            elif len(parts) == 3 and any(str(x).startswith(str(c) + ".") for x in all_clades): return str(c) + "-like"
            return c
        df_geno['Root_Clade'] = df_geno['Clade'].apply(get_root)
    else:
        df_geno['Root_Clade'] = df_geno['Clade'].fillna('Unassigned')

    df_merged = pd.merge(df_geno, df_meta, left_on='SampleID', right_on='ID')
    df_merged['Norm_Loc'] = df_merged['LOCATION'].apply(normalize_str)
    df_loc['Norm_Pob'] = df_loc['Población'].apply(normalize_str)
    
    df = pd.merge(df_merged, df_loc, left_on='Norm_Loc', right_on='Norm_Pob', how='left')
    df['Provincia'] = df['Provincia'].fillna('Sense dades')
    df['Latitud'] = df['Latitud'].fillna(41.8204)
    df['Longitud'] = df['Longitud'].fillna(1.5412)

    if 'DATE' in df.columns:
        df['DATE'] = pd.to_datetime(df['DATE'], errors='coerce')
        iso = df['DATE'].dt.isocalendar()
        s_year = iso.year.where(iso.week >= 40, iso.year - 1)
        df['Season'] = s_year.astype(str) + "-" + (s_year + 1).astype(str)
    else:
        df['Season'] = "All Time"

    seasons = ["All Time"] + sorted([s for s in df['Season'].unique() if pd.notna(s) and s != "Unknown Season"])

    valid_views = [{'label': 'Subtipus (H)', 'id': 'subtypes', 'col': 'H_Subtype', 'filter': None}]
    for h in sorted([h for h in df['H_Subtype'].unique() if h != 'Unknown']):
        if not all(c == "-" for c in df[df['H_Subtype'] == h]['Clade'].unique()):
            valid_views.append({'label': f'Clades ({h})', 'id': f'clades_{h}', 'col': 'Root_Clade', 'filter': h})
    
    if df[(df['Clade'] == '2.3.4.4b') & (df['Genotype'] != '-')].shape[0] > 0:
        valid_views.append({'label': 'Genotips (2.3.4.4b)', 'id': 'genotypes', 'col': 'Genotype', 'filter': '2.3.4.4b_clade'})

    # 3. CONSTRUCCIÓ DEL DASHBOARD
    fig = make_subplots(
        rows=2, cols=1, specs=[[{"type": "mapbox"}], [{"type": "xy"}]],
        subplot_titles=("<b>Distribució Geogràfica de Mostres</b>", "<b>Frequència Relativa per Província (%)</b>"),
        row_heights=[0.6, 0.4], vertical_spacing=0.12
    )

    trace_registry = [] 

    for season in seasons:
        df_season = df if season == "All Time" else df[df['Season'] == season]
        for v in valid_views:
            if v['filter'] == '2.3.4.4b_clade':
                df_v = df_season[df_season['Clade'] == '2.3.4.4b'].copy()
            elif v['filter']:
                df_v = df_season[df_season['H_Subtype'] == v['filter']].copy()
            else:
                df_v = df_season.copy()
            
            if df_v.empty: continue
            df_v = df_v[df_v[v['col']] != "-"]

            labels = sorted(df_v[v['col']].unique())
            
            # Determinació del nom singular de la vista per al hover
            view_label_singular = v['label'].split(' ')[0].rstrip('s').replace('Clades', 'Clade').replace('Genotips', 'Genotip').replace('Subtipus', 'Subtipus')

            # Mapatge de colors
            if v['id'] == 'subtypes':
                color_map = {l: subtype_base_colors.get(l, '#888888') for l in labels}
            elif v['id'] == 'genotypes':
                color_map = {l: plotly_colors.qualitative.Vivid[i % len(plotly_colors.qualitative.Vivid)] for i, l in enumerate(labels)}
            else:
                base_c = subtype_base_colors.get(v['filter'], '#888888')
                shades = generate_shades(base_c, len(labels))
                color_map = {l: shades[i] for i, l in enumerate(labels)}

            if is_colorblind:
                color_map = {l: okabe_ito_colors[i % len(okabe_ito_colors)] for i, l in enumerate(labels)}

            # A) Mapa
            df_map = df_v.groupby(['LOCATION', 'Latitud', 'Longitud', v['col']]).size().reset_index(name='Count')
            for label in labels:
                sub = df_map[df_map[v['col']] == label]
                
                # Customdata per al mapa: [Població, Casos]
                cdata_map = sub[['LOCATION', 'Count']].values
                
                fig.add_trace(go.Scattermapbox(
                    lat=sub['Latitud'], lon=sub['Longitud'], mode='markers',
                    marker=dict(size=sub['Count'] * 6 + 7, color=color_map[label], opacity=0.8),
                    name=str(label), legendgroup=v['id'], showlegend=True,
                    customdata=cdata_map,
                    hovertemplate=(
                        "<b>Població:</b> %{customdata[0]}<br>" +
                        "<b>" + view_label_singular + ":</b> " + str(label) + "<br>" +
                        "<b>Casos:</b> %{customdata[1]}<extra></extra>"
                    ),
                    visible=False
                ), row=1, col=1)
                trace_registry.append({'season': season, 'view': v['id']})

            # B) Barres apilades al 100%
            prov_totals = df_v.groupby('Provincia').size()
            df_bar = df_v.groupby(['Provincia', v['col']]).size().reset_index(name='Count')
            df_bar['Pct'] = df_bar.apply(lambda r: (r['Count'] / prov_totals[r['Provincia']]) * 100, axis=1)
            
            for label in labels:
                sub = df_bar[df_bar[v['col']] == label].copy()
                # Afegim el total de la província per al hover (format x/y)
                sub['TotalProv'] = sub['Provincia'].map(prov_totals)
                
                # Customdata per a barres: [Count, TotalProv]
                cdata_bar = sub[['Count', 'TotalProv']].values
                
                fig.add_trace(go.Bar(
                    x=sub['Provincia'], y=sub['Pct'], name=str(label),
                    marker_color=color_map[label], legendgroup=v['id'], showlegend=False,
                    text=sub['Count'].astype(str), textposition='inside',
                    customdata=cdata_bar,
                    hovertemplate=(
                        "<b>" + view_label_singular + ":</b> " + str(label) + "<br>" +
                        "<b>Occurrence:</b> %{customdata[0]}/%{customdata[1]} samples (%{y:.1f}%)<extra></extra>"
                    ),
                    visible=False
                ), row=2, col=1)
                trace_registry.append({'season': season, 'view': v['id']})

    # 4. MENÚS NATIUS I LAYOUT
    def get_vis_mask(s, v_id):
        return [True if (r['season'] == s and r['view'] == v_id) else False for r in trace_registry]

    initial_vis = get_vis_mask(seasons[0], 'subtypes')
    for i, vis in enumerate(initial_vis):
        fig.data[i].visible = vis

    s_btns = [dict(label=s, method="update", args=[{"visible": get_vis_mask(s, 'subtypes')}]) for s in seasons]
    v_btns = [dict(label=v['label'], method="update", args=[{"visible": get_vis_mask(seasons[0], v['id'])}]) for v in valid_views]

    fig.update_layout(
        mapbox_style="carto-positron", mapbox=dict(center=dict(lat=41.7, lon=1.8), zoom=6.8),
        barmode='stack', height=1000,
        title=dict(text="<b>Monitoratge Geolocalitzat de FluTyper</b>", x=0.5, font=dict(size=22)),
        updatemenus=[
            dict(buttons=s_btns, x=0.1, y=1.12, xanchor="left", direction="down"),
            dict(buttons=v_btns, x=0.9, y=1.12, xanchor="right", direction="down")
        ],
        margin=dict(t=150, b=50, l=80, r=80)
    )
    fig.update_yaxes(title_text="Freqüència (%)", range=[0, 100], row=2, col=1)

    fig.write_html("GeographicReport.html")
    """
}