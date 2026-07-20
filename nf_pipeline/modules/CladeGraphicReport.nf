#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CladeGraphicReport {
    errorStrategy 'ignore'
    debug true

    input:
    path(genotyping_file)
    path(metadata_file) 

    output:
    path("CladeGraphicReport.html"), emit: report
    path("CladeEvolutionReport.html"), emit: evolution_report, optional: true

    script:
    def meta_str = metadata_file ? metadata_file.toString() : ""
    """
    #!/usr/bin/env python3
    import os
    import json
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from plotly import colors
    import colorsys
    import random

    def get_contrast_text_color(hex_str):
        hex_str = str(hex_str).lstrip('#')
        if len(hex_str) == 6:
            try:
                r, g, b = tuple(int(hex_str[i:i+2], 16) for i in (0, 2, 4))
                luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                return '#000000' if luminance > 0.55 else '#ffffff'
            except ValueError:
                return '#ffffff'
        return '#ffffff'

    is_colorblind = "${params.colorblind}".lower() == "true"

    okabe_ito_colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999', '#000000']
    
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
            r, g, b = [int(clean_hex[i:i+2], 16) / 255.0 for i in (0, 2, 4)]
        except ValueError:
            return ['#888888'] * n
        hue, lightness, saturation = colorsys.rgb_to_hls(r, g, b)
        step = 0.8 / (n - 1)
        brightness_levels = [0.1 + (step * i) for i in range(n)]
        random.Random(base_hex).shuffle(brightness_levels)
        shades = []
        for new_lightness in brightness_levels:
            new_r, new_g, new_b = colorsys.hls_to_rgb(hue, new_lightness, saturation)
            shades.append(f"#{int(new_r * 255):02x}{int(new_g * 255):02x}{int(new_b * 255):02x}")
        return shades

    genotyping_df = pd.read_csv("${genotyping_file}")
    genotyping_df['H_Subtype'] = genotyping_df['Subtype'].astype(str).str.extract(r'(H\\d+)', expand=False)
    
    if "${params.protocol}".upper() == "HUMAN":
        genotyping_df['H_Subtype'] = genotyping_df['H_Subtype'].replace({'H1': 'A(H1)pdm09', 'H3': 'A(H3)'})

    # Ensure Clade columns are treated as categorical/string and clean them
    for col in ['Clade', 'Genotype', 'Sub-genotype']:
        genotyping_df[col] = genotyping_df[col].fillna("Unassigned").astype(str).str.strip()
        
    # Safely handle Optional Metadata
    genotyping_df['Age_Group'] = 'Sense dades'
    genotyping_df['Sex'] = 'Sense dades'
    genotyping_df['Season'] = 'Unknown Season'
    
    metadata_param = "${meta_str}"
    if metadata_param and os.path.isfile(metadata_param):
        df_meta = pd.read_csv(metadata_param, skipinitialspace=True)
        # Normalize columns to uppercase to avoid case-sensitivity issues
        df_meta.columns = [c.upper() for c in df_meta.columns]
        
        if 'ID' in df_meta.columns:
            df_meta = df_meta.dropna(subset=['ID']).drop_duplicates(subset=['ID'], keep='first')
            
            # Map optional demographics
            if 'AGE GROUP' in df_meta.columns:
                df_meta = df_meta.rename(columns={'AGE GROUP': 'AGE_GROUP'})
            
            if 'AGE_GROUP' in df_meta.columns:
                df_meta['AGE_GROUP'] = df_meta['AGE_GROUP'].fillna('Sense dades')
                
            if 'SEX' in df_meta.columns:
                df_meta['SEX'] = df_meta['SEX'].fillna('Sense dades')
                
            # Merge Metadata
            genotyping_df = genotyping_df.merge(df_meta, left_on='SampleID', right_on='ID', how='left')
            genotyping_df['Age_Group'] = genotyping_df['AGE_GROUP'].fillna('Sense dades') if 'AGE_GROUP' in genotyping_df.columns else 'Sense dades'
            genotyping_df['Sex'] = genotyping_df['SEX'].fillna('Sense dades') if 'SEX' in genotyping_df.columns else 'Sense dades'
            
            # Calculate Season and WEEK if Date exists
            if 'DATE' in df_meta.columns:
                genotyping_df['DATE'] = pd.to_datetime(genotyping_df['DATE'], errors='coerce')
                iso_cal = genotyping_df['DATE'].dt.isocalendar()
                s_year = iso_cal.year.where(iso_cal.week >= 40, iso_cal.year - 1)
                genotyping_df['Season'] = s_year.apply(lambda x: f"{int(x)}-{int(x)+1}" if pd.notna(x) else "Unknown Season")
                # Create WEEK column as the Monday of each ISO week (needed for evolution charts)
                genotyping_df['WEEK'] = genotyping_df['DATE'].dt.to_period('W').apply(
                    lambda p: p.start_time if pd.notna(p) else pd.NaT
                )

    invalid_clade_values = {'-', 'Unassigned', 'No dataset available'}

    # Clade logic refinement for Avian
    if "${params.protocol}".upper() == "AVIAN":
        genotyping_df['Root_Clade'] = genotyping_df['Clade']
    else:
        # Existing logic for human/other
        all_clades = genotyping_df['Clade'].dropna().unique()
        def get_root_clade(c):
            if pd.isna(c) or c in ["Unassigned", "No dataset available"]: return c
            parts = str(c).split('.')
            return ".".join(parts[:3]) + "-like" if len(parts) > 3 else c
        genotyping_df['Root_Clade'] = genotyping_df['Clade'].apply(get_root_clade)

    raw_seasons = genotyping_df['Season'].dropna().unique()
    seasons_pie = sorted([s for s in raw_seasons if s not in ["All Time", "Unknown Season"]], reverse=True)
    seasons_evo = seasons_pie.copy()
    seasons_evo.append("All Time")

    age_order = {'0-2': 0, '3-4': 1, '5-14': 2, '15-65': 3, '>65': 4}
    age_groups = ['All'] + sorted([a for a in genotyping_df['Age_Group'].unique() if str(a).strip() not in ['nan', '', 'None', 'Sense dades']], key=lambda x: age_order.get(str(x).strip(), 99))
    sexs = ['All'] + sorted([g for g in genotyping_df['Sex'].unique() if str(g).strip() not in ['nan', '', 'None', 'Sense dades']])

    unique_h_subtypes = sorted(genotyping_df['H_Subtype'].dropna().unique())
    valid_h_subtypes = []
    for h in unique_h_subtypes:
        sub = genotyping_df[genotyping_df['H_Subtype'] == h]
        # Only count rows that have a real clade assignment
        valid_clade_rows = sub[~sub['Clade'].isin(invalid_clade_values)]
        if valid_clade_rows.empty:
            continue
        valid_h_subtypes.append(h)
            
    include_genotype_chart = genotyping_df[(genotyping_df['Clade'] == '2.3.4.4b') & (genotyping_df['Genotype'] != '-')].shape[0] > 0
    total_charts = 1 + len(valid_h_subtypes) + (1 if include_genotype_chart else 0)
    
    global_color_map = {}
    h_labels = sorted([h for h in genotyping_df['H_Subtype'].dropna().unique() if h != '-'])
    h_color_map = {}
    h_color_map["Unassigned"] = '#d3d3d3'
    for h in h_labels:
        h_color_map[h] = okabe_ito_colors[len(h_color_map) % len(okabe_ito_colors)] if is_colorblind else subtype_base_colors.get(str(h), '#888888')
    global_color_map['h_subtypes'] = h_color_map
    
    for h in valid_h_subtypes:
        base_c = subtype_base_colors.get(str(h), '#888888')
        h_clades = sorted([c for c in genotyping_df[genotyping_df['H_Subtype'] == h]['Root_Clade'].dropna().unique() if c not in ["Unassigned", "No dataset available", "-"]])
        clade_palette = okabe_ito_colors if is_colorblind else generate_shades(base_c, len(h_clades) + 2)
        clade_map = {}
        for idx, clade in enumerate(h_clades): clade_map[clade] = clade_palette[idx % len(clade_palette)]
        clade_map["Unassigned"] = '#d3d3d3'
        clade_map["Others"] = clade_palette[-1]
        global_color_map[f'clades_{h}'] = clade_map
    
    if include_genotype_chart:
        g_labels = sorted([g for g in genotyping_df[genotyping_df['Clade'] == '2.3.4.4b']['Genotype'].dropna().unique() if g not in ["-", "None", "", "nan", "Unassigned"]])
        genotype_palette = okabe_ito_colors if is_colorblind else colors.qualitative.Vivid
        genotype_map = {}
        for idx, genotype in enumerate(g_labels): genotype_map[genotype] = genotype_palette[idx % len(genotype_palette)]
        global_color_map['genotypes'] = genotype_map

    season_options_evo = "".join([f'<option value="{s}">{s}</option>' for s in seasons_evo])
    season_options_pie = "".join([f'<option value="{s}">{s}</option>' for s in seasons_pie])
    age_options = "".join([f'<option value="{a}">{a}</option>' for a in age_groups])
    sex_options = "".join([f'<option value="{g}">{g}</option>' for g in sexs])

    def generate_ui_html(report_type):
        opts = season_options_evo if report_type == 'evo' else season_options_pie
        return f'''
        <div style="display:flex; gap:15px; justify-content:center; margin-top:20px; font-family:Arial; background:#f9f9f9; padding:15px; border-radius:8px; border:1px solid #ddd; width:fit-content; margin-left:auto; margin-right:auto; box-shadow:0 4px 10px rgba(0,0,0,0.05);">
            <div style="min-width:150px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">SEASON</label><br>
                <select id="sel_season_{report_type}" style="padding:6px; border-radius:4px; width:100%; border:1px solid #ccc; background:white;">{opts}</select>
            </div>
            <div style="min-width:120px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">AGE GROUP</label><br>
                <select id="sel_age_{report_type}" style="padding:6px; border-radius:4px; width:100%; border:1px solid #ccc; background:white;">{age_options}</select>
            </div>
            <div style="min-width:120px;">
                <label style="font-size:10px; font-weight:bold; color:#666;">SEX</label><br>
                <select id="sel_sex_{report_type}" style="padding:6px; border-radius:4px; width:100%; border:1px solid #ccc; background:white;">{sex_options}</select>
            </div>
        </div>
        '''

    # CLADE EVOLUTION BAR CHARTS
    if 'WEEK' in genotyping_df.columns and not genotyping_df['WEEK'].isna().all() and len(valid_h_subtypes) > 0:
        df_all = genotyping_df.dropna(subset=['WEEK']).copy()
        
        if "${params.protocol}".upper() == "HUMAN":
            subplot_titles_evo = ["<b>Influenza A Subtype Weekly Evolution</b>"] + [f"<b>Clade Evolution for {h}</b>" for h in valid_h_subtypes]
        else:
            subplot_titles_evo = ["<b>H Subtype Weekly Evolution</b>"] + [f"<b>Clade Evolution for {h}</b>" for h in valid_h_subtypes]
        if include_genotype_chart:
            subplot_titles_evo.append("<b>Genotype Evolution (Clade 2.3.4.4b)</b>")
        
        fig_evo = make_subplots(rows=total_charts, cols=1, shared_xaxes=True, vertical_spacing=0.08, subplot_titles=subplot_titles_evo)

        def add_stacked_bars(df_subset, group_col, row_num, detail_col=None, color_key='', meta_dict=None):
            if df_subset.empty: return
            weekly_totals = df_subset.groupby('WEEK').size().rename('TotalWeek')
            if detail_col:
                detail_counts = df_subset.groupby(['WEEK', group_col, detail_col]).size().reset_index(name='DetailCount')
                detail_counts = detail_counts[~detail_counts[detail_col].isin(["-", "None", "", "nan", "Unassigned"])]
                if not detail_counts.empty:
                    detail_counts['DetailStr'] = "- " + detail_counts[detail_col].astype(str) + ": " + detail_counts['DetailCount'].astype(str)
                    hover_details = detail_counts.groupby(['WEEK', group_col])['DetailStr'].apply(lambda x: "<br>".join(x)).reset_index(name='Hover_Details')
                else:
                    hover_details = pd.DataFrame(columns=['WEEK', group_col, 'Hover_Details'])
            else:
                hover_details = pd.DataFrame(columns=['WEEK', group_col, 'Hover_Details'])

            group_counts = df_subset.groupby(['WEEK', group_col]).size().reset_index(name='Count')
            group_counts = group_counts.merge(weekly_totals, on='WEEK')
            group_counts = group_counts.merge(hover_details, on=['WEEK', group_col], how='left') if not hover_details.empty else group_counts.assign(Hover_Details=float('nan'))
            group_counts['Hover_Details'] = group_counts['Hover_Details'].fillna("")
            group_counts['Pct'] = (group_counts['Count'] / group_counts['TotalWeek']) * 100

            def make_hover_extra(row):
                val = str(row[group_col])
                details = row['Hover_Details']
                if not details: return ""
                if group_col == 'Final_Label' and val.endswith("-like"): return "<br><br><b>Clade Breakdown:</b><br>" + details
                elif group_col == 'Genotype': return "<br><br><b>Sub-genotypes Breakdown:</b><br>" + details
                return ""

            group_counts['Hover_Extra'] = group_counts.apply(make_hover_extra, axis=1)
            color_map = global_color_map.get(color_key, {})

            for g in group_counts[group_col].unique():
                g_df = group_counts[group_counts[group_col] == g]
                color = color_map.get(str(g), '#888888')
                display_group_col = "Group" if group_col == "Final_Label" else group_col.replace("_", " ")
                
                fig_evo.add_trace(go.Bar(
                    x=g_df['WEEK'], y=g_df['Pct'], name=str(g), showlegend=True, text=g_df['Count'], textposition='inside', insidetextanchor='middle', textangle=0,
                    textfont=dict(color=get_contrast_text_color(color), size=14, family='Arial'), marker_color=color, legend=f"legend{row_num}" if row_num > 1 else "legend",
                    customdata=g_df[['Count', 'TotalWeek', 'Pct', 'Hover_Extra']],
                    hovertemplate=("<b>Week:</b> %{x|%V, %Y}<br><b>" + display_group_col + ":</b> " + str(g) + "<br><b>Occurrence:</b> %{customdata[0]}/%{customdata[1]} samples (%{customdata[2]:.1f}%)%{customdata[3]}<extra></extra>"),
                    visible=False, meta=meta_dict, width=504800000
                ), row=row_num, col=1)
            fig_evo.update_yaxes(title_text="Frequency (%)", range=[0, 100], row=row_num, col=1)

        for age in age_groups:
            for sex in sexs:
                df_view = df_all.copy()
                if age != 'All': df_view = df_view[df_view['Age_Group'] == age]
                if sex != 'All': df_view = df_view[df_view['Sex'] == sex]
                meta_dict = {'age': age, 'sex': sex}
                
                add_stacked_bars(df_view, 'H_Subtype', 1, color_key='h_subtypes', meta_dict=meta_dict)
                current_row = 2
                for h in valid_h_subtypes:
                    sub_df = df_view[df_view['H_Subtype'] == h].copy()
                    sub_df = sub_df[~sub_df['Clade'].isin(invalid_clade_values)]
                    sub_df['Final_Label'] = sub_df['Root_Clade']
                    add_stacked_bars(sub_df, 'Final_Label', current_row, detail_col='Clade', color_key=f'clades_{h}', meta_dict=meta_dict)
                    current_row += 1

                if include_genotype_chart:
                    sub_df = df_view[(df_view['Clade'] == '2.3.4.4b') & (df_view['Genotype'] != '-')].copy()
                    add_stacked_bars(sub_df, 'Genotype', current_row, detail_col='Sub-genotype', color_key='genotypes', meta_dict=meta_dict)

        season_ranges = {}
        for season_val in sorted(df_all['Season'].dropna().unique()):
            if season_val in ["Unknown Season", "All Time"]: continue
            try:
                y1, y2 = int(season_val.split('-')[0]), int(season_val.split('-')[1])
                season_ranges[season_val] = [pd.to_datetime(f'{y1}-W40-1', format='%G-W%V-%u').strftime('%Y-%m-%d'), pd.to_datetime(f'{y2}-W39-7', format='%G-W%V-%u').strftime('%Y-%m-%d')]
            except Exception:
                s_data = df_all[df_all['Season'] == season_val]
                if not s_data.empty: season_ranges[season_val] = [s_data['WEEK'].min().strftime('%Y-%m-%d'), s_data['WEEK'].max().strftime('%Y-%m-%d')]

        legends_layout = {}
        h_domain = (1.0 - (total_charts - 1) * 0.08) / total_charts
        for i in range(1, total_charts + 1):
            legends_layout[f"legend{i}" if i > 1 else "legend"] = dict(y=1.0 - (i - 1) * (h_domain + 0.08), yanchor="top", x=1.02, xanchor="left", title_text=subplot_titles_evo[i-1].replace("<b>","").replace("</b>",""), tracegroupgap=0)

        initial_lbl_evo = seasons_evo[0] if seasons_evo[0] == "All Time" else f"Season {seasons_evo[0]}"
        fig_evo.update_layout(
            barmode='stack', title=dict(text=f"<b>Evolution of Subtypes and Clades - {initial_lbl_evo}</b>", x=0.45 if "${params.protocol}" == "AVIAN" else 0.465, y=0.98, xanchor="center", yanchor="top", font=dict(size=24)),
            height=450 * total_charts, hovermode="closest", margin=dict(t=120, b=80, l=80, r=200), **legends_layout
        )
        fig_evo.update_xaxes(tickformat="Week %V<br>%Y", showticklabels=True)

        js_evo = f'''
        <script>
            var seasonRanges = {json.dumps(season_ranges)};
            var totalCharts = {total_charts};
            function updateEvoPlot() {{
                var s = document.getElementById('sel_season_evo').value;
                var a = document.getElementById('sel_age_evo').value;
                var g = document.getElementById('sel_sex_evo').value;
                var plotDivs = document.getElementsByClassName('plotly-graph-div');
                if (plotDivs.length === 0) return;
                var plotDiv = plotDivs[0];
                var update = {{ visible: [] }};
                
                for (var i = 0; i < plotDiv.data.length; i++) {{
                    var meta = plotDiv.data[i].meta;
                    if (meta && meta.age === a && meta.sex === g) {{
                        var hasData = false;
                        if (s === 'All Time') {{ hasData = plotDiv.data[i].x && plotDiv.data[i].x.length > 0; }}
                        else if (seasonRanges[s]) {{
                            var sStart = new Date(seasonRanges[s][0]), sEnd = new Date(seasonRanges[s][1]);
                            if (plotDiv.data[i].x) {{
                                for (var j = 0; j < plotDiv.data[i].x.length; j++) {{
                                    var xDate = new Date(plotDiv.data[i].x[j]);
                                    if (xDate >= sStart && xDate <= sEnd) {{ hasData = true; break; }}
                                }}
                            }}
                        }}
                        update.visible.push(hasData ? true : false);
                    }} else {{ update.visible.push(false); }}
                }}
                Plotly.restyle(plotDiv, update);
                
                var layoutUpdate = {{}};
                layoutUpdate['title.text'] = '<b>Evolution of Subtypes and Clades - ' + (s === 'All Time' ? s : 'Season ' + s) + '</b>';
                for (var i = 0; i < totalCharts; i++) {{
                    var ax = i === 0 ? 'xaxis' : 'xaxis' + (i + 1);
                    if (s !== 'All Time' && seasonRanges[s]) {{
                        layoutUpdate[ax + '.range'] = seasonRanges[s];
                        layoutUpdate[ax + '.autorange'] = false;
                        layoutUpdate[ax + '.fixedrange'] = true;
                    }} else {{
                        layoutUpdate[ax + '.autorange'] = true;
                        layoutUpdate[ax + '.fixedrange'] = false;
                    }}
                }}
                Plotly.relayout(plotDiv, layoutUpdate);
            }}
            document.getElementById('sel_season_evo').addEventListener('change', updateEvoPlot);
            document.getElementById('sel_age_evo').addEventListener('change', updateEvoPlot);
            document.getElementById('sel_sex_evo').addEventListener('change', updateEvoPlot);
            window.addEventListener('load', updateEvoPlot);
        </script>
        '''
        with open("CladeEvolutionReport.html", "w") as f: f.write(fig_evo.to_html(include_plotlyjs='cdn', full_html=True).replace('<body>', '<body>\\n' + generate_ui_html('evo')).replace('</body>', js_evo + '\\n</body>'))

    # PIE CHARTS
    cols, rows = total_charts, 1

    if "${params.protocol}".upper() == "HUMAN":
        subplot_titles = ["<b>Influenza A Subtype Distribution</b>"] + [f"<b>Clade Distribution for {h}</b>" for h in valid_h_subtypes]
    else:
        subplot_titles = ["<b>H Subtype Distribution</b>"] + [f"<b>Clade Distribution for {h}</b>" for h in valid_h_subtypes]
    if include_genotype_chart: subplot_titles.append("<b>Genotype Distribution (Clade 2.3.4.4b)</b>")
    fig = make_subplots(rows=rows, cols=cols, specs=[[{"type": "domain"} for _ in range(cols)]], subplot_titles=subplot_titles, horizontal_spacing=0.05)

    for season in seasons_pie:
        for age in age_groups:
            for sex in sexs:
                meta_dict = {'season': season, 'age': age, 'sex': sex}
                df_view = genotyping_df[genotyping_df['Season'] == season].copy()
                if age != 'All': df_view = df_view[df_view['Age_Group'] == age]
                if sex != 'All': df_view = df_view[df_view['Sex'] == sex]

                h_counts = df_view['H_Subtype'].value_counts().reset_index()
                h_counts.columns = ['Label', 'Count']
                if h_counts.empty: fig.add_trace(go.Pie(labels=["No Data"], values=[1], name="H Subtypes", textinfo='none', hoverinfo='none', marker=dict(colors=['#f0f0f0']), visible=False, meta=meta_dict), row=1, col=1)
                else:
                    h_counts['Text'] = h_counts['Count'].astype(str) + '/' + str(h_counts['Count'].sum())
                    fig.add_trace(go.Pie(
                        labels=h_counts['Label'], values=h_counts['Count'], name="H Subtypes", text=h_counts['Text'], texttemplate='<b>%{label}</b><br><b>%{text}</b><br><b>%{percent}</b>',
                        textposition='outside', insidetextorientation='horizontal', rotation=270, automargin=True, hole=0.35, marker=dict(colors=[global_color_map['h_subtypes'].get(str(lbl), '#888888') for lbl in h_counts['Label']], line=dict(color='#ffffff', width=2)),
                        hoverlabel=dict(font_size=14), hovertemplate='<b>H Subtype:</b> %{label}<br><b>Count:</b> %{text}<br><b>Percentage:</b> %{percent}<extra></extra>', visible=False, meta=meta_dict
                    ), row=1, col=1)

                for i, h in enumerate(valid_h_subtypes):
                    c_col = i + 2
                    sub_df = df_view[df_view['H_Subtype'] == h]
                    sub_df = sub_df[~sub_df['Clade'].isin(invalid_clade_values)]
                    if len(sub_df) == 0: fig.add_trace(go.Pie(labels=["No Data"], values=[1], name=str(h), textinfo='none', hoverinfo='none', marker=dict(colors=['#f0f0f0']), visible=False, meta=meta_dict), row=1, col=c_col); continue
                    
                    orig_counts = sub_df.groupby(['Clade', 'Root_Clade']).size().reset_index(name='Count')
                    orig_counts['Hover_Detail'] = orig_counts.apply(lambda x: f"- {x['Clade']}: {x['Count']}/{len(sub_df)} ({x['Count']/len(sub_df):.1%})" if x['Count'] > 0 else "", axis=1)
                    root_grouped = orig_counts.groupby('Root_Clade').agg(Root_Count=('Count', 'sum'), Root_Hover_Details=('Hover_Detail', lambda x: "<br>".join([d for d in x if d]))).reset_index()
                    root_grouped['Final_Label'] = root_grouped.apply(lambda x: "Others" if x['Root_Count']/len(sub_df) < 0.01 and str(x['Root_Clade']) not in ["Unassigned", "No dataset available"] else x['Root_Clade'], axis=1)
                    final_grouped = root_grouped.groupby('Final_Label').agg(Final_Count=('Root_Count', 'sum'), Final_Hover_Details=('Root_Hover_Details', lambda x: "<br>".join([d for d in x if d]))).reset_index()
                    final_grouped['Hover_Extra'] = final_grouped.apply(lambda x: "<br><br><b>Clade Breakdown:</b><br>" + x['Final_Hover_Details'] if str(x['Final_Label']).endswith("-like") or str(x['Final_Label']) == "Others" else "", axis=1)
                    final_grouped['Text'] = final_grouped['Final_Count'].astype(str) + '/' + str(len(sub_df))
                    
                    fig.add_trace(go.Pie(
                        labels=final_grouped['Final_Label'], values=final_grouped['Final_Count'], name=str(h), text=final_grouped['Text'], texttemplate='<b>%{label}</b><br><b>%{text}</b><br><b>%{percent}</b>',
                        textposition='outside', rotation=270, automargin=True, insidetextorientation='horizontal', hole=0.35, marker=dict(colors=[global_color_map[f'clades_{h}'].get(str(lbl), '#888888') for lbl in final_grouped['Final_Label']], line=dict(color='#ffffff', width=2)),
                        hoverlabel=dict(font_size=14, align='left'), customdata=final_grouped['Hover_Extra'], hovertemplate='<b>Group:</b> %{label}<br><b>Total Count:</b> %{text}<br><b>Group Percentage:</b> %{percent}%{customdata}<extra></extra>', visible=False, meta=meta_dict
                    ), row=1, col=c_col)

                if include_genotype_chart:
                    sub_df = df_view[(df_view['Clade'] == '2.3.4.4b') & (df_view['Genotype'] != '-')]
                    if len(sub_df) == 0: fig.add_trace(go.Pie(labels=["No Data"], values=[1], name="Genotypes", textinfo='none', hoverinfo='none', marker=dict(colors=['#f0f0f0']), visible=False, meta=meta_dict), row=1, col=total_charts)
                    else:
                        orig_counts = sub_df.groupby(['Genotype', 'Sub-genotype']).size().reset_index(name='Count')
                        orig_counts['Hover_Detail'] = orig_counts.apply(lambda x: f"- {x['Sub-genotype']}: {x['Count']}/{len(sub_df)} ({x['Count']/len(sub_df):.1%})" if str(x['Sub-genotype']) not in ["-", "None", "", "nan", "Unassigned"] else "", axis=1)
                        root_grouped = orig_counts.groupby('Genotype').agg(Final_Count=('Count', 'sum'), Hover_Details=('Hover_Detail', lambda x: "<br>".join([d for d in x if d]))).reset_index()
                        root_grouped['Hover_Extra'] = root_grouped.apply(lambda x: "<br><br><b>Sub-genotypes Breakdown:</b><br>" + x['Hover_Details'] if x['Hover_Details'] else "", axis=1)
                        root_grouped['Text'] = root_grouped['Final_Count'].astype(str) + '/' + str(len(sub_df))

                        fig.add_trace(go.Pie(
                            labels=root_grouped['Genotype'], values=root_grouped['Final_Count'], name="Genotypes 2.3.4.4b", text=root_grouped['Text'], texttemplate='<b>%{label}</b><br><b>%{text}</b><br><b>%{percent}</b>',
                            textposition='outside', rotation=270, automargin=True, insidetextorientation='horizontal', hole=0.35, marker=dict(colors=[global_color_map['genotypes'].get(str(lbl), '#888888') for lbl in root_grouped['Genotype']], line=dict(color='#ffffff', width=2)),
                            hoverlabel=dict(font_size=14, align='left'), customdata=root_grouped['Hover_Extra'], hovertemplate='<b>Genotype:</b> %{label}<br><b>Total Count:</b> %{text}<br><b>Percentage:</b> %{percent}%{customdata}<extra></extra>', visible=False, meta=meta_dict
                        ), row=1, col=total_charts)

    if 'annotations' in fig['layout']:
        for annotation in fig['layout']['annotations']: annotation['y'] += 0.1
    
    fig.update_layout(title=dict(text=f"<b>Subtype and Clade Report - Season {seasons_pie[0] if seasons_pie else 'No Data'}</b>", x=0.5, y=0.98, xanchor="center", yanchor="top", font=dict(size=24)), height=680, showlegend=False, hovermode="closest", margin=dict(t=180, b=80, l=40, r=40), uniformtext=dict(minsize=10, mode='show'))

    js_pie = f'''
    <script>
        function updatePiePlot() {{
            var s = document.getElementById('sel_season_pie').value;
            var a = document.getElementById('sel_age_pie').value;
            var g = document.getElementById('sel_sex_pie').value;
            var plotDivs = document.getElementsByClassName('plotly-graph-div');
            if (plotDivs.length === 0) return;
            var plotDiv = plotDivs[0];
            var update = {{ visible: [] }};
            for (var i = 0; i < plotDiv.data.length; i++) {{
                var meta = plotDiv.data[i].meta;
                update.visible.push(meta && meta.season === s && meta.age === a && meta.sex === g);
            }}
            Plotly.restyle(plotDiv, update);
            Plotly.relayout(plotDiv, {{ 'title.text': '<b>Subtype and Clade Report - Season ' + s + '</b>' }});
        }}
        document.getElementById('sel_season_pie').addEventListener('change', updatePiePlot);
        document.getElementById('sel_age_pie').addEventListener('change', updatePiePlot);
        document.getElementById('sel_sex_pie').addEventListener('change', updatePiePlot);
        window.addEventListener('load', updatePiePlot);
    </script>
    '''
    with open("CladeGraphicReport.html", "w") as f: f.write(fig.to_html(include_plotlyjs='cdn', full_html=True).replace('<body>', '<body>\\n' + generate_ui_html('pie')).replace('</body>', js_pie + '\\n</body>'))
    """
}