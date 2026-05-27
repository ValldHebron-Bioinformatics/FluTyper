#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CladeGraphicReport {
    errorStrategy 'ignore'
    debug true

    input:
    path(genotyping_file)
    path(metadata_file) // Nextflow safely ignores this if it receives an empty collection

    output:
    path("CladeGraphicReport.html"), emit: report
    path("CladeEvolutionReport.html"), emit: evolution_report, optional: true

    script:
    // Evaluate the Groovy variable. If it's an empty list, pass an empty string to Python.
    def meta_str = metadata_file ? metadata_file.toString() : ""
    """
    #!/usr/bin/env python3
    import os
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from plotly import colors
    import colorsys
    import random

    # Define text color (black or white) based on background color for readability
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

    # Okabe-Ito Colors Palette (colorblind-friendly)
    okabe_ito_colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999', '#000000']
    
    # 18 Base Colors mapping for H Subtypes
    subtype_base_colors = {
        'H1': '#1F77B4','H2': '#D62728','H3': '#FF7F0E','H4': '#2CA02C','H5': '#9467BD','H6': '#8C564B', 
        'H7': '#E377C2','H8': '#7F7F7F','H9': '#BCBD22','H10': '#17BECF','H11': '#393B79','H12': '#637939',
        'H13': '#8C6D31','H14': '#843C39','H15': '#7B4173','H16': '#5254A3','H17': '#8CA252','H18': '#BD9E39' 
    }

    # Generate a randomized monochromatic palette
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

    # Dataframe preparation
    genotyping_df = pd.read_csv("${genotyping_file}")
    genotyping_df['H_Subtype'] = genotyping_df['Subtype'].astype(str).str.extract(r'(H\\d+)', expand=False)

        
    # Standardize Clade, Genotype, and Sub-genotype columns
    for col in ['Clade', 'Genotype', 'Sub-genotype']:
        genotyping_df[col] = genotyping_df[col].fillna("Unassigned").astype(str).str.strip()
        genotyping_df.loc[genotyping_df[col].str.lower().str.contains('unassigned'), col] = 'Unassigned'
        
    
    genotyping_df.loc[genotyping_df['Clade'] == '-', 'Clade'] = 'No dataset available'


    # Dynamic Root Grouping Logic
    if "${params.protocol}".upper() == "AVIAN":
        genotyping_df['Root_Clade'] = genotyping_df['Clade']
    else:
        all_clades = genotyping_df['Clade'].dropna().unique()
        def get_root_clade(c):
            if pd.isna(c) or c in ["Unassigned", "No dataset available"]: return c
            parts = str(c).split('.')
            if len(parts) > 3: return ".".join(parts[:3]) + "-like"
            elif len(parts) == 3 and any(str(x).startswith(str(c) + ".") for x in all_clades): return str(c) + "-like"
            return c
        genotyping_df['Root_Clade'] = genotyping_df['Clade'].apply(get_root_clade)

    # Inject the safely evaluated Groovy string
    metadata_param = "${meta_str}"

    # Process metadata purely inside Python
    if metadata_param and os.path.isfile(metadata_param):
        df_meta = pd.read_csv(metadata_param, skipinitialspace=True)
        df_meta['DATE'] = pd.to_datetime(df_meta['DATE'], format='%Y-%m-%d')
        df_meta['WEEK'] = df_meta['DATE'].dt.to_period('W').dt.to_timestamp()
        
        # Merge utilizing SampleID
        genotyping_df = genotyping_df.merge(df_meta[['ID', 'DATE', 'WEEK']], left_on='SampleID', right_on='ID', how='left')

    # Calculate Season based on the merged DATE column
    if 'DATE' in genotyping_df.columns:
        # Force conversion to datetime to catch any empty strings or garbage as pure nulls
        genotyping_df['DATE'] = pd.to_datetime(genotyping_df['DATE'], errors='coerce')
        
        iso_cal = genotyping_df['DATE'].dt.isocalendar()
        s_year = iso_cal.year.where(iso_cal.week >= 40, iso_cal.year - 1)
        valid_dates = s_year.notna()
        
        if valid_dates.any():
            # If we have some valid dates, initialize with pd.NA so the missing ones vanish
            genotyping_df['Season'] = pd.NA
            genotyping_df.loc[valid_dates, 'Season'] = (
                s_year[valid_dates].astype(int).astype(str) + "-" + 
                (s_year[valid_dates].astype(int) + 1).astype(str)
            )
        else:
            # The column exists but absolutely zero valid dates were found
            genotyping_df['Season'] = "All Time"
    else:
        genotyping_df['Season'] = "All Time"

    # Extract unique seasons, automatically wiping out the pd.NA nulls
    seasons = sorted(genotyping_df['Season'].dropna().unique())
    
    # Ultimate failsafe: if the list is totally empty, default the view to All Time
    if len(seasons) == 0:
        genotyping_df['Season'] = "All Time"
        seasons = ["All Time"]

    seasons = sorted(genotyping_df['Season'].dropna().unique())
    unique_h_subtypes = sorted(genotyping_df['H_Subtype'].dropna().unique())
    
    # Filter out subtypes 
    valid_h_subtypes = []
    for h in unique_h_subtypes:
        clades_for_h = genotyping_df[genotyping_df['H_Subtype'] == h]['Clade'].unique()
        if not (len(clades_for_h) == 1 and clades_for_h[0] == "No dataset available"):
            valid_h_subtypes.append(h)
    
    # Check if there are any valid Genotypes for clade 2.3.4.4b globally
    include_genotype_chart = genotyping_df[(genotyping_df['Clade'] == '2.3.4.4b') & (genotyping_df['Genotype'] != '-')].shape[0] > 0
    
    # Consistent color mapping across all seasons
    global_color_map = {}
    
    # H Subtype colors
    h_labels = sorted([h for h in genotyping_df['H_Subtype'].dropna().unique() if h != '-'])
    h_color_map = {}
    for h in h_labels:
        h_color_map[h] = okabe_ito_colors[len(h_color_map) % len(okabe_ito_colors)] if is_colorblind else subtype_base_colors.get(str(h), '#888888')
    global_color_map['h_subtypes'] = h_color_map
    
    # Clade colors for each H subtype
    for h in valid_h_subtypes:
        base_c = subtype_base_colors.get(str(h), '#888888')
        h_clades = sorted([c for c in genotyping_df[genotyping_df['H_Subtype'] == h]['Root_Clade'].dropna().unique() if c not in ["Unassigned", "No dataset available", "-"]])
        
        if is_colorblind:
            clade_palette = okabe_ito_colors
        else:
            clade_palette = generate_shades(base_c, len(h_clades) + 2)
        
        clade_map = {}
        for idx, clade in enumerate(h_clades):
            clade_map[clade] = clade_palette[idx % len(clade_palette)]
        clade_map["Unassigned"] = clade_palette[-2] if len(clade_palette) > 1 else '#888888'
        clade_map["Others"] = clade_palette[-1]
        global_color_map[f'clades_{h}'] = clade_map
    
    # Genotype colors for 2.3.4.4b
    if include_genotype_chart:
        g_labels = sorted([g for g in genotyping_df[genotyping_df['Clade'] == '2.3.4.4b']['Genotype'].dropna().unique() if g not in ["-", "None", "", "nan", "Unassigned"]])
        genotype_palette = okabe_ito_colors if is_colorblind else colors.qualitative.Vivid
        genotype_map = {}
        for idx, genotype in enumerate(g_labels):
            genotype_map[genotype] = genotype_palette[idx % len(genotype_palette)]
        global_color_map['genotypes'] = genotype_map
            
    total_charts = 1 + len(valid_h_subtypes) + (1 if include_genotype_chart else 0)

    # CLADE EVOLUTION BAR CHARTS
    if 'WEEK' in genotyping_df.columns and not genotyping_df['WEEK'].isna().all():
        df_all = genotyping_df.dropna(subset=['WEEK']).copy()
        
        subplot_titles_evo = ["<b>H Subtype Weekly Evolution</b>"] + [f"<b>Clade Evolution for {h}</b>" for h in valid_h_subtypes]
        if include_genotype_chart:
            subplot_titles_evo.append("<b>Genotype Evolution (Clade 2.3.4.4b)</b>")
        
        fig_evo = make_subplots(rows=total_charts, cols=1, shared_xaxes=True, vertical_spacing=0.08, subplot_titles=subplot_titles_evo)
        # Registry to track which subtype/clade/genotype corresponds to which trace for visibility toggling
        trace_registry = []

        def add_stacked_bars(df_subset, group_col, row_num, detail_col=None, palette_type='clade', base_color='#888888', color_key=''):
            if df_subset.empty:
                return
                
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
            
            if not hover_details.empty:
                group_counts = group_counts.merge(hover_details, on=['WEEK', group_col], how='left')
            else:
                group_counts['Hover_Details'] = float('nan')
                
            group_counts['Hover_Details'] = group_counts['Hover_Details'].fillna("")
            group_counts['Pct'] = (group_counts['Count'] / group_counts['TotalWeek']) * 100

            def make_hover_extra(row):
                val = str(row[group_col])
                details = row['Hover_Details']
                if not details:
                    return ""
                
                if group_col == 'Final_Label' and val.endswith("-like"):
                    return "<br><br><b>Clade Breakdown:</b><br>" + details
                elif group_col == 'Genotype':
                    return "<br><br><b>Sub-genotypes Breakdown:</b><br>" + details
                return ""

            group_counts['Hover_Extra'] = group_counts.apply(make_hover_extra, axis=1)

            unique_groups = group_counts[group_col].unique()
            
            # Get color mapping from global_color_map
            if color_key:
                color_map = global_color_map.get(color_key, {})
            else:
                color_map = {}

            for idx, g in enumerate(unique_groups):
                g_df = group_counts[group_counts[group_col] == g]
                
                if color_key:
                    color = color_map.get(str(g), '#888888')
                elif palette_type == 'h_subtype':
                    color = okabe_ito_colors[idx % len(okabe_ito_colors)] if is_colorblind else subtype_base_colors.get(str(g), '#888888')
                else:
                    color = '#888888'
                
                display_group_col = "Group" if group_col == "Final_Label" else group_col.replace("_", " ")
                
                fig_evo.add_trace(
                    go.Bar(
                        x=g_df['WEEK'], 
                        y=g_df['Pct'],
                        name=str(g),
                        text=g_df['Count'],
                        textposition='inside',
                        insidetextanchor='middle',
                        textangle=0,
                        textfont=dict(color=get_contrast_text_color(color), size=14, family='Arial'),
                        marker_color=color,
                        legendgroup=f"row_{row_num}",
                        legendgrouptitle_text=subplot_titles_evo[row_num-1].replace("<b>","").replace("</b>",""),
                        customdata=g_df[['Count', 'TotalWeek', 'Pct', 'Hover_Extra']],
                        hovertemplate=(
                            "<b>Week:</b> %{x|%V, %Y}<br>"
                            "<b>" + display_group_col + ":</b> " + str(g) + "<br>"
                            "<b>Occurrence:</b> %{customdata[0]}/%{customdata[1]} samples (%{customdata[2]:.1f}%)"
                            "%{customdata[3]}<extra></extra>"
                        )
                    ), 
                    row=row_num, col=1
                )
                trace_registry.append({'col': group_col, 'val': str(g)})

            fig_evo.update_yaxes(title_text="Frequency (%)", range=[0, 100], row=row_num, col=1)

        add_stacked_bars(df_all, 'H_Subtype', 1, palette_type='h_subtype', color_key='h_subtypes')

        current_row = 2
        for h in valid_h_subtypes:
            sub_df = df_all[df_all['H_Subtype'] == h].copy()
            sub_df['Final_Label'] = sub_df['Root_Clade']
            add_stacked_bars(sub_df, 'Final_Label', current_row, detail_col='Clade', palette_type='clade', base_color=subtype_base_colors.get(str(h), '#888888'), color_key=f'clades_{h}')
            current_row += 1

        if include_genotype_chart:
            sub_df = df_all[(df_all['Clade'] == '2.3.4.4b') & (df_all['Genotype'] != '-')].copy()
            add_stacked_bars(sub_df, 'Genotype', current_row, detail_col='Sub-genotype', palette_type='genotype', color_key='genotypes')

        all_time_start = df_all['WEEK'].min() - pd.Timedelta(days=7)
        all_time_end = df_all['WEEK'].max() + pd.Timedelta(days=14)
        all_time_range = [all_time_start.strftime('%Y-%m-%d'), all_time_end.strftime('%Y-%m-%d')] if pd.notnull(all_time_start) else None

        season_ranges = {}
        for season_val in sorted(df_all['Season'].dropna().unique()):
            if season_val == "Unknown Season" or season_val == "All Time":
                continue
            try:
                y1 = int(season_val.split('-')[0])
                y2 = int(season_val.split('-')[1])
                s_start = pd.to_datetime(f'{y1}-W40-1', format='%G-W%V-%u')
                s_end = pd.to_datetime(f'{y2}-W39-7', format='%G-W%V-%u')
                season_ranges[season_val] = [s_start.strftime('%Y-%m-%d'), s_end.strftime('%Y-%m-%d')]
            except Exception:
                s_data = df_all[df_all['Season'] == season_val]
                if not s_data.empty:
                    s_start = s_data['WEEK'].min()
                    s_end = s_data['WEEK'].max()
                    season_ranges[season_val] = [s_start.strftime('%Y-%m-%d'), s_end.strftime('%Y-%m-%d')]

        dropdown_buttons_evo = []
        if all_time_range:
            all_time_args = {f"xaxis{i+1 if i>0 else ''}.range": all_time_range for i in range(total_charts)}
            # For "All Time", we want all traces to be visible regardless of subtype/clade/genotype, so we set all to True
            all_time_mask = [True] * len(trace_registry)
            dropdown_buttons_evo.append(dict(args=[{"visible": all_time_mask}, all_time_args], label="All Time", method="update"))
            
        for season_val, s_range in season_ranges.items():
            season_args = {f"xaxis{i+1 if i>0 else ''}.range": s_range for i in range(total_charts)}
            # Determine which traces should be visible for this season based on the active subtypes/clades/genotypes in the data for that season
            s_data = df_all[df_all['Season'] == season_val]
            active_groups = {
                'H_Subtype': set(s_data['H_Subtype'].dropna().astype(str).unique()),
                'Final_Label': set(s_data['Root_Clade'].dropna().astype(str).unique()) if 'Root_Clade' in s_data.columns else set(),
                'Genotype': set(s_data['Genotype'].dropna().astype(str).unique()) if 'Genotype' in s_data.columns else set()
            }
            
            # Build a True/False mask for every drawn trace by checking if its specific value exists in our active labels
            season_visibility_mask = [
                (trace['col'] in active_groups and trace['val'] in active_groups[trace['col']])
                for trace in trace_registry
            ]
            
            dropdown_buttons_evo.append(dict(args=[{"visible": season_visibility_mask}, season_args], label=f"Season {season_val}", method="update"))

        if "${params.protocol}" == "AVIAN":
            center_x = 0.45
        else:
            center_x = 0.465

        fig_evo.update_layout(
            barmode='stack',
            title=dict(text="<b>Evolution of Subtypes and Clades</b>", x=center_x, y=0.98, xanchor="center", yanchor="top", font=dict(size=24)),
            updatemenus=[dict(active=0, buttons=dropdown_buttons_evo, x=0.5, xanchor="center", y=1.07, yanchor="top", direction="down", showactive=True)],
            height=450 * total_charts,
            hovermode="closest",
            margin=dict(t=160, b=80, l=80, r=200),
            legend=dict(tracegroupgap=30, y=1, yanchor="top")
        )
        
        fig_evo.update_xaxes(tickformat="Week %V<br>%Y", showticklabels=True)
        fig_evo.write_html("CladeEvolutionReport.html")

    # Set up a dynamic grid for the subplots (1 column wide)
    cols = 1
    rows = total_charts
    
    # Define the subplot type as 'domain' for pie charts
    specs = [[{"type": "domain"}] for _ in range(rows)]
    
    # Subplot titles in bold
    subplot_titles = ["<b>H Subtype Distribution</b>"] + [f"<b>Clade Distribution for {h}</b>" for h in valid_h_subtypes]
    if include_genotype_chart:
        subplot_titles.append("<b>Genotype Distribution (Clade 2.3.4.4b)</b>")
    
    # Create the subplot figure with extra vertical spacing to allow titles to shift up
    fig = make_subplots(rows=rows, cols=cols, specs=specs, subplot_titles=subplot_titles, vertical_spacing=0.15)

    traces_per_season = total_charts
    total_traces = len(seasons) * traces_per_season

    for s_idx, season in enumerate(seasons):
        season_df = genotyping_df[genotyping_df['Season'] == season]
        is_visible = (s_idx == 0)

        # H Subtype pie chart
        h_counts = season_df['H_Subtype'].value_counts().reset_index()
        h_counts.columns = ['Label', 'Count']
        
        if h_counts.empty:
            fig.add_trace(go.Pie(labels=["No Data"], values=[1], name="H Subtypes", textinfo='none', hoverinfo='none', visible=is_visible, marker=dict(colors=['#f0f0f0'])), row=1, col=1)
        else:
            # Calculate the total to display as part of the customized text (e.g., "5/20")
            total_h = h_counts['Count'].sum()
            h_counts['Text'] = h_counts['Count'].astype(str) + '/' + str(total_h)
            
            h_color_map = global_color_map['h_subtypes']
            pie_colors = [h_color_map.get(str(lbl), '#888888') for lbl in h_counts['Label']]

            fig.add_trace(
                go.Pie(
                    labels=h_counts['Label'], 
                    values=h_counts['Count'], 
                    name="H Subtypes", 
                    text=h_counts['Text'],
                    texttemplate='<b>%{label}</b><br><b>%{text}</b><br><b>%{percent}</b>',
                    textposition='outside',
                    insidetextorientation='horizontal',
                    rotation=90, 
                    automargin=True,
                    hole=0.35,
                    marker=dict(colors=pie_colors, line=dict(color='#ffffff', width=2)),
                    hoverlabel=dict(font_size=14),
                    hovertemplate='<b>H Subtype:</b> %{label}<br><b>Count:</b> %{text}<br><b>Percentage:</b> %{percent}<extra></extra>',
                    visible=is_visible
                ),
                row=1, col=1
            )

        # Clade pie chart for each valid specific H Subtype
        for i, h in enumerate(valid_h_subtypes):
            # Calculate the proper row placement (adding 2 because row 1 is the overall chart)
            r = i + 2
            
            # Filter the dataframe for the specific subtype and get clade counts
            sub_df = season_df[season_df['H_Subtype'] == h]
            total_c = len(sub_df)
            
            if total_c == 0:
                fig.add_trace(go.Pie(labels=["No Data"], values=[1], name=str(h), textinfo='none', hoverinfo='none', visible=is_visible, marker=dict(colors=['#f0f0f0'])), row=r, col=1)
                continue

            orig_counts = sub_df.groupby(['Clade', 'Root_Clade']).size().reset_index(name='Count')
            orig_counts['Orig_Pct'] = orig_counts['Count'] / total_c
            orig_counts['Hover_Detail'] = orig_counts.apply(lambda x: f"- {x['Clade']}: {x['Count']}/{total_c} ({x['Orig_Pct']:.1%})" if x['Count'] > 0 else "", axis=1)
            
            root_grouped = orig_counts.groupby('Root_Clade').agg(Root_Count=('Count', 'sum'), Root_Hover_Details=('Hover_Detail', lambda x: "<br>".join([detail for detail in x if detail])), Num_Clades=('Clade', 'nunique')).reset_index()
            root_grouped['Root_Pct'] = root_grouped['Root_Count'] / total_c
            
            root_grouped['Final_Label'] = root_grouped.apply(
                lambda x: "Others" if x['Root_Pct'] < 0.02 and str(x['Root_Clade']) not in ["Unassigned", "No dataset available"] else x['Root_Clade'], 
                axis=1
            )
            
            final_grouped = root_grouped.groupby('Final_Label').agg(Final_Count=('Root_Count', 'sum'), Final_Hover_Details=('Root_Hover_Details', lambda x: "<br>".join([detail for detail in x if detail])), Total_Unique_Clades=('Num_Clades', 'sum')).reset_index()
            
            final_grouped['Hover_Extra'] = final_grouped.apply(lambda x: "<br><br><b>Clade Breakdown:</b><br>" + x['Final_Hover_Details'] if str(x['Final_Label']).endswith("-like") or str(x['Final_Label']) == "Others" else "", axis=1)
            
            # Calculate the total clades for this specific subtype for the customized text
            final_grouped['Text'] = final_grouped['Final_Count'].astype(str) + '/' + str(total_c)
            
            clade_color_map = global_color_map[f'clades_{h}']
            pie_colors = [clade_color_map.get(str(lbl), '#888888') for lbl in final_grouped['Final_Label']]

            fig.add_trace(
                go.Pie(
                    labels=final_grouped['Final_Label'], 
                    values=final_grouped['Final_Count'], 
                    name=str(h), 
                    text=final_grouped['Text'], 
                    texttemplate='<b>%{label}</b><br><b>%{text}</b><br><b>%{percent}</b>',
                    textposition='outside',
                    rotation=90,
                    automargin=True,
                    insidetextorientation='horizontal',
                    hole=0.35,
                    marker=dict(colors=pie_colors, line=dict(color='#ffffff', width=2)),
                    hoverlabel=dict(font_size=14, align='left'),
                    customdata=final_grouped['Hover_Extra'],
                    hovertemplate='<b>Group:</b> %{label}<br><b>Total Count:</b> %{text}<br><b>Group Percentage:</b> %{percent}%{customdata}<extra></extra>',
                    visible=is_visible
                ),
                row=r, col=1
            )
            
        # Genotype chart for Clade 2.3.4.4b
        if include_genotype_chart:
            r = total_charts
            sub_df = season_df[(season_df['Clade'] == '2.3.4.4b') & (season_df['Genotype'] != '-')]
            total_g = len(sub_df)
            
            if total_g == 0:
                fig.add_trace(go.Pie(labels=["No Data"], values=[1], name="Genotypes", textinfo='none', hoverinfo='none', visible=is_visible, marker=dict(colors=['#f0f0f0'])), row=r, col=1)
            else:
                orig_counts = sub_df.groupby(['Genotype', 'Sub-genotype']).size().reset_index(name='Count')
                orig_counts['Orig_Pct'] = orig_counts['Count'] / total_g
                orig_counts['Hover_Detail'] = orig_counts.apply(
                    lambda x: f"- {x['Sub-genotype']}: {x['Count']}/{total_g} ({x['Orig_Pct']:.1%})" if str(x['Sub-genotype']) not in ["-", "None", "", "nan", "Unassigned"] else "", axis=1
                )

                root_grouped = orig_counts.groupby('Genotype').agg(
                    Final_Count=('Count', 'sum'),
                    Hover_Details=('Hover_Detail', lambda x: "<br>".join([d for d in x if d]))
                ).reset_index()

                root_grouped['Hover_Extra'] = root_grouped.apply(
                    lambda x: "<br><br><b>Sub-genotypes Breakdown:</b><br>" + x['Hover_Details'] if x['Hover_Details'] else "", axis=1
                )
                root_grouped['Text'] = root_grouped['Final_Count'].astype(str) + '/' + str(total_g)

                genotype_color_map = global_color_map['genotypes']
                pie_colors = [genotype_color_map.get(str(lbl), '#888888') for lbl in root_grouped['Genotype']]

                fig.add_trace(
                    go.Pie(
                        labels=root_grouped['Genotype'],
                        values=root_grouped['Final_Count'],
                        name="Genotypes 2.3.4.4b",
                        text=root_grouped['Text'],
                        texttemplate='<b>%{label}</b><br><b>%{text}</b><br><b>%{percent}</b>',
                        textposition='outside',
                        rotation=90,
                        automargin=True,
                        insidetextorientation='horizontal',
                        hole=0.35,
                        marker=dict(colors=pie_colors, line=dict(color='#ffffff', width=2)),
                        hoverlabel=dict(font_size=14, align='left'),
                        customdata=root_grouped['Hover_Extra'],
                        hovertemplate='<b>Genotype:</b> %{label}<br><b>Total Count:</b> %{text}<br><b>Percentage:</b> %{percent}%{customdata}<extra></extra>',
                        visible=is_visible
                    ),
                    row=r, col=1
                )

    dropdown_buttons = []
    for s_idx, season in enumerate(seasons):
        visibility_array = [False] * total_traces
        start_idx = s_idx * traces_per_season
        end_idx = start_idx + traces_per_season
        for j in range(start_idx, end_idx):
            visibility_array[j] = True
            
        # Prevent the word 'Season' from attaching to 'All Time'
        display_label = season if season == "All Time" else f"Season {season}"
            
        dropdown_buttons.append(dict(
            args=[
                {"visible": visibility_array}, 
                {"title": dict(text=f"<b>Subtype and Clade Report - {display_label}</b>", x=0.5, y=0.98, xanchor="center", yanchor="top")}
            ], 
            label=display_label, 
            method="update"
        ))
        
    first_lbl = seasons[0] if seasons[0] == "All Time" else f"Season {seasons[0]}"
    default_title = f"<b>Subtype and Clade Report - {first_lbl}</b>"

    # Push all subplot titles upward by adjusting their y coordinate
    for annotation in fig['layout']['annotations']:
        annotation['y'] += 0.02

    # Final layout adjustments
    fig.update_layout(
        title=dict(text=default_title, x=0.5, y=0.98, xanchor="center", yanchor="top", font=dict(size=24)),
        updatemenus=[dict(active=0, buttons=dropdown_buttons, x=0.5, xanchor="center", y=1.07, yanchor="top", direction="down", showactive=True)],
        height=750 * rows,
        showlegend=False,
        hovermode="closest",
        margin=dict(t=220, b=80, l=120, r=120),
        uniformtext=dict(minsize=10, mode='show') 
    )

    fig.write_html("CladeGraphicReport.html")
    """
}