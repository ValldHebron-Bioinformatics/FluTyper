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

    # Define text color (black or white) based on background color for readability
    def get_contrast_text_color(hex_str):
        hex_str = hex_str.lstrip('#')
        if len(hex_str) == 6:
            r, g, b = tuple(int(hex_str[i:i+2], 16) for i in (0, 2, 4))
            luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
            return '#000000' if luminance > 0.55 else '#ffffff'
        return '#ffffff'

    # Okabe-Ito Colors Palette (colorblind-friendly)
    okabe_ito_colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999', '#000000']
    # Cris Colors Palette (she hates Barça colors)
    cris_colors = ['#F9DC5C', '#CD733D', '#C84630', '#94B0DA', '#676F86', '#3A2D32', '#F2B5D4', "#20AA71", "#376AAD"]
    
    if "${params.colorblind}".lower() == "true":
        color_dict= okabe_ito_colors
    else:
        color_dict = cris_colors
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
        iso_cal = genotyping_df['DATE'].dt.isocalendar()
        s_year = iso_cal.year.where(iso_cal.week >= 40, iso_cal.year - 1)
        genotyping_df['Season'] = s_year.astype(str) + "-" + (s_year + 1).astype(str)
        genotyping_df['Season'] = genotyping_df['Season'].fillna("Unknown Season")
    else:
        genotyping_df['Season'] = "All Time"

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
            
    total_charts = 1 + len(valid_h_subtypes) + (1 if include_genotype_chart else 0)

    # CLADE EVOLUTION BAR CHARTS
    if 'WEEK' in genotyping_df.columns and not genotyping_df['WEEK'].isna().all():
        df_all = genotyping_df.dropna(subset=['WEEK']).copy()
        
        subplot_titles_evo = ["<b>H Subtype Weekly Evolution</b>"] + [f"<b>Clade Evolution for {h}</b>" for h in valid_h_subtypes]
        if include_genotype_chart:
            subplot_titles_evo.append("<b>Genotype Evolution (Clade 2.3.4.4b)</b>")
        
        fig_evo = make_subplots(rows=total_charts, cols=1, shared_xaxes=True, vertical_spacing=0.08, subplot_titles=subplot_titles_evo)

        def add_stacked_bars(df_subset, group_col, row_num, detail_col=None):
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
                
                if group_col == 'Final_Label' and (val.endswith("-like") or val == "Others"):
                    return "<br><br><b>Clade Breakdown:</b><br>" + details
                elif group_col == 'Genotype':
                    return "<br><br><b>Sub-genotypes Breakdown:</b><br>" + details
                return ""

            group_counts['Hover_Extra'] = group_counts.apply(make_hover_extra, axis=1)

            unique_groups = group_counts[group_col].unique()
            for idx, g in enumerate(unique_groups):
                g_df = group_counts[group_counts[group_col] == g]
                color = color_dict[idx % len(color_dict)]
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
            fig_evo.update_yaxes(title_text="Frequency (%)", range=[0, 100], row=row_num, col=1)

        add_stacked_bars(df_all, 'H_Subtype', 1)

        current_row = 2
        for h in valid_h_subtypes:
            sub_df = df_all[df_all['H_Subtype'] == h].copy()
            
            orig_counts = sub_df.groupby(['Clade', 'Root_Clade']).size().reset_index(name='Count')
            total_c = len(sub_df)
            orig_counts['Orig_Pct'] = orig_counts['Count'] / total_c if total_c > 0 else 0
            root_grouped_evo = orig_counts.groupby('Root_Clade').agg(Root_Count=('Count', 'sum')).reset_index()
            root_grouped_evo['Root_Pct'] = root_grouped_evo['Root_Count'] / total_c if total_c > 0 else 0
            
            sub_df['Final_Label'] = sub_df['Root_Clade'].apply(
                lambda x: "Others" if x in root_grouped_evo[root_grouped_evo['Root_Pct'] < 0.02]['Root_Clade'].values and x not in ["Unassigned", "No dataset available"] else x
            )
            
            add_stacked_bars(sub_df, 'Final_Label', current_row, detail_col='Clade')
            current_row += 1

        if include_genotype_chart:
            sub_df = df_all[(df_all['Clade'] == '2.3.4.4b') & (df_all['Genotype'] != '-')].copy()
            add_stacked_bars(sub_df, 'Genotype', current_row, detail_col='Sub-genotype')

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
                # REMOVED: pd.Timedelta(days=7) padding so it zooms strictly from Week 40 to Week 39
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
            dropdown_buttons_evo.append(dict(args=[all_time_args], label="All Time", method="relayout"))
            
        for season_val, s_range in season_ranges.items():
            season_args = {f"xaxis{i+1 if i>0 else ''}.range": s_range for i in range(total_charts)}
            dropdown_buttons_evo.append(dict(args=[season_args], label=f"Season {season_val}", method="relayout"))

        fig_evo.update_layout(
            barmode='stack',
            title=dict(text="<b>Evolution of Subtypes and Clades</b>", x=0.45, y=0.98, xanchor="center", yanchor="top", font=dict(size=24)),
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
                    marker=dict(colors=color_dict, line=dict(color='#ffffff', width=2)),
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
                    marker=dict(colors=color_dict, line=dict(color='#ffffff', width=2)),
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
                        marker=dict(colors=color_dict, line=dict(color='#ffffff', width=2)),
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
            
        dropdown_buttons.append(dict(args=[{"visible": visibility_array}, {"title": dict(text=f"<b>Subtype and Clade Report - Season {season}</b>", x=0.5, y=0.98, xanchor="center", yanchor="top")}], label=f"Season {season}", method="update"))
    default_title = f"<b>Subtype and Clade Report - Season {seasons[0]}</b>" if seasons else "<b>Subtype and Clade Report</b>"

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