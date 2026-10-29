#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsGraphicReport {
    errorStrategy 'ignore'
    debug true
    // This process generates an interactive HTML report visualizing the mutations detected in the samples.

    input:
    path(full_mutations)
    path(metadata_file)

    output:
    path("MutationsReport.html"), emit: report

    script:
    def meta_str = metadata_file ? metadata_file.toString() : ""
    """
    #!/usr/bin/env python3
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import json

    # Load the mutations dataset and guarantee clean string IDs immediately
    df = pd.read_excel("${full_mutations}", keep_default_na=False)
    if 'SAMPLE_ID' in df.columns:
        df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str).str.strip()

    lengths_df = pd.read_csv("${params.protocols[params.protocol].resources}/annotations.csv")
    lengths_dict = dict(zip(lengths_df['Protein'].astype(str), lengths_df['Length']))

    # Load metadata and prepare the Season, Age_Group, and Sex columns
    df_meta = pd.read_csv("${meta_str}", skipinitialspace=True) if "${meta_str}" else pd.DataFrame()
    
    if not df_meta.empty:
        # Strip invisible characters and standardize to uppercase
        df_meta.columns = [str(c).replace('\\ufeff', '').strip().upper() for c in df_meta.columns]
        
        # Standardize the ID column to guarantee the merge executes
        if 'ID' in df_meta.columns and 'SAMPLE_ID' not in df_meta.columns:
            df_meta = df_meta.rename(columns={'ID': 'SAMPLE_ID'})
            
        if 'SAMPLE_ID' in df_meta.columns:
            df_meta['SAMPLE_ID'] = df_meta['SAMPLE_ID'].astype(str).str.strip()
            df_meta = df_meta.dropna(subset=['SAMPLE_ID'])
            df_meta = df_meta.drop_duplicates(subset=['SAMPLE_ID'], keep='first')

            # Extract demographic columns safely
            if 'AGE GROUP' in df_meta.columns:
                df_meta['Age_Group'] = df_meta['AGE GROUP'].fillna('Sense dades').astype(str).str.strip()
            elif 'AGE_GROUP' in df_meta.columns:
                df_meta['Age_Group'] = df_meta['AGE_GROUP'].fillna('Sense dades').astype(str).str.strip()
            else:
                df_meta['Age_Group'] = 'Sense dades'

            if 'SEX' in df_meta.columns:
                df_meta['Sex'] = df_meta['SEX'].fillna('Sense dades').astype(str).str.strip()
            else:
                df_meta['Sex'] = 'Sense dades'

            # Identify columns to merge and safely include DATE if present
            merge_cols = ['SAMPLE_ID', 'Age_Group', 'Sex']
            if 'DATE' in df_meta.columns:
                merge_cols.append('DATE')
                
            df = pd.merge(df, df_meta[merge_cols], on='SAMPLE_ID', how='left')

    # Fill missing demographic columns so the rest of the script can always reference them
    if 'Age_Group' not in df.columns:
        df['Age_Group'] = 'Sense dades'
    if 'Sex' not in df.columns:
        df['Sex'] = 'Sense dades'
    df['Age_Group'] = df['Age_Group'].fillna('Sense dades')
    df['Sex']       = df['Sex'].fillna('Sense dades')

    # Calculate the ISO season using a robust function to prevent NaN float errors
    if 'DATE' in df.columns:
        def get_season(date_val):
            if pd.isna(date_val):
                return "Unknown Season"
            iso_week = date_val.isocalendar().week
            iso_year = date_val.isocalendar().year
            s_year = iso_year if iso_week >= 40 else iso_year - 1
            return f"{s_year}-{s_year + 1}"

        dates = pd.to_datetime(df['DATE'], errors='coerce')
        df['Season'] = dates.apply(get_season)
    else:
        df['Season'] = "Unknown Season"

    # Build the full cross-product of (season x age x sex) views so every
    # combination of filters has its own set of traces in the figure.
    df_all_time = df.copy()
    df_all_time['Season'] = 'All Time'
    df_expanded = pd.concat([df_all_time, df], ignore_index=True)
    df_expanded = df_expanded[df_expanded['Season'].notna()]

    # Standardize missing values and replace pipes with a line break + spaces for indentation
    df_expanded['EFFECT'] = df_expanded['EFFECT'].replace('', 'Unknown').fillna('Unknown').astype(str).str.replace(' | ', '<br>                 ')
    df_expanded['SUBTYPE'] = df_expanded['SUBTYPE'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df_expanded['REF_SUBTYPE'] = df_expanded['REF_SUBTYPE'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df_expanded['FOUND_IN'] = df_expanded['FOUND_IN'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df_expanded['POSITION_REF'] = df_expanded['POSITION_REF'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df_expanded['POSITION'] = pd.to_numeric(df_expanded['POSITION'], errors='coerce')

    # Define the coloring logic based on mutation type or marker status
    def get_mutation_category(row):
        '''
        Categorizes each mutation based on its type or marker status.
        '''
        if str(row.get('MARKER', 'No')) == 'Yes':
            return 'Marker'
        return str(row.get('MUTATION_TYPE', 'Unknown'))

    df_expanded['Color_Category'] = df_expanded.apply(get_mutation_category, axis=1)

    # Define color mapping 
    if "${params.colorblind}".lower() == "true":
        color_map = {'Marker': '#D55E00', 'Substitution': '#0072B2', 'Deletion': '#000000', 'Insertion': '#CC79A7'}
    else:
        color_map = {'Marker': '#C84630', 'Substitution': '#94B0DA', 'Deletion': '#3A2D32', 'Insertion': '#F9DC5C'}
    df_expanded['ColorCode'] = df_expanded['Color_Category'].map(lambda x: color_map.get(x, '#aaaaaa'))

    def get_plot_group(row):
        '''
        Determines the plot group for each mutation based on the protein and subtype.
        For HUMAN protocol, all proteins are strictly separated by subtype to avoid mixing.
        For AVIAN protocol, only HA and NA surface proteins are separated by subtype.
        '''
        protein = str(row.get('PROTEIN', 'Unknown'))
        subtype = str(row.get('REF_SUBTYPE', 'Unknown'))
        
        if "${params.protocol}" == "HUMAN":
            return f"{protein} - {subtype}"
        else:
            if protein in ['HA1', 'HA2', 'NA ']:
                return f"{protein} - {subtype}"
            return protein

    df_expanded['Plot_Group'] = df_expanded.apply(get_plot_group, axis=1).astype(str)

    # Canonical key linking each mutation to its Frequency-Evolution report.
    # MUST replicate get_plot_name()'s branching in DateGraphicReport.nf EXACTLY:
    # HUMAN protocol -> every protein gets a subtype suffix.
    # AVIAN protocol -> only HA1/HA2/NA (surface proteins) get a subtype suffix;
    #                   internal proteins (PB2, PB1, PA, NP, M1/M2, NS1/NS2) are
    #                   pooled across subtypes into one report, e.g. "PB2" not "PB2_H5N1".
    def make_evo_key(protein, subtype):
        clean_subtype = str(subtype).replace('/', '_').strip()
        if "${params.protocol}" == "HUMAN" or protein in ['HA1', 'HA2', 'NA ']:
            return f"{protein}_{clean_subtype}"
        return protein

    df_expanded['EVO_KEY'] = df_expanded.apply(
        lambda r: make_evo_key(str(r.get('PROTEIN', 'Unknown')), str(r.get('REF_SUBTYPE', 'Unknown'))),
        axis=1
    )

    # Build the complete list of filter values including Sense dades
    age_order = {'0-2': 0, '3-4': 1, '5-14': 2, '15-65': 3, '>65': 4}
    age_groups = ['All'] + sorted(
        [a for a in df_expanded['Age_Group'].unique() if str(a).strip() not in ['nan', '', 'None']],
        key=lambda x: age_order.get(str(x).strip(), 99)
    )
    sexs = ['All'] + sorted(
        [g for g in df_expanded['Sex'].unique() if str(g).strip() not in ['nan', '', 'None']]
    )

    # Calculate total unique samples based on Subtype rather than individual Plot Group
    def totals_for_view(age_val, sex_val):
        view = df_expanded.copy()
        if age_val != 'All':
            view = view[view['Age_Group'] == age_val]
        if sex_val != 'All':
            view = view[view['Sex'] == sex_val]
            
        t = view.groupby(['Plot_Group', 'Season'])['SAMPLE_ID'].nunique().reset_index()
        t.rename(columns={'SAMPLE_ID': 'Total_Group_Samples'}, inplace=True)
        
        t['Age_Filter'] = age_val
        t['Sex_Filter'] = sex_val
        return t[['Plot_Group', 'Season', 'Total_Group_Samples', 'Age_Filter', 'Sex_Filter']]

    all_totals = pd.concat(
        [totals_for_view(a, g) for a in age_groups for g in sexs],
        ignore_index=True
    )

    group_cols = ['Plot_Group', 'Season', 'POSITION', 'POSITION_REF', 'AA_MUTATION', 'Color_Category', 'ColorCode', 'EVO_KEY']
    
    def list_unique_items(data_column, joiner=', '):
        '''
        Returns a string of unique, non-empty items from a pandas Series.
        '''
        valid_items = []
        for item in data_column.unique():
            if str(item).strip() != '':
                valid_items.append(str(item))
        return joiner.join(valid_items)

    def grouped_for_view(age_val, sex_val):
        view = df_expanded.copy()
        if age_val != 'All':
            view = view[view['Age_Group'] == age_val]
        if sex_val != 'All':
            view = view[view['Sex'] == sex_val]

        g = view.groupby(group_cols, dropna=False).agg(
            Sample_Count=('SAMPLE_ID', 'nunique'),
            Sample_IDs=('SAMPLE_ID', list_unique_items),
            Subtypes=('SUBTYPE', list_unique_items),
            EFFECT=('EFFECT', lambda x: list_unique_items(x, '<br>                 ')),
            FOUND_IN=('FOUND_IN', list_unique_items)
        ).reset_index()

        totals = all_totals[
            (all_totals['Age_Filter'] == age_val) &
            (all_totals['Sex_Filter'] == sex_val)
        ][['Plot_Group', 'Season', 'Total_Group_Samples']]

        g = pd.merge(g, totals, on=['Plot_Group', 'Season'], how='left')
        g['Percentage'] = (g['Sample_Count'] / g['Total_Group_Samples'] * 100).round(2)
        g['Age_Filter']    = age_val
        g['Sex_Filter'] = sex_val
        return g

    df_grouped = pd.concat(
        [grouped_for_view(a, ge) for a in age_groups for ge in sexs],
        ignore_index=True
    )

    segment_mapping = {
        'PB2': 1, 'PB1': 2, 'PB1-F2': 2, 'PA': 3, 'PA-X': 3, 'HA1': 4,
        'HA2': 4, 'NP': 5, 'NA ': 6, 'M1': 7, 'M2': 7, 'NS1': 8, 'NS2': 8,
    }

    def custom_sort_key(group_name):
        '''
        Custom sorting key for plot groups based on biological segment order.
        '''
        base_protein = group_name.split(' - ')[0]
        return (segment_mapping.get(base_protein, 99), group_name)

    groups = sorted(df_grouped['Plot_Group'].unique(), key=custom_sort_key)
    rows_count = len(groups)
    
    row_height = 400
    vert_spacing = 80
    total_figure_height = max(row_height * rows_count, 600)
    spacing = vert_spacing / total_figure_height if rows_count > 1 else 0

    fig = make_subplots(rows=rows_count, cols=1, subplot_titles=groups, vertical_spacing=spacing)

    epitope_definitions = {
        'HA1 - A(H3N2)': [
            {'name': 'RBD', 'positions': [98, 152, 153, 154, 155, 156], 'color': '#E41A1C'},
            {'name': 'RBD-130LOOP', 'positions': [131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148], 'color': '#377EB8'},
            {'name': 'RBD-180LOOP', 'positions': [183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195], 'color': '#4DAF4A'},
            {'name': 'RBD-220LOOP', 'positions': [218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230], 'color': '#984EA3'},
            {'name': 'A site', 'positions': [122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168], 'color': '#FF7F00'},
            {'name': 'B site', 'positions': [128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198], 'color': '#F781BF'},
            {'name': 'C site', 'positions': [44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312], 'color': '#A65628'},
            {'name': 'D site', 'positions': [96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248], 'color': '#FFFF33'},
            {'name': 'E site', 'positions': [57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265], 'color': '#00CED1'},
        ],
        'HA1 - A(H1N1)pdm09': [
            {'name': 'Cb', 'positions': [70, 71, 72, 73, 74, 75], 'color': '#E41A1C'},
            {'name': 'Sa', 'positions': [124, 125, 153, 154, 155, 156, 157, 159, 160, 161, 162, 163, 164], 'color': '#377EB8'},
            {'name': 'RBD', 'positions': [91, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 143, 144, 145, 149, 150, 151, 152, 180, 181, 182, 183, 215, 216, 217, 218, 219, 220, 223, 224, 225, 226, 227], 'color': '#4DAF4A'},
            {'name': 'Ca2', 'positions': [137, 138, 139, 140, 141, 142, 221, 222], 'color': '#984EA3'},
            {'name': 'Ca1', 'positions': [166, 167, 168, 169, 170, 203, 204, 205, 235, 236, 237], 'color': '#FF7F00'},
            {'name': 'Sb', 'positions': [184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195], 'color': '#F781BF'},
        ]
    }

    available_seasons = df_expanded['Season'].unique()
    valid_seasons = [str(s) for s in available_seasons if str(s) not in ['All Time', 'Unknown Season', 'nan']]
    sorted_seasons = sorted(valid_seasons, reverse=True)
    if 'Unknown Season' in available_seasons:
        sorted_seasons.append('Unknown Season')
    sorted_seasons.append('All Time')
    default_season = sorted_seasons[0]

    for i, group in enumerate(groups, start=1):

        if "${params.protocol}" == "HUMAN" and group in epitope_definitions:
            for epitope in epitope_definitions[group]:
                if 'positions' in epitope:
                    x_vals     = epitope['positions']
                    y_vals     = [115] * len(x_vals)
                    width_vals = [1]   * len(x_vals)
                elif 'start' in epitope and 'end' in epitope:
                    width_val  = epitope['end'] - epitope['start']
                    x_vals     = [epitope['start'] + (width_val / 2)]
                    y_vals     = [115]
                    width_vals = [width_val]
                else:
                    continue

                fig.add_trace(
                    go.Bar(
                        x=x_vals, y=y_vals, width=width_vals,
                        marker=dict(color=epitope['color'], line=dict(width=0)),
                        opacity=0.4, hoverinfo='skip', showlegend=False,
                        meta="Any"
                    ),
                    row=i, col=1
                )

        group_df = df_grouped[df_grouped['Plot_Group'] == group]

        for season in available_seasons:
            for age_val in age_groups:
                for sex_val in sexs:
                    season_df = group_df[
                        (group_df['Season']        == season) &
                        (group_df['Age_Filter']    == age_val) &
                        (group_df['Sex_Filter'] == sex_val)
                    ]
                    if season_df.empty:
                        continue

                    for mut_type in season_df['Color_Category'].unique():
                        mut_df = season_df[season_df['Color_Category'] == mut_type]

                        hover_data = mut_df[['Sample_IDs','Subtypes','AA_MUTATION','EFFECT','Sample_Count','Percentage','FOUND_IN','POSITION_REF','Total_Group_Samples','EVO_KEY']].values

                        if mut_type == 'Marker':
                            scatter_mode   = 'markers+text'
                            scatter_text   = ["<b>" + str(x) + "</b>" for x in mut_df['AA_MUTATION']]
                            text_pos_array = ['top center' if idx % 2 == 0 else 'bottom center' for idx in range(len(mut_df))]
                        else:
                            scatter_mode   = 'markers'
                            scatter_text   = None
                            text_pos_array = None

                        if "${params.protocol}" == "AVIAN":
                            hover_template_str = (
                                "<b>Position:</b> %{x}<br>"
                                "<b>Reference Position (H5N1 numbering):</b> %{customdata[7]}<br>"
                                "<b>Mutation:</b> %{customdata[2]}<br>"
                                "<b>Effect(s):</b> %{customdata[3]}<br>"
                                "                 <b>Found in:</b>  %{customdata[6]}<br>"
                                "<b>Occurrence:</b> %{customdata[4]}/%{customdata[8]} sample(s) (%{customdata[5]}%)<br>"
                                "<extra></extra>"
                            )
                        else:
                            hover_template_str = (
                                "<b>Position:</b> %{x}<br>"
                                "<b>Mutation:</b> %{customdata[2]}<br>"
                                "<b>Effect(s):</b> %{customdata[3]}<br>"
                                "                 <b>Found in:</b>  %{customdata[6]}<br>"
                                "<b>Occurrence:</b> %{customdata[4]}/%{customdata[8]} sample(s) (%{customdata[5]}%)<br>"
                                "<extra></extra>"
                            )

                        trace_meta = json.dumps({
                            'season': str(season),
                            'age':    age_val,
                            'sex': sex_val
                        })
                        trace_visibility = (
                            str(season) == default_season and
                            age_val     == 'All' and
                            sex_val  == 'All'
                        )

                        fig.add_trace(
                            go.Scatter(
                                x=mut_df['POSITION'],
                                y=mut_df['Percentage'],
                                mode=scatter_mode,
                                text=scatter_text,
                                textposition=text_pos_array,
                                textfont=dict(size=11, color="black"),
                                name=mut_type,
                                marker=dict(
                                    color=mut_df['ColorCode'].tolist(),
                                    size=12,
                                    line=dict(width=1, color='DarkSlateGrey')
                                ),
                                customdata=hover_data,
                                hovertemplate=hover_template_str,
                                legendgroup=mut_type,
                                showlegend=False,
                                visible=trace_visibility,
                                meta=trace_meta
                            ),
                            row=i, col=1
                        )

        base_protein = group.split(' - ')[0]
        max_length = lengths_dict.get(group, lengths_dict.get(base_protein, None))
        
        if max_length:
            fig.update_xaxes(range=[0, max_length+5], title_text="Position", row=i, col=1)
        else:
            fig.update_xaxes(title_text="Position", row=i, col=1)
            
        fig.update_yaxes(range=[0, 115], title_text="Frequency (%)", row=i, col=1)
    
    fig.update_layout(
        barmode='overlay', 
        height=total_figure_height, 
        showlegend=False, 
        hovermode="closest",
        hoverlabel=dict(align="left"), 
        margin=dict(t=40, b=80, l=80, r=80),
        plot_bgcolor='#ececec', 
    )

    graph_html = fig.to_html(full_html=False, include_plotlyjs='cdn', div_id="plotly-graphs")
    default_val = ${params.threshold}*100

    legend_html = '<div style="display: flex; justify-content: center; flex-wrap: wrap; gap: 20px; margin-top: 15px; font-size: 14px; padding-top: 10px; border-top: 1px solid #eaeaea;">'
    for mut_type, color in color_map.items():
        legend_html += f'<div style="display: flex; align-items: center;"><span style="display: inline-block; width: 14px; height: 14px; background-color: {color}; border-radius: 50%; margin-right: 6px; border: 1px solid #555;"></span>{mut_type}</div>'
    legend_html += '</div>'

    if "${params.protocol}" == "HUMAN":
        legend_html += '<div style="display: flex; justify-content: center; gap: 40px; margin-top: 10px; padding-top: 10px; border-top: 1px dashed #eaeaea;">'
        
        legend_html += '<div style="display: flex; flex-direction: column; align-items: center;">'
        legend_html += '<div style="font-weight: bold; font-size: 13px; margin-bottom: 5px; color: #444;">Epítops A(H1N1)pdm09</div>'
        legend_html += '<div style="display: flex; flex-wrap: wrap; justify-content: center; gap: 10px; font-size: 11px; color: #333; max-width: 300px;">'
        for ep in epitope_definitions.get('HA1 - A(H1N1)pdm09', []):
            legend_html += f'<div style="display: flex; align-items: center;"><span style="display: inline-block; width: 12px; height: 12px; background-color: {ep["color"]}; opacity: 0.7; margin-right: 4px; border: 1px solid #999;"></span>{ep["name"]}</div>'
        legend_html += '</div></div>'
        
        legend_html += '<div style="display: flex; flex-direction: column; align-items: center;">'
        legend_html += '<div style="font-weight: bold; font-size: 13px; margin-bottom: 5px; color: #444;">Epítops A(H3N2)</div>'
        legend_html += '<div style="display: flex; flex-wrap: wrap; justify-content: center; gap: 10px; font-size: 11px; color: #333; max-width: 400px;">'
        for ep in epitope_definitions.get('HA1 - A(H3N2)', []):
            legend_html += f'<div style="display: flex; align-items: center;"><span style="display: inline-block; width: 12px; height: 12px; background-color: {ep["color"]}; opacity: 0.7; margin-right: 4px; border: 1px solid #999;"></span>{ep["name"]}</div>'
        legend_html += '</div></div>'
        
        legend_html += '</div>'

    current_protocol = "${params.protocol}".upper()
    if current_protocol == "AVIAN":
        subtitle_text = "Markers are always displayed."
        js_marker_bypass = "dataSeries.name === 'Marker'"
    else:
        subtitle_text = "Markers are filtered by the selected frequency threshold."
        js_marker_bypass = "false"

    # Season range selectors. FROM: 'All Time', then every dated season in chronological
    # order, plus 'Unknown Season' (single-season view only) when undated samples exist.
    # TO: every dated season in chronological order. Both start on the default season so the
    # initial view is the same single-season view as before.
    range_seasons = sorted(valid_seasons)
    from_seasons = ['All Time'] + range_seasons
    if 'Unknown Season' in available_seasons:
        from_seasons.append('Unknown Season')
    default_to = default_season if default_season in range_seasons else (range_seasons[-1] if range_seasons else '')
    season_from_options = "".join([f'<option value="{s}"{" selected" if s == default_season else ""}>{s}</option>' for s in from_seasons])
    season_to_options   = "".join([f'<option value="{s}"{" selected" if s == default_to else ""}>{s}</option>' for s in range_seasons])
    season_to_disabled  = "" if default_season in range_seasons else " disabled"
    age_options    = "".join([f'<option value="{a}">{a}</option>' for a in age_groups])
    sex_options = "".join([f'<option value="{g}">{g}</option>' for g in sexs])

    html_template = f'''
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Mutations Summary</title>
        <style>
            body {{ font-family: arial; text-align: center; margin: 0; padding: 0; }}
            .sticky-header {{
                position: sticky; top: 0; background-color: rgba(255, 255, 255, 0.96);
                padding: 15px 20px; z-index: 1000; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); border-bottom: 1px solid #eaeaea;
            }}
            .controls-container {{
                display: flex; justify-content: center; align-items: flex-end; gap: 30px;
                margin: 15px auto 5px auto; flex-wrap: wrap; max-width: 1000px;
            }}
            .control-group {{
                display: flex; flex-direction: column; align-items: center;
            }}
            .control-group label {{ font-size: 10px; font-weight: bold; color: #666; margin-bottom: 4px; }}
            .control-group select {{
                padding: 6px; border-radius: 4px; border: 1px solid #ccc;
                background: white; font-size: 14px; min-width: 130px;
            }}
            .control-group.slider-group {{ min-width: 250px; }}
            .graph-container {{ padding: 20px; }}
        </style>
    </head>
    <body>

        <div class="sticky-header">
            <h2 id="report-title" style="margin: 0 0 5px 0;">Mutation Summary per Protein - Season {default_season}</h2>
            <p style="color: gray; font-size: 14px; margin: 0;">{subtitle_text}</p>

            <div class="controls-container">
                <div class="control-group">
                    <label>SEASON FROM</label>
                    <select id="seasonSel" onchange="onSeasonFromChange()" style="min-width:160px;">
                        {season_from_options}
                    </select>
                </div>

                <div class="control-group">
                    <label>SEASON TO</label>
                    <select id="seasonToSel" onchange="onSeasonToChange()" style="min-width:160px;"{season_to_disabled}>
                        {season_to_options}
                    </select>
                </div>

                <div class="control-group">
                    <label>AGE GROUP</label>
                    <select id="ageSel" onchange="applyFilters()">
                        {age_options}
                    </select>
                </div>

                <div class="control-group">
                    <label>SEX</label>
                    <select id="sexSel" onchange="applyFilters()">
                        {sex_options}
                    </select>
                </div>

                <div class="control-group slider-group">
                    <label>MINIMUM FREQUENCY THRESHOLD: <span id="sliderValue">{default_val}%</span></label>
                    <input type="range" id="freqSlider" min="0" max="100" value="{default_val}"
                           oninput="applyFilters()" style="width: 100%; margin-top: 8px;">
                    <p style="color: gray; font-size: 11px; margin-top: 4px; margin-bottom: 0;">
                        Percentage relative to sequences in the selected season(s)
                    </p>
                </div>
            </div>
            
            {legend_html}
        </div>

        <div class="graph-container">
            {graph_html}
        </div>

        <script>
            // setInterval creates a loop that checks every 200ms if the Plotly graph has been rendered and contains data.
            var checkGraphReady = setInterval(function() {{
                var graphContainer = document.getElementById('plotly-graphs');
                
                if (graphContainer && graphContainer.data && graphContainer.data.length > 0) {{
                    clearInterval(checkGraphReady);

                    graphContainer.parsedMeta       = [];
                    graphContainer.originalYValues  = [];
                    for (var si = 0; si < graphContainer.data.length; si++) {{
                        var raw = graphContainer.data[si].meta;
                        try {{ graphContainer.parsedMeta.push(JSON.parse(raw)); }}
                        catch(e) {{ graphContainer.parsedMeta.push(raw); }}
                        var yArr = graphContainer.data[si].y;
                        graphContainer.originalYValues.push(yArr ? Array.from(yArr) : null);
                    }}

                    buildSeasonIndex(graphContainer);
                    applyFilters();

                    // Marker click -> ask the parent dashboard (index.html) to open the matching
                    // Frequency Evolution report with this mutation isolated.
                    graphContainer.on('plotly_click', function(evt) {{
                        var pt = evt.points && evt.points[0];
                        if (!pt || pt.data.name !== 'Marker' || !pt.customdata) return;

                        var mutation = pt.customdata[2];
                        var evoKey   = pt.customdata[9];
                        if (!evoKey) return;

                        window.parent.postMessage(
                            {{ type: 'openEvolution', evoKey: evoKey, mutation: mutation }},
                            '*'
                        );
                    }});
                }}
            }}, 200);

            // ---- Season range (SEASON FROM / SEASON TO) ----
            // Single-season and All Time views reuse the pre-built per-season traces unchanged.
            // A multi-season range reuses the 'All Time' traces as templates (they hold every
            // mutation point) and recomputes each point from the per-season traces:
            //   frequency = sum(per-season sample counts) / sum(per-season protein totals).
            // Seasons are disjoint sample sets, so summing per-season unique counts is exact.
            var EFFECT_JOIN = '<br>                 ';

            // Index of a season inside the chronological SEASON TO list (-1 for All Time / Unknown Season).
            function rangeIndexOf(season) {{
                var opts = document.getElementById('seasonToSel').options;
                for (var i = 0; i < opts.length; i++) {{
                    if (opts[i].value === season) return i;
                }}
                return -1;
            }}

            function onSeasonFromChange() {{
                var toSel   = document.getElementById('seasonToSel');
                var fromIdx = rangeIndexOf(document.getElementById('seasonSel').value);
                if (fromIdx < 0) {{
                    toSel.disabled = true;
                }} else {{
                    toSel.disabled = false;
                    if (toSel.selectedIndex < fromIdx) toSel.selectedIndex = fromIdx;
                }}
                applyFilters();
            }}

            function onSeasonToChange() {{
                var toSel   = document.getElementById('seasonToSel');
                var fromIdx = rangeIndexOf(document.getElementById('seasonSel').value);
                if (fromIdx >= 0 && toSel.selectedIndex < fromIdx) toSel.selectedIndex = fromIdx;
                applyFilters();
            }}

            // Chronological list of seasons when a multi-season range is active, otherwise null.
            function activeSeasonRange() {{
                var toSel   = document.getElementById('seasonToSel');
                var fromIdx = rangeIndexOf(document.getElementById('seasonSel').value);
                if (fromIdx < 0 || toSel.disabled || toSel.selectedIndex <= fromIdx) return null;
                var seasons = [];
                for (var i = fromIdx; i <= toSel.selectedIndex; i++) seasons.push(toSel.options[i].value);
                return seasons;
            }}

            function totalsKey(ds, meta, season) {{
                return [ds.xaxis || 'x', meta.age, meta.sex, season].join('||');
            }}

            function pointKey(ds, meta, pi) {{
                var cd = ds.customdata[pi];
                return [ds.xaxis || 'x', meta.age, meta.sex, ds.name, ds.x[pi], cd[2], cd[7], cd[9]].join('||');
            }}

            // One-time index of the per-season traces: protein totals and per-point rows.
            function buildSeasonIndex(gd) {{
                gd.originalCustom = [];
                gd.seasonTotals   = {{}};
                gd.seasonPoints   = {{}};
                gd.templateIdx    = [];
                gd.rangeCache     = {{}};
                gd.appliedCdKey   = 'original';
                for (var si = 0; si < gd.data.length; si++) {{
                    var ds   = gd.data[si];
                    var meta = gd.parsedMeta[si];
                    gd.originalCustom.push(ds.customdata ? ds.customdata.map(function(r) {{ return Array.from(r); }}) : null);
                    if (!meta || typeof meta !== 'object' || !ds.customdata || !ds.x) continue;
                    if (meta.season === 'All Time') {{
                        gd.templateIdx.push(si);
                        continue;
                    }}
                    for (var pi = 0; pi < ds.customdata.length; pi++) {{
                        var row = ds.customdata[pi];
                        gd.seasonTotals[totalsKey(ds, meta, meta.season)] = Number(row[8]);
                        var pk = pointKey(ds, meta, pi);
                        if (!gd.seasonPoints[pk]) gd.seasonPoints[pk] = {{}};
                        gd.seasonPoints[pk][meta.season] = row;
                    }}
                }}
            }}

            // Same rounding as pandas .round(2) (numpy: scale by 100, round half to even).
            function roundHalfEven2(value) {{
                var scaled  = value * 100;
                var rounded = Math.round(scaled);
                if (Math.abs(scaled % 1) === 0.5 && rounded % 2 !== 0) rounded -= 1;
                return rounded / 100;
            }}

            function addUnique(list, text, sep) {{
                if (text === null || text === undefined) return;
                String(text).split(sep).forEach(function(item) {{
                    var clean = item.trim();
                    if (clean !== '' && list.indexOf(clean) < 0) list.push(clean);
                }});
            }}

            // Aggregated rows ({{pct, cd}} or null per point) for every template trace of a range.
            function getRangeRows(gd, seasons) {{
                var cacheKey = seasons.join(',');
                if (gd.rangeCache[cacheKey]) return gd.rangeCache[cacheKey];
                var result = {{}};
                gd.templateIdx.forEach(function(si) {{
                    var ds   = gd.data[si];
                    var meta = gd.parsedMeta[si];
                    var orig = gd.originalCustom[si];
                    var rows = [];
                    var cds  = [];
                    var total = 0;
                    seasons.forEach(function(s) {{ total += gd.seasonTotals[totalsKey(ds, meta, s)] || 0; }});
                    for (var pi = 0; pi < orig.length; pi++) {{
                        var perSeason = gd.seasonPoints[pointKey(ds, meta, pi)] || {{}};
                        var count = 0, ids = [], subtypes = [], effects = [], foundIn = [];
                        seasons.forEach(function(s) {{
                            var r = perSeason[s];
                            if (!r) return;
                            count += Number(r[4]);
                            addUnique(ids, r[0], ',');
                            addUnique(subtypes, r[1], ',');
                            addUnique(effects, r[3], '<br>');
                            addUnique(foundIn, r[6], ',');
                        }});
                        if (count === 0 || total === 0) {{
                            rows.push(null);
                            cds.push(orig[pi]);
                            continue;
                        }}
                        var pct = roundHalfEven2(count / total * 100);
                        rows.push(pct);
                        cds.push([ids.join(', '), subtypes.join(', '), orig[pi][2], effects.join(EFFECT_JOIN),
                                  count, pct, foundIn.join(', '), orig[pi][7], total, orig[pi][9]]);
                    }}
                    result[si] = {{ pct: rows, cd: cds }};
                }});
                gd.rangeCache[cacheKey] = result;
                return result;
            }}

            function applyFilters() {{
                var minimumFrequency = parseFloat(document.getElementById('freqSlider').value);
                var seasonRange      = activeSeasonRange();
                // A range is drawn on the 'All Time' template traces with recomputed values.
                var activeSeason     = seasonRange ? 'All Time' : document.getElementById('seasonSel').value;
                var activeAge        = document.getElementById('ageSel').value;
                var activeSex     = document.getElementById('sexSel').value;

                document.getElementById('sliderValue').innerText = minimumFrequency + '%';

                var titleLabel = seasonRange
                    ? 'Seasons ' + seasonRange[0] + ' to ' + seasonRange[seasonRange.length - 1]
                    : (activeSeason === 'All Time' ? 'All Time' : 'Season ' + activeSeason);
                document.getElementById('report-title').innerText =
                    'Mutation Summary per Protein - ' + titleLabel;

                var graphContainer = document.getElementById('plotly-graphs');
                if (!graphContainer || !graphContainer.originalYValues) return;

                var rangeRows = seasonRange ? getRangeRows(graphContainer, seasonRange) : null;

                // Swap hover data of the template traces between original and range values.
                var cdKey = seasonRange ? seasonRange.join(',') : 'original';
                if (cdKey !== graphContainer.appliedCdKey && graphContainer.templateIdx.length > 0) {{
                    var cdValues = graphContainer.templateIdx.map(function(si) {{
                        return rangeRows ? rangeRows[si].cd : graphContainer.originalCustom[si];
                    }});
                    Plotly.restyle(graphContainer, {{ customdata: cdValues }}, graphContainer.templateIdx);
                    graphContainer.appliedCdKey = cdKey;
                }}

                var newY          = [];
                var newVisibility = [];

                for (var si = 0; si < graphContainer.data.length; si++) {{
                    var dataSeries  = graphContainer.data[si];
                    var meta        = graphContainer.parsedMeta[si];
                    var baselineY   = graphContainer.originalYValues[si];

                    if (!baselineY) {{
                        newY.push(null);
                        newVisibility.push(false);
                        continue;
                    }}

                    if (meta === 'Any') {{
                        newVisibility.push(true);
                        newY.push(baselineY);
                        continue;
                    }}

                    var seasonMatch = meta.season === activeSeason;
                    var ageMatch    = meta.age    === activeAge;
                    var sexMatch = meta.sex === activeSex;

                    if (!seasonMatch || !ageMatch || !sexMatch) {{
                        newVisibility.push(false);
                        newY.push(baselineY); 
                        continue;
                    }}

                    newVisibility.push(true);

                    if (rangeRows && rangeRows[si]) {{
                        var rangePct = rangeRows[si].pct;
                        var bypass   = {js_marker_bypass};
                        newY.push(rangePct.map(function(p) {{
                            if (p === null) return null;
                            return (bypass || p >= minimumFrequency) ? p : null;
                        }}));
                    }} else if ({js_marker_bypass}) {{
                        newY.push(baselineY);
                    }} else {{
                        var filteredY = [];
                        for (var pi = 0; pi < baselineY.length; pi++) {{
                            if (dataSeries.customdata && dataSeries.customdata[pi]) {{
                                var pct = dataSeries.customdata[pi][5];
                                filteredY.push(pct >= minimumFrequency ? baselineY[pi] : null);
                            }} else {{
                                filteredY.push(null);
                            }}
                        }}
                        newY.push(filteredY);
                    }}
                }}

                Plotly.restyle(graphContainer, {{ y: newY, visible: newVisibility }});
            }}
        </script>
    </body>
    </html>
    '''
    
    with open("MutationsReport.html", "w", encoding="utf-8") as f:
        f.write(html_template)
    """
}