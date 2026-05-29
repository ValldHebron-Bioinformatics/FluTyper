#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsGraphicReport {
    errorStrategy 'ignore'
    debug true

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

    # Load the mutations dataset, only first sheet (All_Proteins)
    df = pd.read_excel("${full_mutations}", keep_default_na=False)
    lengths_df = pd.read_csv("${params.protocols[params.protocol].resources}/annotations.csv")
    lengths_dict = dict(zip(lengths_df['Protein'].astype(str), lengths_df['Length']))

    # Load metadata and prepare the Season column
    df_meta = pd.read_csv("${meta_str}", skipinitialspace=True) if "${meta_str}" else pd.DataFrame()
    
    if not df_meta.empty:
        df_meta.columns = [str(c).strip().upper() for c in df_meta.columns]
        if 'ID' in df_meta.columns:
            # Merge mutations with metadata to get dates
            df = pd.merge(df, df_meta, left_on='SAMPLE_ID', right_on='ID', how='left')

    if 'DATE' in df.columns:
        df['DATE'] = pd.to_datetime(df['DATE'], errors='coerce')
        iso = df['DATE'].dt.isocalendar()
        s_year = iso.year.where(iso.week >= 40, iso.year - 1)
        df['Season'] = s_year.astype(str) + "-" + (s_year + 1).astype(str)
    else:
        df['Season'] = "Unknown Season"

    # Create an explicit 'All Time' dataset to guarantee global calculations remain accurate
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
        if str(row.get('MARKER', 'No')) == 'Yes':
            return 'Marker'
        return str(row.get('MUTATION_TYPE', 'Unknown'))

    df_expanded['Color_Category'] = df_expanded.apply(get_mutation_category, axis=1)

    if "${params.colorblind}".lower() == "true":
        # Colorblind-friendly palette
        color_map = {
            'Marker': '#D55E00',       
            'Substitution': '#0072B2', 
            'Deletion': '#000000',     
            'Insertion': '#CC79A7'    
        }
    else:
        # Cris Colors Palette
        color_map = {
            'Marker': '#C84630',       
            'Substitution': '#94B0DA', 
            'Deletion': '#3A2D32',     
            'Insertion': '#F9DC5C'    
        }
    df_expanded['ColorCode'] = df_expanded['Color_Category'].map(lambda x: color_map.get(x, '#aaaaaa'))

    # Define plot groups based on the selected Nextflow protocol
    def get_plot_group(row):
        protein = str(row.get('PROTEIN', 'Unknown'))
        subtype = str(row.get('REF_SUBTYPE', 'Unknown'))
        
        if "${params.protocol}" == "HUMAN":
            # For HUMAN, strictly separate all proteins by subtype to avoid H1N1/H3N2 mixing
            return f"{protein} - {subtype}"
        else:
            # For AVIAN, keep original behavior: only separate HA and NA surface proteins
            if protein in ['HA1', 'HA2', 'NA']:
                return f"{protein} - {subtype}"
            return protein

    df_expanded['Plot_Group'] = df_expanded.apply(get_plot_group, axis=1).astype(str)

    # Calculate total unique samples per protein AND per season
    total_samples_per_group = df_expanded.groupby(['Plot_Group', 'Season'])['SAMPLE_ID'].nunique().reset_index()
    total_samples_per_group.rename(columns={'SAMPLE_ID': 'Total_Group_Samples'}, inplace=True)

    # Define the relevant columns for grouping and aggregation
    group_cols = ['Plot_Group', 'Season', 'POSITION', 'POSITION_REF', 'AA_MUTATION', 'Color_Category', 'ColorCode']
    
    # Function to list unique items in a column
    def list_unique_items(data_column, joiner=', '):
        valid_items = []
        for item in data_column.unique():
            if str(item).strip() != '':
                valid_items.append(str(item))
        return joiner.join(valid_items)

    # Group the data and calculate metrics
    df_grouped = df_expanded.groupby(group_cols, dropna=False).agg(
        Sample_Count=('SAMPLE_ID', 'nunique'),
        Sample_IDs=('SAMPLE_ID', list_unique_items),
        Subtypes=('SUBTYPE', list_unique_items),
        EFFECT=('EFFECT', lambda x: list_unique_items(x, '<br>                 ')),
        FOUND_IN=('FOUND_IN', list_unique_items)
    ).reset_index()

    # Bring in the total sample numbers to calculate percentages
    df_grouped = pd.merge(df_grouped, total_samples_per_group, on=['Plot_Group', 'Season'], how='left')    
    # Calculate the percentage
    df_grouped['Percentage'] = (df_grouped['Sample_Count'] / df_grouped['Total_Group_Samples']) * 100
    df_grouped['Percentage'] = df_grouped['Percentage'].round(2)
    
    # Define biological segment mapping for sorting
    segment_mapping = {
        'PB2': 1,
        'PB1': 2,
        'PB1-F2': 2,
        'PA': 3,
        'PA-X': 3,
        'HA1': 4,
        'HA2': 4,
        'NP': 5,
        'NA': 6,
        'M1': 7,
        'M2': 7,
        'NS1': 8,
        'NS2': 8,
    }

    def custom_sort_key(group_name):
        base_protein = group_name.split(' - ')[0]
        segment_num = segment_mapping.get(base_protein, 99)
        return (segment_num, group_name)

    # Prepare subplots using the custom segment order
    groups = sorted(df_grouped['Plot_Group'].unique(), key=custom_sort_key)
    rows_count = len(groups)
    
    # Calculate vertical spacing to avoid overlap
    row_height = 400
    vert_spacing = 80
    total_figure_height = max(row_height * rows_count, 600)
    spacing = vert_spacing / total_figure_height if rows_count > 1 else 0

    fig = make_subplots(rows=rows_count, cols=1, subplot_titles=groups, vertical_spacing=spacing)

    epitope_definitions = {
        'HA1 - H3N2': [
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
        'HA1 - H1N1': [
            {'name': 'Cb', 'positions': [70, 71, 72, 73, 74, 75], 'color': '#E41A1C'},
            {'name': 'Sa', 'positions': [124, 125, 153, 154, 155, 156, 157, 159, 160, 161, 162, 163, 164], 'color': '#377EB8'},
            {'name': 'RBD', 'positions': [91, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 143, 144, 145, 149, 150, 151, 152, 180, 181, 182, 183, 215, 216, 217, 218, 219, 220, 223, 224, 225, 226, 227], 'color': '#4DAF4A'},
            {'name': 'Ca2', 'positions': [137, 138, 139, 140, 141, 142, 221, 222], 'color': '#984EA3'},
            {'name': 'Ca1', 'positions': [166, 167, 168, 169, 170, 203, 204, 205, 235, 236, 237], 'color': '#FF7F00'},
            {'name': 'Sb', 'positions': [184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195], 'color': '#F781BF'},
        ]
    }

    # Determine unique valid seasons and sort to ensure the newest season is default
    available_seasons = df_expanded['Season'].unique()
    valid_seasons = [str(s) for s in available_seasons if str(s) not in ['All Time', 'Unknown Season', 'nan']]
    sorted_seasons = sorted(valid_seasons, reverse=True)
    sorted_seasons.append('All Time')
    
    default_season = sorted_seasons[0]

    # Make scatter plots for each group, iterating through available seasons
    for i, group in enumerate(groups, start=1):
        
        # Add epitope backgrounds as native Bar charts to avoid background z-index hiding and JS failures
        if "${params.protocol}" == "HUMAN" and group in epitope_definitions:
            for epitope in epitope_definitions[group]:
                if 'positions' in epitope:
                    x_vals = epitope['positions']
                    y_vals = [115] * len(x_vals)
                    width_vals = [1] * len(x_vals)
                elif 'start' in epitope and 'end' in epitope:
                    width_val = epitope['end'] - epitope['start']
                    x_vals = [epitope['start'] + (width_val / 2)]
                    y_vals = [115]
                    width_vals = [width_val]
                else:
                    continue
                
                fig.add_trace(
                    go.Bar(
                        x=x_vals,
                        y=y_vals,
                        width=width_vals,
                        marker=dict(color=epitope['color'], line=dict(width=0)),
                        opacity=0.4,
                        hoverinfo='skip',
                        showlegend=False,
                        meta="Any"
                    ),
                    row=i, col=1
                )

        group_df = df_grouped[df_grouped['Plot_Group'] == group]
        
        for season in available_seasons:
            season_df = group_df[group_df['Season'] == season]
            if season_df.empty: continue
            
            for mut_type in season_df['Color_Category'].unique():
                mut_df = season_df[season_df['Color_Category'] == mut_type]
                
                # Pack aggregated data
                hover_data = mut_df[['Sample_IDs', 'Subtypes', 'AA_MUTATION', 'EFFECT', 'Sample_Count', 'Percentage', 'FOUND_IN', 'POSITION_REF', 'Total_Group_Samples']].values
                
                # Set mode and text for Markers only
                if mut_type == 'Marker':
                    scatter_mode = 'markers+text'
                    scatter_text = ["<b>" + str(x) + "</b>" for x in mut_df['AA_MUTATION']]
                    
                    # Create an alternating array to stagger text up and down
                    text_pos_array = ['top center' if idx % 2 == 0 else 'bottom center' for idx in range(len(mut_df))]
                else:
                    scatter_mode = 'markers'
                    scatter_text = None
                    text_pos_array = None

                # Construeix la targeta emergent dinàmicament en funció del protocol
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

                # Traces are initially hidden unless they belong to the default season
                trace_visibility = True if str(season) == default_season else False

                fig.add_trace(
                    go.Scatter(
                        x=mut_df['POSITION'],
                        y=mut_df['Percentage'],
                        mode=scatter_mode,
                        text=scatter_text,
                        textposition=text_pos_array,
                        textfont=dict(size=11, color="black"),
                        name=mut_type,
                        marker=dict(color=mut_df['ColorCode'].tolist(), size=12, line=dict(width=1, color='DarkSlateGrey')),
                        customdata=hover_data,
                        hovertemplate=hover_template_str,
                        legendgroup=mut_type,
                        showlegend=False,
                        visible=trace_visibility,
                        meta=season 
                    ),
                    row=i, col=1
                )

        # Force the X-axis range based on the reference lengths file
        base_protein = group.split(' - ')[0]
        max_length = lengths_dict.get(group, lengths_dict.get(base_protein, None))
        
        if max_length:
            fig.update_xaxes(range=[0, max_length+5], title_text="Position", row=i, col=1)
        else:
            fig.update_xaxes(title_text="Position", row=i, col=1)
            
        fig.update_yaxes(range=[0, 115], title_text="Frequency (%)", row=i, col=1)
    
    # Graph layout adjustments
    fig.update_layout(
        barmode='overlay', # Ensures that overlapping epitope regions seamlessly blend colors instead of shifting
        height=total_figure_height, 
        showlegend=False, 
        hovermode="closest",
        hoverlabel=dict(align="left"), 
        margin=dict(t=40, b=80, l=80, r=80),
        plot_bgcolor='#ececec', 
    )

    # Plotly graph to HTML
    graph_html = fig.to_html(full_html=False, include_plotlyjs='cdn', div_id="plotly-graphs")
    default_val = ${params.threshold}*100

    # Sticky legend in html
    legend_html = '<div style="display: flex; justify-content: center; flex-wrap: wrap; gap: 20px; margin-top: 15px; font-size: 14px; padding-top: 10px; border-top: 1px solid #eaeaea;">'
    for mut_type, color in color_map.items():
        legend_html += f'<div style="display: flex; align-items: center;"><span style="display: inline-block; width: 14px; height: 14px; background-color: {color}; border-radius: 50%; margin-right: 6px; border: 1px solid #555;"></span>{mut_type}</div>'
    legend_html += '</div>'

    # Split sticky legend for epitopes dynamically if HUMAN protocol is active
    if "${params.protocol}" == "HUMAN":
        legend_html += '<div style="display: flex; justify-content: center; gap: 40px; margin-top: 10px; padding-top: 10px; border-top: 1px dashed #eaeaea;">'
        
        # H1N1 Epitopes
        legend_html += '<div style="display: flex; flex-direction: column; align-items: center;">'
        legend_html += '<div style="font-weight: bold; font-size: 13px; margin-bottom: 5px; color: #444;">Epítops H1N1</div>'
        legend_html += '<div style="display: flex; flex-wrap: wrap; justify-content: center; gap: 10px; font-size: 11px; color: #333; max-width: 300px;">'
        for ep in epitope_definitions.get('HA1 - H1N1', []):
            legend_html += f'<div style="display: flex; align-items: center;"><span style="display: inline-block; width: 12px; height: 12px; background-color: {ep["color"]}; opacity: 0.7; margin-right: 4px; border: 1px solid #999;"></span>{ep["name"]}</div>'
        legend_html += '</div></div>'
        
        # H3N2 Epitopes
        legend_html += '<div style="display: flex; flex-direction: column; align-items: center;">'
        legend_html += '<div style="font-weight: bold; font-size: 13px; margin-bottom: 5px; color: #444;">Epítops H3N2</div>'
        legend_html += '<div style="display: flex; flex-wrap: wrap; justify-content: center; gap: 10px; font-size: 11px; color: #333; max-width: 400px;">'
        for ep in epitope_definitions.get('HA1 - H3N2', []):
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

    # Gather clean season names for the dropdown
    season_options = "".join([f'<option value="{s}">{s}</option>' for s in sorted_seasons])

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
                display: flex; justify-content: center; align-items: center; gap: 40px; margin: 15px auto 5px auto; max-width: 800px;
            }}
            .control-group {{
                display: flex; flex-direction: column; align-items: center; width: 100%;
            }}
            .graph-container {{ padding: 20px; }}
        </style>
    </head>
    <body>

        <div class="sticky-header">
            <h2 style="margin: 0 0 5px 0;">Mutation Summary per Protein</h2>
            <p style="color: gray; font-size: 14px; margin: 0;">{subtitle_text}</p>

            <div class="controls-container">
                <div class="control-group">
                    <label><b>Season Filter:</b></label>
                    <select id="seasonSel" onchange="applyFilters()" style="margin-top:8px; padding:5px; border-radius:4px; font-size:14px; width:100%; max-width: 250px;">
                        {season_options}
                    </select>
                </div>

                <div class="control-group">
                    <label><b>Minimum Frequency Threshold:</b> <span id="sliderValue">{default_val}%</span></label>
                    <input type="range" id="freqSlider" min="0" max="100" value="{default_val}" oninput="applyFilters()" style="width: 100%; margin-top:8px;">
                    <p style="color: gray; font-size: 11px; margin-top: 4px; margin-bottom: 0;">
                        Percentage relative to sequences in the active season
                    </p>
                </div>
            </div>
            
            {legend_html}
        </div>

        <div class="graph-container">
            {graph_html}
        </div>

        <script>
            var checkGraphReady = setInterval(function() {{
                var graphContainer = document.getElementById('plotly-graphs');
                
                // If the graph is loaded and has data, we can store the original Y values.
                if (graphContainer && graphContainer.data && graphContainer.data.length > 0) {{
                    clearInterval(checkGraphReady); 
                    
                    // Save the original Y values to a custom property so we don't lose them when filtering
                    graphContainer.originalYValues = [];
                    for (var seriesIndex = 0; seriesIndex < graphContainer.data.length; seriesIndex++) {{
                        var dataSeries = graphContainer.data[seriesIndex];
                        if (dataSeries.y) {{
                            graphContainer.originalYValues.push(Array.from(dataSeries.y));
                        }} else {{
                            graphContainer.originalYValues.push(null);
                        }}
                    }}
                    
                    applyFilters();
                }}
            }}, 200);

            function applyFilters() {{
                var minimumFrequency = document.getElementById('freqSlider').value;
                var activeSeason = document.getElementById('seasonSel').value;
                document.getElementById('sliderValue').innerText = minimumFrequency + '%';
                var graphContainer = document.getElementById('plotly-graphs');
                
                // If the graph hasn't been properly saved to memory yet, exit silently to avoid errors
                if (!graphContainer || !graphContainer.originalYValues) return;
                
                var newVerticalCoordinates = [];
                var newVisibility = [];
                
                for(var seriesIndex = 0; seriesIndex < graphContainer.data.length; seriesIndex++) {{
                    var dataSeries = graphContainer.data[seriesIndex];
                    var baselineYValues = graphContainer.originalYValues[seriesIndex];
                    
                    if (!baselineYValues) {{
                        newVerticalCoordinates.push(null);
                        newVisibility.push(false);
                        continue;
                    }}

                    // First filter: Check if the trace belongs to the selected season or is an epitope shape ("Any")
                    if (dataSeries.meta !== activeSeason && dataSeries.meta !== "Any") {{
                        newVisibility.push(false);
                        newVerticalCoordinates.push(baselineYValues); // Data hidden, Y structure maintained
                        continue;
                    }} else {{
                        newVisibility.push(true);
                    }}

                    // Second filter: Apply frequency threshold ignoring "Any" structural objects
                    if (dataSeries.meta === 'Any' || {js_marker_bypass}) {{
                        newVerticalCoordinates.push(baselineYValues);
                    }} else {{
                        var filteredYValues = [];
                        for(var pointIndex = 0; pointIndex < baselineYValues.length; pointIndex++) {{
                            // Safely check that the extra custom data exists before trying to read the percentage
                            if (dataSeries.customdata && dataSeries.customdata[pointIndex]) {{
                                var pointPercentage = dataSeries.customdata[pointIndex][5];
                                if (pointPercentage >= minimumFrequency) {{
                                    filteredYValues.push(baselineYValues[pointIndex]);
                                }} else {{
                                    filteredYValues.push(null);
                                }}
                            }} else {{
                                filteredYValues.push(null);
                            }}
                        }}
                        newVerticalCoordinates.push(filteredYValues);
                    }}
                }}
                
                // Package restyle updates for a single smooth execution
                Plotly.restyle(graphContainer, {{y: newVerticalCoordinates, visible: newVisibility}});
            }}
        </script>
    </body>
    </html>
    '''

    with open("MutationsReport.html", "w", encoding="utf-8") as f:
        f.write(html_template)
    """
}