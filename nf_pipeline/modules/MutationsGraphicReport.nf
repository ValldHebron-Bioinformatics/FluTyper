#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsGraphicReport {
    // This process generates a comprehensive HTML report visualizing mutation data using Plotly.
    // It reads mutation data from an Excel file, processes it, and creates interactive scatter 
    // plots for each protein group.
    errorStrategy 'ignore'
    debug true

    input:
    path(full_mutations)

    output:
    path("MutationsReport.html"), emit: report

    script:
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

    # Standardize missing values and replace pipes with a line break + spaces for indentation
    df['EFFECT'] = df['EFFECT'].replace('', 'Unknown').fillna('Unknown').astype(str).str.replace(' | ', '<br>                 ')
    df['SUBTYPE'] = df['SUBTYPE'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['REF_SUBTYPE'] = df['REF_SUBTYPE'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['FOUND_IN'] = df['FOUND_IN'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['POSITION_REF'] = df['POSITION_REF'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['POSITION'] = pd.to_numeric(df['POSITION'], errors='coerce')

    # Define the coloring logic based on mutation type or marker status
    def get_mutation_category(row):
        '''
        Categorizes each mutation based on its type or marker status.
        '''
        if str(row.get('MARKER', 'No')) == 'Yes':
            return 'Marker'
        return str(row.get('MUTATION_TYPE', 'Unknown'))

    df['Color_Category'] = df.apply(get_mutation_category, axis=1)

    # Define color mapping 
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
    df['ColorCode'] = df['Color_Category'].map(lambda x: color_map.get(x, '#aaaaaa'))

    # Define plot groups based on the selected Nextflow protocol
    def get_plot_group(row):
        '''
        Determines the plot group for each mutation based on the protein and subtype.
        For HUMAN protocol, all proteins are strictly separated by subtype to avoid mixing.
        For AVIAN protocol, only HA and NA surface proteins are separated by subtype.
        '''
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

    df['Plot_Group'] = df.apply(get_plot_group, axis=1).astype(str)

    # Calculate total unique samples per protein(Plot_Group)
    total_samples_per_group = df.groupby('Plot_Group')['SAMPLE_ID'].nunique().reset_index()
    total_samples_per_group.rename(columns={'SAMPLE_ID': 'Total_Group_Samples'}, inplace=True)

    # Define the relevant columns for grouping and aggregation
    group_cols = ['Plot_Group', 'POSITION', 'POSITION_REF', 'AA_MUTATION', 'Color_Category', 'ColorCode']
    

    def list_unique_items(data_column, joiner=', '):
        '''
        Returns a string of unique, non-empty items from a pandas Series.
        '''
        valid_items = []
        for item in data_column.unique():
            if str(item).strip() != '':
                valid_items.append(str(item))
        return joiner.join(valid_items)

    # Group the data and calculate metrics
    df_grouped = df.groupby(group_cols, dropna=False).agg(
        Sample_Count=('SAMPLE_ID', 'nunique'),
        Sample_IDs=('SAMPLE_ID', list_unique_items),
        Subtypes=('SUBTYPE', list_unique_items),
        EFFECT=('EFFECT', lambda x: list_unique_items(x, '<br>                 ')),
        FOUND_IN=('FOUND_IN', list_unique_items)
    ).reset_index()

    # Bring in the total sample numbers to calculate percentages
    df_grouped = pd.merge(df_grouped, total_samples_per_group, on='Plot_Group', how='left')
    
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
        '''
        Custom sorting key for plot groups based on biological segment order.
        '''
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
    
    if rows_count > 1:
        spacing = vert_spacing / total_figure_height
    else:
        spacing = 0

    fig = make_subplots(rows=rows_count, cols=1, subplot_titles=groups, vertical_spacing=spacing)

    # Epitope definitions for the HUMAN protocol
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

    # Make scatter plots for each group
    for i, group in enumerate(groups, start=1):
        
        # Add epitope background bars for HUMAN protocol before rendering the actual mutations
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
        
        for mut_type in group_df['Color_Category'].unique():
            mut_df = group_df[group_df['Color_Category'] == mut_type]
            
            # Pack aggregated data
            hover_data = mut_df[['Sample_IDs', 'Subtypes', 'AA_MUTATION', 'EFFECT', 'Sample_Count', 'Percentage', 'FOUND_IN', 'POSITION_REF', 'Total_Group_Samples']].values
            
            # Set mode and text for Markers only to add labels above/below the points
            if mut_type == 'Marker':
                scatter_mode = 'markers+text'
                scatter_text = ["<b>" + str(x) + "</b>" for x in mut_df['AA_MUTATION']]
                # Create an alternating array to stagger text up and down
                text_pos_array = ['top center' if idx % 2 == 0 else 'bottom center' for idx in range(len(mut_df))]
            else:
                scatter_mode = 'markers'
                scatter_text = None
                text_pos_array = None

            if "${params.protocol}" == "AVIAN":
                ref = "H5N1"
            elif "${params.protocol}" == "HUMAN":
                ref = "H1N1"
            else:
                ref = "Unknown"

            # Plotly scatter trace for the mutation type
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
                    hovertemplate=(
                        "<b>Position:</b> %{x}<br>"
                        "<b>Reference Position (" + ref + " numbering):</b> %{customdata[7]}<br>"
                        "<b>Mutation:</b> %{customdata[2]}<br>"
                        "<b>Effect(s):</b> %{customdata[3]}<br>"
                        "                  <b>Found in:</b>  %{customdata[6]}<br>"
                        "<b>Occurrence:</b> %{customdata[4]}/%{customdata[8]} sample(s) (%{customdata[5]}%)<br>"
                        "<extra></extra>"
                    ),
                    legendgroup=mut_type,
                    showlegend=False
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
    
    # Graph layout adjustments (showlegend=False because we use an html sticky legend)
    fig.update_layout(
        height=total_figure_height, 
        showlegend=False, 
        hovermode="closest",
        hoverlabel=dict(align="left"), 
        margin=dict(t=40, b=80, l=80, r=80) 
    )

    # Plotly graph to HTML
    graph_html = fig.to_html(full_html=False, include_plotlyjs='cdn', div_id="plotly-graphs")
    default_val = ${params.threshold}*100

    # Sticky legend in html
    legend_html = '<div style="display: flex; justify-content: center; flex-wrap: wrap; gap: 20px; margin-top: 15px; font-size: 14px; padding-top: 10px; border-top: 1px solid #eaeaea;">'
    for mut_type, color in color_map.items():
        legend_html += f'<div style="display: flex; align-items: center;"><span style="display: inline-block; width: 14px; height: 14px; background-color: {color}; border-radius: 50%; margin-right: 6px; border: 1px solid #555;"></span>{mut_type}</div>'
    legend_html += '</div>'

    # Determine JavaScript filtering logic and subtitle text based on protocol
    current_protocol = "${params.protocol}".upper()
    if current_protocol == "AVIAN":
        subtitle_text = "Markers are always displayed."
        js_marker_bypass = "dataSeries.name === 'Marker'"
    else:
        subtitle_text = "Markers are filtered by the selected frequency threshold."
        js_marker_bypass = "false"

    # Create the full HTML template with embedded graph and slider
    html_template = f'''
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Mutations Summary</title>
        <style>
            /* CSS handles the visual presentation. The sticky-header keeps the controls visible when scrolling down a large graph. */
            body {{
                font-family: arial; 
                text-align: center; 
                margin: 0; 
                padding: 0;
            }}
            // This class makes the header with the slider and legend stick to the top
            .sticky-header {{ 
                position: sticky;
                top: 0;
                background-color: rgba(255, 255, 255, 0.96); 
                padding: 15px 20px;
                z-index: 1000; /* Ensures the header stays on top of the graph elements */
                box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); 
                border-bottom: 1px solid #eaeaea;
            }}
            // This class centers the slider
            .slider-container {{
                margin: 15px auto 5px auto; 
                max-width: 600px;
            }}
            .graph-container {{
                padding: 20px;
            }}
        </style>
    </head>
    <body>
        <div class="sticky-header">
            <h2 style="margin: 0 0 5px 0;">Mutation Summary per Protein</h2>
            <p style="color: gray; font-size: 14px; margin: 0;">{subtitle_text}</p>
            <div class="slider-container">
                <label><b>Minimum Frequency Threshold:</b> <span id="sliderValue">{default_val}%</span></label>
                <br><br>
                <input type="range" id="freqSlider" min="0" max="100" value="{default_val}" oninput="applyFrequencyFilter(this.value)" style="width: 80%;">
                
                <p style="color: gray; font-size: 12px; margin-top: 8px; margin-bottom: 0;">
                    Percentage of frequency of mutations in the sequences assessed
                </p>
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
                
                // If the graph is loaded and has data, we can store the original Y values.
                if (graphContainer && graphContainer.data && graphContainer.data.length > 0) {{
                    clearInterval(checkGraphReady); // Stop the checking loop once we confirm the graph exists
                    
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
                    
                    // Initial plot update to apply the default threshold upon loading
                    applyFrequencyFilter(document.getElementById('freqSlider').value);
                }}
            }}, 200);

            function applyFrequencyFilter(minimumFrequency) {{
                // Update the text label next to the slider to show the current number
                document.getElementById('sliderValue').innerText = minimumFrequency + '%';
                var graphContainer = document.getElementById('plotly-graphs');
                
                // If the graph hasn't been properly saved to memory yet, exit silently to avoid browser errors
                if (!graphContainer || !graphContainer.originalYValues) return;
                
                var newVerticalCoordinates = [];
                
                for(var seriesIndex = 0; seriesIndex < graphContainer.data.length; seriesIndex++) {{
                    var dataSeries = graphContainer.data[seriesIndex];
                    var baselineYValues = graphContainer.originalYValues[seriesIndex];
                    
                    if (!baselineYValues) {{
                        newVerticalCoordinates.push(null);
                        continue;
                    }}

                    // Apply conditional bypass logic based on the Nextflow protocol injection
                    if ({js_marker_bypass}) {{
                        // If bypass is true, just push the original data back without filtering
                        newVerticalCoordinates.push(baselineYValues);
                    }} else {{
                        var filteredYValues = [];
                        for(var pointIndex = 0; pointIndex < baselineYValues.length; pointIndex++) {{
                            // Safely check that the extra custom data (injected by Python) exists before trying to read the percentage
                            if (dataSeries.customdata && dataSeries.customdata[pointIndex]) {{
                                // Index [5] corresponds to the exact column where you stored the percentage in your Python dataframe logic
                                var pointPercentage = dataSeries.customdata[pointIndex][5];
                                if (pointPercentage >= minimumFrequency) {{
                                    filteredYValues.push(baselineYValues[pointIndex]); // Keep the point
                                }} else {{
                                    filteredYValues.push(null); // Setting it to null tells Plotly to visually hide this specific point
                                }}
                            }} else {{
                                filteredYValues.push(null);
                            }}
                        }}
                        newVerticalCoordinates.push(filteredYValues);
                    }}
                }}
                
                // Send all Y-coordinate updates to the graph at once for a smooth visual transition.
                Plotly.restyle(graphContainer, {{y: newVerticalCoordinates}});
            }}
        </script>
    </body>
    </html>
    '''

    with open("MutationsReport.html", "w", encoding="utf-8") as f:
        f.write(html_template)
    """
}
