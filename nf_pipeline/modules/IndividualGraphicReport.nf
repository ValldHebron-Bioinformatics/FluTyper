#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process IndividualGraphicReport {
    errorStrategy 'ignore'
    debug true

    input:
    tuple val(sample_id), path(individual_mutations)

    output:
    path("samples/${sample_id}/${sample_id}_MutationsReport.html"), emit: report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import os

    # Load the mutations dataset
    df = pd.read_csv("${individual_mutations}", keep_default_na=False)
    lengths_df = pd.read_csv("${params.protocols[params.protocol].resources}/annotations.csv")
    lengths_dict = dict(zip(lengths_df['Protein'].astype(str), lengths_df['Length']))

    # Standardize missing values and replace pipes with a line break + spaces for indentation
    df['EFFECT'] = df['EFFECT'].replace('', 'Unknown').fillna('Unknown').astype(str).str.replace(' | ', '<br>                 ')
    df['SUBTYPE'] = df['SUBTYPE'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['REF_SUBTYPE'] = df['REF_SUBTYPE'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['FOUND_IN'] = df['FOUND_IN'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['POSITION_REF'] = df['POSITION_REF'].replace('', 'Unknown').fillna('Unknown').astype(str)

    # Extract coordinates for ranges and single lines
    def extract_pos(pos_str):
        try:
            parts = str(pos_str).split('-')
            start = float(parts[0].strip())
            end = float(parts[-1].strip())
            return start, end
        except:
            return 0.0, 0.0

    df[['PLOT_START', 'PLOT_END']] = df['POSITION'].apply(lambda x: pd.Series(extract_pos(x)))
    df['PLOT_CENTER'] = (df['PLOT_START'] + df['PLOT_END']) / 2
    df['IS_RANGE']    = df['PLOT_END'] > df['PLOT_START']

    # Define categories
    def get_mutation_category(row):
        if str(row.get('MARKER', 'No')) == 'Yes':
            return 'Marker'
        return str(row.get('MUTATION_TYPE', 'Substitution'))

    df['Color_Category'] = df.apply(get_mutation_category, axis=1)

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

    def get_plot_group(row):
        protein = str(row.get('PROTEIN', 'Unknown'))
        if protein in ['HA1', 'HA2', 'NA']:
            subtype = str(row.get('REF_SUBTYPE', 'Unknown'))
            return f"{protein} - {subtype}"
        return protein

    df['Plot_Group'] = df.apply(get_plot_group, axis=1).astype(str)
    df = df.sort_values(['Plot_Group', 'PLOT_START']).reset_index(drop=True)

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
    groups = sorted(df['Plot_Group'].unique(), key=custom_sort_key)
    rows_count = len(groups)
    
    row_height = 180
    vert_spacing = 140
    # Not automatic to achieve better visibility
    total_figure_height = max(500, row_height * rows_count + 150)
    spacing = vert_spacing / total_figure_height if rows_count > 1 else 0

    fig = make_subplots(rows=rows_count, cols=1, subplot_titles=[f"<br><b>{group}</b><br> <br>" for group in groups], vertical_spacing=spacing)
    
    # Global legend
    for mut_type, color in color_map.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(symbol='square', color=color, size=14),
                name=mut_type, legendgroup=mut_type, showlegend=True
            ),
            row=1, col=1
        )

    # Make scatter plots for each group
    for i, group in enumerate(groups, start=1):
        group_df = df[df['Plot_Group'] == group]
        
        for mut_type in group_df['Color_Category'].unique():
            mut_df = group_df[group_df['Color_Category'] == mut_type].copy()
            color = color_map.get(mut_type, '#000000')
            
            # Draw rectangles for ranges and vertical lines for single positions
            xs_line, ys_line = [], []
            for _, r in mut_df.iterrows():
                if r['IS_RANGE']:
                    fig.add_shape(
                        type='rect', x0=r['PLOT_START'], x1=r['PLOT_END'],
                        y0=0, y1=1, fillcolor=color,
                        line_width=0, row=i, col=1
                    )
                else:
                    xs_line.extend([r['PLOT_CENTER'], r['PLOT_CENTER'], None])
                    ys_line.extend([0, 1, None])

            if xs_line:
                fig.add_trace(
                    go.Scatter(
                        x=xs_line, y=ys_line, mode='lines',
                        line=dict(color=color, width=2.0 if mut_type == 'Marker' else 1.5),
                        hoverinfo='skip', showlegend=False, legendgroup=mut_type
                    ), row=i, col=1
                )

            # Hover data for markers and mutations
            hover_data = mut_df[['SUBTYPE', 'AA_MUTATION', 'EFFECT', 'FOUND_IN', 'POSITION_REF', 'POSITION']].values
            # Set mode and text for Markers only
            if mut_type == 'Marker':
                scatter_mode = 'markers+text'
                scatter_text = ["<b>" + str(x) + "</b>" for x in mut_df['AA_MUTATION']]
                x_pos_array = ['top center' if idx % 2 == 0 else 'bottom center' for idx in range(len(mut_df))]
                y_pos_array = [1 if idx % 2 == 0 else 0 for idx in range(len(mut_df))]
            else:
                scatter_mode = 'markers'
                scatter_text = None
                x_pos_array = None
                y_pos_array = [0.5] * len(mut_df)
            
            if "${params.protocol}" == "AVIAN":
                ref = "H5N1"
            elif "${params.protocol}" == "HUMAN":
                ref = "H1N1"
            else:
                ref = "Unknown"

            fig.add_trace(
                go.Scatter(
                    x=mut_df['PLOT_CENTER'],
                    y=y_pos_array,
                    mode=scatter_mode,
                    text=scatter_text,
                    textposition=x_pos_array,
                    textfont=dict(size=11, color="black"),
                    name=mut_type,
                    marker=dict(color=color, size=10, opacity=0),
                    cliponaxis=False,
                    customdata=hover_data,
                    hovertemplate=(
                        "<b>Position:</b> %{customdata[5]}<br>"
                        "<b>Reference Position (" + ref + " numbering):</b> %{customdata[4]}<br>"
                        "<b>Mutation:</b> %{customdata[1]}<br>"
                        "<b>Effect(s):</b> %{customdata[2]}<br>"
                        "                 <b>Found in:</b>  %{customdata[3]}<br>"
                        "<extra></extra>"
                    ),
                    legendgroup=mut_type,
                    showlegend=False
                ), row=i, col=1
            )

        # Set x and y axes properties
        base_protein = group.split(' - ')[0]
        max_length = lengths_dict.get(group, lengths_dict.get(base_protein, 800))
        
        fig.update_xaxes(
                    range=[-1, max_length+10], 
                    title=dict(text="Amino acid position", standoff=20), 
                    ticks="outside", 
                    ticklen=15, 
                    tickcolor="rgba(0,0,0,0)", 
                    showgrid=False, 
                    zeroline=False, 
                    row=i, 
                    col=1
                )        
        fig.update_yaxes(range=[0, 1], showticklabels=False, showgrid=False, zeroline=False, row=i, col=1)

    fig.update_layout(
        title_text="<b>Genomic Barcode Profile - Sample ${sample_id}</b>",
        height=total_figure_height,
        showlegend=True,
        hovermode="closest",
        hoverlabel=dict(align="left"), 
        margin=dict(t=100, b=80, l=80, r=80),
        plot_bgcolor='#e2dfdf', 
    )

    os.makedirs("samples/${sample_id}", exist_ok=True)
    fig.write_html("samples/${sample_id}/${sample_id}_MutationsReport.html")
    """
}