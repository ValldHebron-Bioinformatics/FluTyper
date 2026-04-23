#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsGraphicReport {
    errorStrategy 'ignore'
    debug true

    input:
    path(filtered_mutations)

    output:
    path("MutationsReport.html"), emit: report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Load the mutations dataset, only first sheet (All_Proteins)
    df = pd.read_excel("${filtered_mutations}", keep_default_na=False)
    lengths_df = pd.read_csv("${params.protocols[params.protocol].resources}/annotations.csv")
    lengths_dict = dict(zip(lengths_df['Protein'].astype(str), lengths_df['Length']))

    # Standardize missing values and replace pipes with a line break + spaces for indentation
    df['EFFECT'] = df['EFFECT'].replace('', 'Unknown').fillna('Unknown').astype(str).str.replace(' | ', '<br>                 ')
    df['SUBTYPE'] = df['SUBTYPE'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['REF_SUBTYPE'] = df['REF_SUBTYPE'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['FOUND_IN'] = df['FOUND_IN'].replace('', 'Unknown').fillna('Unknown').astype(str)
    df['POSITION_REF'] = df['POSITION_REF'].replace('', 'Unknown').fillna('Unknown').astype(str)

    # Define the coloring logic based on mutation type or marker status
    def get_mutation_category(row):
        if str(row.get('MARKER', 'No')) == 'Yes':
            return 'Marker'
        return str(row.get('MUTATION_TYPE', 'Unknown'))

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

    # Define plot groups, extracting from REF_SUBTYPE specifically for HA and NA
    def get_plot_group(row):
        protein = str(row.get('PROTEIN', 'Unknown'))
        if protein in ['HA1', 'HA2', 'NA']:
            subtype = str(row.get('REF_SUBTYPE', 'Unknown'))
            return f"{protein} - {subtype}"
        return protein

    df['Plot_Group'] = df.apply(get_plot_group, axis=1).astype(str)

    # Calculate total unique samples per protein(Plot_Group)
    total_samples_per_group = df.groupby('Plot_Group')['SAMPLE_ID'].nunique().reset_index()
    total_samples_per_group.rename(columns={'SAMPLE_ID': 'Total_Group_Samples'}, inplace=True)

    # Define the relevant columns for grouping and aggregation
    group_cols = ['Plot_Group', 'POSITION', 'POSITION_REF', 'AA_MUTATION', 'Color_Category', 'ColorCode']
    
    # Function to list unique items in a column
    def list_unique_items(data_column, joiner=', '):
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
        base_protein = group_name.split(' - ')[0]
        segment_num = segment_mapping.get(base_protein, 99)
        return (segment_num, group_name)

    # Prepare subplots using the custom segment order
    groups = sorted(df_grouped['Plot_Group'].unique(), key=custom_sort_key)
    rows_count = len(groups)
    
    # Calculate vertical spacing to avoid overlap
    row_height = 400
    vert_spacing = 80
    total_figure_height = row_height * rows_count
    
    if rows_count > 1:
        spacing = vert_spacing / total_figure_height
    else:
        spacing = 0

    fig = make_subplots(rows=rows_count, cols=1, subplot_titles=groups, vertical_spacing=spacing)

    # Global legend
    for mut_type, color in color_map.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode='markers',
                name=mut_type,
                marker=dict(color=color, size=12, line=dict(width=1, color='DarkSlateGrey')),
                legendgroup=mut_type,
                showlegend=True
            ),
            row=1, col=1
        )

    # Make scatter plots for each group
    for i, group in enumerate(groups, start=1):
        group_df = df_grouped[df_grouped['Plot_Group'] == group]
        
        for mut_type in group_df['Color_Category'].unique():
            mut_df = group_df[group_df['Color_Category'] == mut_type]
            
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
                        "<b>Reference Position (H5N1 numbering):</b> %{customdata[7]}<br>"
                        "<b>Mutation:</b> %{customdata[2]}<br>"
                        "<b>Effect(s):</b> %{customdata[3]}<br>"
                        "                 <b>Found in:</b>  %{customdata[6]}<br>"
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

    fig.update_layout(
        title_text=(
            "Mutation Summary per Protein<br>"
            "<span style='font-size:14px; color:gray; font-weight:normal;'>"
            "Displaying markers and mutations occurring at a frequency exceeding ${params.threshold * 100}% among samples per protein group."
            "</span>"
        ),
        height=total_figure_height,
        showlegend=True,
        hovermode="closest",
        hoverlabel=dict(align="left"), 
        margin=dict(t=100, b=80, l=80, r=80)
    )

    fig.write_html("MutationsReport.html")
    """
}