#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MutationsGraphicReport {
    errorStrategy 'ignore'
    debug true

    input:
    path(relevant_mutations)

    output:
    path("MutationsReport.html"), emit: report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Load the mutations dataset, only first sheet (All_Proteins)
    df = pd.read_excel("${relevant_mutations}", keep_default_na=False)
    lengths_df = pd.read_csv("${params.protocols[params.protocol].resources}/annotations.csv")
    lengths_dict = dict(zip(lengths_df['Protein'].astype(str), lengths_df['Length']))

    # Standardize missing values and stringify key columns
    df['MARKER_ID'] = df['MARKER_ID'].replace('', 'None').fillna('None').astype(str)
    df['EFFECT'] = df['EFFECT'].replace('', 'Unknown').fillna('Unknown').astype(str)
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

    color_map = {
        'Marker': '#fe0000',      
        'Substitution': '#2243f5',
        'Deletion': '#000000',    
        'Insertion': '#00ff73'    
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
    group_cols = ['Plot_Group', 'POSITION', 'POSITION_REF', 'AA_MUTATION', 'MARKER_ID', 'EFFECT', 'Color_Category', 'ColorCode', 'FOUND_IN']
    
    # Function to list unique items in a column
    def list_unique_items(data_column):
        valid_items = []
        for item in data_column.unique():
            # Only add the item if it is not a blank space
            if str(item).strip() != '':
                valid_items.append(str(item))
        # Join everything nicely with a comma
        return ', '.join(valid_items)

    # Group the data and calculate metrics
    df_grouped = df.groupby(group_cols, dropna=False).agg(
        Sample_Count=('SAMPLE_ID', 'nunique'),
        Sample_IDs=('SAMPLE_ID', list_unique_items),
        Subtypes=('SUBTYPE', list_unique_items)
    ).reset_index()

    # Bring in the total sample numbers to calculate percentages
    df_grouped = pd.merge(df_grouped, total_samples_per_group, on='Plot_Group', how='left')
    
    # Calculate the percentage in two simple mathematical steps
    df_grouped['Percentage'] = (df_grouped['Sample_Count'] / df_grouped['Total_Group_Samples']) * 100
    df_grouped['Percentage'] = df_grouped['Percentage'].round(2)

    # Prepare subplots
    groups = sorted(df_grouped['Plot_Group'].unique())
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
            
            # Pack aggregated data including Total_Group_Samples at index 9
            hover_data = mut_df[['Sample_IDs', 'Subtypes', 'AA_MUTATION', 'EFFECT', 'MARKER_ID', 'Sample_Count', 'Percentage', 'FOUND_IN', 'POSITION_REF', 'Total_Group_Samples']].values
            
            fig.add_trace(
                go.Scatter(
                    x=mut_df['POSITION'],
                    y=mut_df['Percentage'],
                    mode='markers',
                    name=mut_type,
                    marker=dict(
                        color=mut_df['ColorCode'].tolist(), 
                        size=12, 
                        line=dict(width=1, color='DarkSlateGrey')
                    ),
                    customdata=hover_data,
                    hovertemplate=(
                        "<b>Position:</b> %{x}<br>"
                        "<b>Reference Position(H5N1 numbering):</b> %{customdata[8]}<br>"
                        "<b>Mutation:</b> %{customdata[2]}<br>"
                        "<b>Effect:</b> %{customdata[3]} [Found In: %{customdata[7]}]<br>"
                        "<b>Occurrence:</b> %{customdata[5]}/%{customdata[9]} sample(s) (%{customdata[6]}%)<br>"
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
            fig.update_xaxes(range=[0, max_length], title_text="Position", row=i, col=1)
        else:
            fig.update_xaxes(title_text="Position", row=i, col=1)
            
        # Update Y-axis to standard linear percentage scale 0-100
        fig.update_yaxes(range=[0, 105], title_text="Frequency (%)", row=i, col=1)

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
        margin=dict(t=100, b=80, l=80, r=80)
    )

    fig.write_html("MutationsReport.html")
    """
}