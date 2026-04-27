#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CladeGraphicReport {
    errorStrategy 'ignore'
    debug true

    input:
    path(genotyping_file)

    output:
    path("CladeGraphicReport.html"), emit: report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from plotly import colors

    # Okabe-Ito Colors Palette (colorblind-friendly)
    okabe_ito_colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999', '#000000']
    # Cris Colors Palette (she hates Barça colors)
    cris_colors = ['#F9DC5C', '#CD733D', '#C84630', '#94B0DA', '#676F86', '#3A2D32']
    
    if "${params.colorblind}".lower() == "true":
        color_dict= okabe_ito_colors
    else:
        color_dict = cris_colors
    # Dataframe preparation
    genotyping_df = pd.read_csv("${genotyping_file}")
    genotyping_df['H_Subtype'] = genotyping_df['Subtype'].astype(str).str[:2]

    unique_h_subtypes = sorted(genotyping_df['H_Subtype'].dropna().unique())
    
    # Filter out subtypes 
    valid_h_subtypes = []
    for h in unique_h_subtypes:
        clades_for_h = genotyping_df[genotyping_df['H_Subtype'] == h]['Clade'].unique()
        if not (len(clades_for_h) == 1 and clades_for_h[0] == "-"):
            valid_h_subtypes.append(h)
    genotyping_df['Clade'] = genotyping_df['Clade'].replace("unassigned", "Unassigned").replace("-", "No dataset available")
    
    total_charts = 1 + len(valid_h_subtypes)

    # Set up a dynamic grid for the subplots (1 column wide)
    cols = 1
    rows = total_charts
    
    # Define the subplot type as 'domain' for pie charts
    specs = [[{"type": "domain"}] for _ in range(rows)]
    
    # Subplot titles in bold
    subplot_titles = ["<b>H Subtype Distribution</b>"] + [f"<b>Clade Distribution for {h}</b>" for h in valid_h_subtypes]
    
    # Create the subplot figure
    fig = make_subplots(rows=rows, cols=cols, specs=specs, subplot_titles=subplot_titles,vertical_spacing=0.1)

    # H Subtype pie chart
    h_counts = genotyping_df['H_Subtype'].value_counts().reset_index()
    h_counts.columns = ['Label', 'Count']
    
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
            textposition='auto',
            insidetextorientation='horizontal',
            automargin=True,
            hole=0.35,
            marker=dict(colors=color_dict, line=dict(color='#ffffff', width=2)),
            hoverlabel=dict(font_size=14),
            hovertemplate='<b>H Subtype:</b> %{label}<br><b>Count:</b> %{text}<br><b>Percentage:</b> %{percent}<extra></extra>'
        ),
        row=1, col=1
    )

    # Clade pie chart for each valid specific H Subtype
    for i, h in enumerate(valid_h_subtypes):
        # Calculate the proper row placement (adding 2 because row 1 is the overall chart)
        r = i + 2
        
        # Filter the dataframe for the specific subtype and get clade counts
        sub_df = genotyping_df[genotyping_df['H_Subtype'] == h]
        c_counts = sub_df['Clade'].value_counts().reset_index()
        c_counts.columns = ['Label', 'Count']
        
        # Calculate the total clades for this specific subtype for the customized text
        total_c = c_counts['Count'].sum()
        c_counts['Text'] = c_counts['Count'].astype(str) + '/' + str(total_c)
        
        fig.add_trace(
            go.Pie(
                labels=c_counts['Label'], 
                values=c_counts['Count'], 
                name=str(h), 
                text=c_counts['Text'], 
                texttemplate='<b>%{label}</b><br><b>%{text}</b><br><b>%{percent}</b>',
                textposition='auto',
                insidetextorientation='horizontal',
                automargin=True,
                hole=0.35,
                marker=dict(colors=color_dict, line=dict(color='#ffffff', width=2)),
                hoverlabel=dict(font_size=14),
                hovertemplate='<b>Clade:</b> %{label}<br><b>Count:</b> %{text}<br><b>Percentage:</b> %{percent}<extra></extra>'

            ),
            row=r, col=1
        )

    # Final layout adjustments
    fig.update_layout(
        title_text="<b>Subtype and Clade Report</b>",
        height=600 * rows,
        showlegend=False,
        hovermode="closest",
        margin=dict(t=80, b=80, l=80, r=80)
    )

    fig.write_html("CladeGraphicReport.html")
    """
}