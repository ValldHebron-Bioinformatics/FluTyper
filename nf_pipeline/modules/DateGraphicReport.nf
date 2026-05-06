#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DateGraphicReport {
    errorStrategy 'ignore'
    debug true
    input:
    path excel_file
    path meta_file

    output:
    path "FrequencyEvolution/**/*.html", emit: metadata

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import os

    def generate_plots(mut_file, meta_file):
        df_mut = pd.read_excel(mut_file, na_filter=False)
        df_meta = pd.read_csv(meta_file, skipinitialspace=True)

        df_meta['DATE'] = pd.to_datetime(df_meta['DATE'], format='%Y-%m-%d')
        df_meta['WEEK'] = df_meta['DATE'].dt.to_period('W').dt.to_timestamp()

        global_weeks = pd.DataFrame({'WEEK': sorted(df_meta['WEEK'].unique())})

        df_all = pd.merge(df_mut, df_meta[['ID', 'WEEK']], left_on='SAMPLE_ID', right_on='ID')
        df_all['REF_SUBTYPE'] = df_all.get('REF_SUBTYPE', 'Unknown').replace('', 'Unknown').fillna('Unknown').astype(str)
        
        # Updated logic to prevent mixing subtypes in HUMAN protocol timelines
        def get_plot_name(row):
            prot = str(row.get('PROTEIN', 'Unknown'))
            subtype = str(row.get('REF_SUBTYPE', 'Unknown')).replace('/', '_')
            
            if "${params.protocol}" == "HUMAN":
                return f"{prot}_{subtype}"
            else:
                if prot in ['HA1', 'HA2', 'NA']:
                    return f"{prot}_{subtype}"
                return prot
            
        df_all['PLOT_NAME'] = df_all.apply(get_plot_name, axis=1)
        
        df_markers = df_all[df_all['MUTATION_TYPE'] == 'Marker'].copy()

        hover_info = df_markers[['AA_MUTATION', 'POSITION_REF', 'EFFECT', 'FOUND_IN']].drop_duplicates(subset=['AA_MUTATION']).replace('', 'Unknown').fillna('Unknown')
        hover_info['EFFECT'] = hover_info['EFFECT'].astype(str).str.replace(' | ', '<br>                 ')

        prot_dict = {
            "HA":  ["HA1", "HA2"], "NA":  ["NA"], "PB2": ["PB2"],
            "PB1": ["PB1", "PB1-F2"], "PA":  ["PA", "PA-X"], "NP":  ["NP"],
            "MP":  ["M1", "M2"], "NS":  ["NS1", "NS2"]
        }

        if "${params.colorblind}".lower() == "true":
            color_list = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999', '#000000']
        else:
            color_list = ['#F9DC5C', '#CD733D', '#C84630', '#94B0DA', '#676F86', '#3A2D32']

        for segment, proteins in prot_dict.items():
            os.makedirs(f"FrequencyEvolution/{segment}", exist_ok=True)
            
            for protein in proteins:
                df_base_markers = df_markers[df_markers['PROTEIN'] == protein]
                if df_base_markers.empty: continue

                for plot_name in df_base_markers['PLOT_NAME'].unique():
                    df_prot = df_base_markers[df_base_markers['PLOT_NAME'] == plot_name]
                    
                    df_plot_all = df_all[df_all['PLOT_NAME'] == plot_name]
                    weekly_totals = df_plot_all.groupby('WEEK')['SAMPLE_ID'].nunique().reset_index(name='Samples (Week)')
                    
                    weekly_totals = pd.merge(global_weeks, weekly_totals, on='WEEK', how='left').fillna({'Samples (Week)': 0})
                    weekly_totals['Samples (Cumulative)'] = weekly_totals['Samples (Week)'].cumsum()
                    
                    weekly_mut_counts = df_prot.groupby(['WEEK', 'AA_MUTATION']).size().unstack(fill_value=0)
                    weekly_mut_counts = weekly_mut_counts.reindex(global_weeks['WEEK'], fill_value=0)
                    cum_mut_counts = weekly_mut_counts.cumsum()
                    
                    cum_long = cum_mut_counts.reset_index().melt(id_vars='WEEK', value_name='Markers (Cumulative)')
                    weekly_long = weekly_mut_counts.reset_index().melt(id_vars='WEEK', value_name='Markers (Week)')
                    
                    plot_df = pd.merge(cum_long, weekly_long, on=['WEEK', 'AA_MUTATION'])
                    plot_df = pd.merge(plot_df, weekly_totals, on='WEEK')
                    
                    plot_df['Freq_Weekly'] = (plot_df['Markers (Week)'] / plot_df['Samples (Week)'].replace(0, np.nan)) * 100
                    plot_df['Freq_Cum'] = (plot_df['Markers (Cumulative)'] / plot_df['Samples (Cumulative)'].replace(0, np.nan)) * 100
                    
                    plot_df = pd.merge(plot_df, hover_info, on='AA_MUTATION', how='left')

                    display_name = plot_name.replace('_', ' - ')
                    
                    # Add dual y-axes with shared x-axis for weekly and cumulative frequencies, and total samples as secondary y-axis
                    fig = make_subplots(
                        rows=2, cols=1,
                        subplot_titles=("<b>Weekly Frequency</b>", "<b>Cumulative Frequency</b>"),
                        shared_yaxes=True,
                        vertical_spacing=0.12,
                        shared_xaxes=True,
                        specs=[[{"secondary_y": True}], [{"secondary_y": True}]]
                    )

                    # Bars for total samples per week on the first subplot
                    fig.add_trace(
                        go.Bar(
                            x=weekly_totals['WEEK'], 
                            y=weekly_totals['Samples (Week)'],
                            name="Total Samples (n)",
                            marker_color='rgba(180,180,180,0.5)',
                            marker_line_width=0,
                            legend="legend2",
                            hovertemplate="<b>Week:</b> %{x|%V, %Y}<br><b>Total sequenced:</b> %{y} samples<extra></extra>"
                        ), row=1, col=1, secondary_y=True
                    )

                    # Add the gray area for cumulative samples
                    fig.add_trace(
                        go.Scatter(
                            x=weekly_totals['WEEK'], 
                            y=weekly_totals['Samples (Cumulative)'],
                            name="Total Samples (n)",
                            marker_color='rgba(180,180,180,0.3)',
                            fillcolor='rgba(180,180,180,0.3)',
                            fill='tozeroy',
                            legend="legend2",
                            hovertemplate="<b>Week:</b> %{x|%V, %Y}<br><b>Total sequenced:</b> %{y} samples<extra></extra>"
                        ), row=2, col=1, secondary_y=True
                    )

                    unique_mutations = plot_df['AA_MUTATION'].unique()
                    
                    for idx, mut in enumerate(unique_mutations):
                        mut_df = plot_df[plot_df['AA_MUTATION'] == mut]
                        color = color_list[idx % len(color_list)]
                        
                        if "${params.protocol}" == "AVIAN":
                            ref = "H5N1"
                        elif "${params.protocol}" == "HUMAN":
                            ref = "H1N1"
                        else:
                            ref = "Unknown"

                        fig.add_trace(
                            go.Scatter(
                                x=mut_df['WEEK'], y=mut_df['Freq_Weekly'],
                                name=mut, mode='lines+markers',
                                connectgaps=True,
                                line=dict(color=color),
                                legendgroup=mut,
                                customdata=mut_df[['POSITION_REF', 'EFFECT', 'FOUND_IN', 'Markers (Week)', 'Samples (Week)', 'Freq_Weekly', 'AA_MUTATION']],
                                hovertemplate=(
                                    "<b>Week:</b> %{x|%V, %Y}<br>"
                                    "<b>Reference Position (" + ref + " numbering):</b> %{customdata[0]}<br>"
                                    "<b>Mutation:</b> %{customdata[6]}<br>"
                                    "<b>Effect(s):</b> %{customdata[1]}<br>"
                                    "                 <b>Found in:</b>  %{customdata[2]}<br>"
                                    "<b>Occurrence (Weekly):</b> %{customdata[3]}/%{customdata[4]} samples (%{customdata[5]:.2f}%)<br>"
                                    "<extra></extra>"
                                )
                            ), row=1, col=1, secondary_y=False
                        )

                        fig.add_trace(
                            go.Scatter(
                                x=mut_df['WEEK'], y=mut_df['Freq_Cum'],
                                name=mut, mode='lines+markers',
                                connectgaps=True,
                                line=dict(color=color),
                                legendgroup=mut, showlegend=False,
                                customdata=mut_df[['POSITION_REF', 'EFFECT', 'FOUND_IN', 'Markers (Cumulative)', 'Samples (Cumulative)', 'Freq_Cum', 'AA_MUTATION']],
                                hovertemplate=(
                                    "<b>Week:</b> %{x|%V, %Y}<br>"
                                    "<b>Reference Position (" + ref + " numbering):</b> %{customdata[0]}<br>"
                                    "<b>Mutation:</b> %{customdata[6]}<br>"
                                    "<b>Effect(s):</b> %{customdata[1]}<br>"
                                    "                 <b>Found in:</b>  %{customdata[2]}<br>"
                                    "<b>Occurrence (Cumulative):</b> %{customdata[3]}/%{customdata[4]} samples (%{customdata[5]:.2f}%)<br>"
                                    "<extra></extra>"
                                )
                            ), row=2, col=1, secondary_y=False
                        )

                    fig.update_yaxes(title_text="Frequency (%)", range=[-5, 105], row=1, col=1, secondary_y=False)
                    fig.update_yaxes(title_text="Cumulative Frequency (%)", range=[-5, 105], row=2, col=1, secondary_y=False)
                    # Set secondary y-axes                    
                    fig.update_yaxes(title_text="Total Samples (N)", showgrid=False, rangemode='nonnegative', fixedrange=True, row=1, col=1, secondary_y=True)
                    fig.update_yaxes(title_text="Total Samples (N)", showgrid=False, rangemode='nonnegative', fixedrange=True, row=2, col=1, secondary_y=True)
                    
                    fig.update_xaxes(tickformat="Week %V<br>%Y", showticklabels=True)
                    
                    fig.update_layout(
                        title=dict(
                            text=f"<b>Frequency in Time Report: {display_name}</b>",
                            font=dict(size=22)
                        ),
                        height=900,
                        hovermode='closest',
                        template="plotly_white",
                        barmode='overlay',
                        legend=dict(
                            title=dict(text="<b>Markers</b>", font=dict(size=14)),
                            x=1.02, 
                            y=1.0, 
                            xanchor="left",
                            yanchor="top",
                        ),
                        legend2=dict(
                            title=dict(text="<b>Total Samples</b>", font=dict(size=14)),
                            orientation="h", 
                            x=0.5, 
                            y=-0.1, 
                            xanchor="center",
                            yanchor="top",
                        )
                    )

                    fig.write_html(f"FrequencyEvolution/{segment}/evolution_{plot_name}.html")

    if __name__ == "__main__":
        generate_plots('${excel_file}', '${meta_file}')
    """
}