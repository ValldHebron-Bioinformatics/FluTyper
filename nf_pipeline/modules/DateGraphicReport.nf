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
    import re

    def generate_plots(mut_file, meta_file):
        df_mut = pd.read_excel(mut_file, na_filter=False)
        df_meta = pd.read_csv(meta_file, skipinitialspace=True)

        df_meta['DATE'] = pd.to_datetime(df_meta['DATE'], format='%Y-%m-%d')
        df_meta['WEEK'] = df_meta['DATE'].dt.to_period('W').dt.to_timestamp()

        # Calculate Seasons
        iso_cal = df_meta['DATE'].dt.isocalendar()
        s_year = iso_cal.year.where(iso_cal.week >= 40, iso_cal.year - 1)
        df_meta['Season'] = s_year.astype(str) + "-" + (s_year + 1).astype(str)
        df_meta['Season'] = df_meta['Season'].fillna("Unknown Season")

        # Determine global x-axis range
        all_time_start = df_meta['WEEK'].min() - pd.Timedelta(days=7)
        all_time_end = df_meta['WEEK'].max() + pd.Timedelta(days=14)
        all_time_range = [all_time_start.strftime('%Y-%m-%d'), all_time_end.strftime('%Y-%m-%d')] if pd.notnull(all_time_start) else None

        # Determine true 52-week seasonal x-axis ranges regardless of data gaps
        season_ranges = {}
        for season in sorted(df_meta['Season'].dropna().unique()):
            if season == "Unknown Season":
                continue
            try:
                y1 = int(season.split('-')[0])
                y2 = int(season.split('-')[1])
                # Start: Year1 Week 40 Monday. End: Year2 Week 39 Sunday
                s_start = pd.to_datetime(f'{y1}-W40-1', format='%G-W%V-%u')
                s_end = pd.to_datetime(f'{y2}-W39-7', format='%G-W%V-%u')
                season_ranges[season] = [s_start.strftime('%Y-%m-%d'), s_end.strftime('%Y-%m-%d')]
            except Exception:
                # Fallback if season format is unexpected
                s_data = df_meta[df_meta['Season'] == season]
                if not s_data.empty:
                    s_start = s_data['WEEK'].min()
                    s_end = s_data['WEEK'].max()
                    season_ranges[season] = [s_start.strftime('%Y-%m-%d'), s_end.strftime('%Y-%m-%d')]

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
            color_list = ['#a0850c', '#CD733D', '#C84630','#A196B0', '#3561a3', '#676F86', '#3A2D32', '#381983', '#009E73' ]

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
                    
                    # Build season-specific recalculated datasets

                    season_plot_data = {}

                    # All-time data
                    season_plot_data["All Time"] = (
                        plot_df.copy(),
                        weekly_totals.copy()
                    )

                    # Seasonal recalculated data
                    for season, s_range in season_ranges.items():

                        s_start = pd.to_datetime(s_range[0])
                        s_end = pd.to_datetime(s_range[1])

                        # Seasonal totals
                        season_totals = weekly_totals[
                            (weekly_totals['WEEK'] >= s_start) &
                            (weekly_totals['WEEK'] <= s_end)
                        ].copy()

                        # RESET cumulative sample counts
                        season_totals['Samples (Cumulative)'] = \
                            season_totals['Samples (Week)'].cumsum()

                        # Seasonal mutation data
                        season_mut = df_prot[
                            (df_prot['WEEK'] >= s_start) &
                            (df_prot['WEEK'] <= s_end)
                        ]

                        season_week_counts = (
                            season_mut
                            .groupby(['WEEK','AA_MUTATION'])
                            .size()
                            .unstack(fill_value=0)
                        )

                        season_week_counts = season_week_counts.reindex(
                            season_totals['WEEK'],
                            fill_value=0
                        )

                        # RESET mutation cumulative counts
                        season_cum_counts = season_week_counts.cumsum()

                        season_cum_long = (
                            season_cum_counts
                            .reset_index()
                            .melt(
                                id_vars='WEEK',
                                value_name='Markers (Cumulative)'
                            )
                        )

                        season_week_long = (
                            season_week_counts
                            .reset_index()
                            .melt(
                                id_vars='WEEK',
                                value_name='Markers (Week)'
                            )
                        )

                        season_df = pd.merge(
                            season_cum_long,
                            season_week_long,
                            on=['WEEK','AA_MUTATION']
                        )

                        season_df = pd.merge(
                            season_df,
                            season_totals,
                            on='WEEK'
                        )

                        season_df['Freq_Weekly'] = (
                            season_df['Markers (Week)'] /
                            season_df['Samples (Week)'].replace(0,np.nan)
                        ) * 100

                        season_df['Freq_Cum'] = (
                            season_df['Markers (Cumulative)'] /
                            season_df['Samples (Cumulative)'].replace(0,np.nan)
                        ) * 100

                        season_df = pd.merge(
                            season_df,
                            hover_info,
                            on='AA_MUTATION',
                            how='left'
                        )

                        season_plot_data[season] = (
                            season_df,
                            season_totals
                        )
                    # Calculate fixed maximum bounds to lock the secondary y-axes
                    max_weekly_samples = weekly_totals['Samples (Week)'].max()
                    max_cum_samples = weekly_totals['Samples (Cumulative)'].max()
                    y2_max = max_weekly_samples * 1.1 if pd.notnull(max_weekly_samples) and max_weekly_samples > 0 else 10
                    y4_max = max_cum_samples * 1.1 if pd.notnull(max_cum_samples) and max_cum_samples > 0 else 10

                    display_name = plot_name.replace('_', ' - ')
                    
                    # Add dual y-axes with shared x-axis
                    fig = make_subplots(
                        rows=2, cols=1,
                        subplot_titles=("<b>Weekly Frequency</b>", "<b>Cumulative Frequency</b>"),
                        shared_yaxes=True,
                        vertical_spacing=0.12,
                        shared_xaxes=True,
                        specs=[[{"secondary_y": True}], [{"secondary_y": True}]]
                    )

                    # Bars for total samples per week
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

                    # Function to extract the number from the mutation string for sorting
                    def extract_mutation_number(mut_string):
                        match = re.search(r'\\d+', str(mut_string))
                        return int(match.group()) if match else float('inf')
                    
                    # Sort unique mutations numerically based on the extracted number
                    unique_mutations = sorted(plot_df['AA_MUTATION'].unique(), key=extract_mutation_number)
                    
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
                                customdata=mut_df[['POSITION_REF', 'EFFECT', 'FOUND_IN', 'Markers (Week)', 'Samples (Week)', 'Freq_Weekly', 'AA_MUTATION']].to_numpy().tolist(),
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
                                customdata=mut_df[['POSITION_REF', 'EFFECT', 'FOUND_IN', 'Markers (Cumulative)', 'Samples (Cumulative)', 'Freq_Cum', 'AA_MUTATION']].to_numpy().tolist(),
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

                    # Dropdown with FULL seasonal recalculation

                    dropdown_buttons = []

                    for season_name, (s_plot_df, s_totals) in season_plot_data.items():

                        if season_name == "All Time":
                            label = "All Time"
                            x_range = all_time_range
                        else:
                            label = f"Season {season_name}"
                            x_range = season_ranges[season_name]

                        # Seasonal y-axis scaling
                        s_week_max = s_totals['Samples (Week)'].max()
                        s_cum_max = s_totals['Samples (Cumulative)'].max()

                        s_y2 = (
                            s_week_max * 1.1
                            if pd.notnull(s_week_max) and s_week_max > 0
                            else 10
                        )

                        s_y4 = (
                            s_cum_max * 1.1
                            if pd.notnull(s_cum_max) and s_cum_max > 0
                            else 10
                        )

                        visible = [True, True]
                        x_updates = [
                            s_totals['WEEK'],
                            s_totals['WEEK']
                        ]

                        y_updates = [
                            s_totals['Samples (Week)'],
                            s_totals['Samples (Cumulative)']
                        ]

                        custom_updates = [
                            None,
                            None
                        ]

                        for mut in unique_mutations:

                            mut_df = s_plot_df[
                                s_plot_df['AA_MUTATION'] == mut
                            ]

                            present = (
                                mut_df['Markers (Week)'].sum() > 0 or
                                mut_df['Markers (Cumulative)'].sum() > 0
                            )

                            visible.extend([present, present])

                            x_updates.extend([
                                mut_df['WEEK'],
                                mut_df['WEEK']
                            ])

                            y_updates.extend([
                                mut_df['Freq_Weekly'],
                                mut_df['Freq_Cum']
                            ])

                            custom_updates.extend([
                                mut_df[[
                                    'POSITION_REF',
                                    'EFFECT',
                                    'FOUND_IN',
                                    'Markers (Week)',
                                    'Samples (Week)',
                                    'Freq_Weekly',
                                    'AA_MUTATION'
                                ]].to_numpy().tolist(),

                                mut_df[[
                                    'POSITION_REF',
                                    'EFFECT',
                                    'FOUND_IN',
                                    'Markers (Cumulative)',
                                    'Samples (Cumulative)',
                                    'Freq_Cum',
                                    'AA_MUTATION'
                                ]].to_numpy().tolist()
                            ])

                        dropdown_buttons.append(
                            dict(
                                label=label,
                                method="update",
                                args=[
                                    {
                                        "visible": visible,
                                        "x": x_updates,
                                        "y": y_updates,
                                        "customdata": custom_updates
                                    },
                                    {
                                        "title.text":
                                            f"<b>Frequency in Time Report: {display_name} - {label}</b>",
                                        "xaxis.range": x_range,
                                        "xaxis2.range": x_range,
                                        "yaxis2.range": [0, s_y2],
                                        "yaxis4.range": [0, s_y4]
                                    }
                                ]
                            )
                        )

                    # Explicitly anchor both Primary (Frequency) and Secondary (Sample Count) Y-axes
                    fig.update_yaxes(title_text="Frequency (%)", range=[-5, 105], fixedrange=True, row=1, col=1, secondary_y=False)
                    fig.update_yaxes(title_text="Cumulative Frequency (%)", range=[-5, 105], fixedrange=True, row=2, col=1, secondary_y=False)
                    
                    fig.update_yaxes(title_text="Total Samples (N)", range=[0, y2_max], showgrid=False, fixedrange=True, row=1, col=1, secondary_y=True)
                    fig.update_yaxes(title_text="Total Samples (N)", range=[0, y4_max], showgrid=False, fixedrange=True, row=2, col=1, secondary_y=True)
                    
                    fig.update_xaxes(tickformat="Week %V<br>%Y", showticklabels=True)
                    
                    fig.update_layout(
                        title=dict(
                            text=f"<b>Frequency in Time Report: {display_name}</b>",
                            font=dict(size=22),
                            x=0.47, y=0.98, xanchor="center", yanchor="top"
                        ),
                        updatemenus=[dict(
                            active=0,
                            buttons=dropdown_buttons,
                            x=0.47,
                            xanchor="center",
                            y=1.15,
                            yanchor="top",
                            direction="down",
                            showactive=True
                        )],
                        height=900,
                        hovermode='closest',
                        template="plotly_white",
                        barmode='overlay',
                        margin=dict(t=160, b=80, l=80, r=80),
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