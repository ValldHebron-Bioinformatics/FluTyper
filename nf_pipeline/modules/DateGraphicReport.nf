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
import json

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):  return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.bool_):    return bool(obj)
        if isinstance(obj, np.ndarray):  return obj.tolist()
        return super().default(obj)

def generate_plots(mut_file, meta_file):
    # Load the mutations dataset and guarantee clean string IDs immediately
    df_mut = pd.read_excel(mut_file, na_filter=False)
    if 'SAMPLE_ID' in df_mut.columns:
        df_mut['SAMPLE_ID'] = df_mut['SAMPLE_ID'].astype(str).str.strip()

    df_meta = pd.read_csv(meta_file, skipinitialspace=True)
    
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

    df_meta['DATE'] = pd.to_datetime(df_meta.get('DATE'), errors='coerce')
    df_meta['WEEK'] = df_meta['DATE'].dt.to_period('W').dt.to_timestamp()

    # Calculate Seasons
    iso_cal = df_meta['DATE'].dt.isocalendar()
    s_year = iso_cal.year.where(iso_cal.week >= 40, iso_cal.year - 1)
    df_meta['Season'] = s_year.astype(str) + "-" + (s_year + 1).astype(str)
    df_meta['Season'] = df_meta['Season'].replace('nan-nan', 'Unknown Season').fillna("Unknown Season")

    # Age Group column
    if 'AGE GROUP' in df_meta.columns:
        df_meta['Age_Group'] = df_meta['AGE GROUP'].fillna('Sense dades').astype(str).str.strip()
    elif 'AGE_GROUP' in df_meta.columns:
        df_meta['Age_Group'] = df_meta['AGE_GROUP'].fillna('Sense dades').astype(str).str.strip()
    else:
        df_meta['Age_Group'] = 'Sense dades'

    # Sex column standardization
    if 'SEX' in df_meta.columns:
        df_meta['Sex'] = df_meta['SEX'].astype(str).str.strip()
        df_meta['Sex'] = df_meta['Sex'].replace(
            {'Sense dades': 'Undetermined', 'Unknown': 'Undetermined', 
             'nan': 'Undetermined', 'None': 'Undetermined', '': 'Undetermined'}
        )
    else:
        df_meta['Sex'] = 'Undetermined'

    # Determine global x-axis range
    all_time_start = df_meta['WEEK'].min() - pd.Timedelta(days=7)
    all_time_end = df_meta['WEEK'].max() + pd.Timedelta(days=14)
    all_time_range = [all_time_start.strftime('%Y-%m-%d'), all_time_end.strftime('%Y-%m-%d')] if pd.notnull(all_time_start) else None

    # Determine true 52-week seasonal x-axis ranges
    season_ranges = {}
    raw_seasons = df_meta['Season'].dropna().unique()
    seasons = sorted([s for s in raw_seasons if s != "Unknown Season"], reverse=True)
    
    for season in seasons:
        try:
            y1 = int(season.split('-')[0])
            y2 = int(season.split('-')[1])
            s_start = pd.to_datetime(f'{y1}-W40-1', format='%G-W%V-%u')
            s_end = pd.to_datetime(f'{y2}-W39-7', format='%G-W%V-%u')
            season_ranges[season] = [s_start.strftime('%Y-%m-%d'), s_end.strftime('%Y-%m-%d')]
        except Exception:
            s_data = df_meta[df_meta['Season'] == season]
            if not s_data.empty:
                s_start = s_data['WEEK'].min()
                s_end = s_data['WEEK'].max()
                season_ranges[season] = [s_start.strftime('%Y-%m-%d'), s_end.strftime('%Y-%m-%d')]

    global_weeks = pd.DataFrame({'WEEK': sorted(df_meta['WEEK'].dropna().unique())})

    # Build filter value lists
    age_order  = {'0-2': 0, '3-4': 1, '5-14': 2, '15-65': 3, '>65': 4}
    age_groups = ['All'] + sorted(
        [a for a in df_meta['Age_Group'].unique()
         if str(a).strip() not in ['nan', '', 'None', 'Sense dades']],
        key=lambda x: age_order.get(str(x).strip(), 99)
    )
    
    sex_unique = [g for g in df_meta['Sex'].unique() if str(g).strip() not in ['nan', '', 'None']]
    if 'Undetermined' not in sex_unique:
        sex_unique.append('Undetermined')
    sexes = ['All'] + sorted(sex_unique)

    merge_cols = ['SAMPLE_ID', 'WEEK', 'Age_Group', 'Sex']
    df_all = pd.merge(df_mut, df_meta[merge_cols], on='SAMPLE_ID', how='left')
    df_all['REF_SUBTYPE'] = df_all.get('REF_SUBTYPE', 'Unknown').replace('', 'Unknown').fillna('Unknown').astype(str)
    
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
        color_list = ['#a0850c', '#CD733D', '#C84630','#A196B0', '#3561a3', '#676F86', '#3A2D32', '#381983', '#009E73']

    def build_season_plot_data(df_prot, df_plot_all, unique_mutations):
        weekly_totals = df_plot_all.groupby('WEEK')['SAMPLE_ID'].nunique().reset_index(name='Samples (Week)')
        weekly_totals = pd.merge(global_weeks, weekly_totals, on='WEEK', how='left').fillna({'Samples (Week)': 0})
        weekly_totals['Samples (Cumulative)'] = weekly_totals['Samples (Week)'].cumsum()

        if df_prot.empty:
            dummy_cols = ['WEEK','AA_MUTATION','Markers (Week)','Markers (Cumulative)',
                          'Samples (Week)','Samples (Cumulative)','Freq_Weekly','Freq_Cum',
                          'POSITION_REF','EFFECT','FOUND_IN']
            empty_df = pd.DataFrame(columns=dummy_cols)
            result = {}
            for season in seasons:
                s_range = season_ranges.get(season)
                if not s_range: continue
                s_start = pd.to_datetime(s_range[0])
                s_end   = pd.to_datetime(s_range[1])
                s_tot = weekly_totals[
                    (weekly_totals['WEEK'] >= s_start) & (weekly_totals['WEEK'] <= s_end)
                ].copy()
                s_tot['Samples (Cumulative)'] = s_tot['Samples (Week)'].cumsum()
                result[season] = (empty_df.copy(), s_tot)
            result["All Time"] = (empty_df.copy(), weekly_totals.copy())
            return result, weekly_totals

        weekly_mut_counts = df_prot.groupby(['WEEK', 'AA_MUTATION']).size().unstack(fill_value=0)
        weekly_mut_counts = weekly_mut_counts.reindex(global_weeks['WEEK'], fill_value=0)
        cum_mut_counts = weekly_mut_counts.cumsum()

        cum_long    = cum_mut_counts.reset_index().melt(id_vars='WEEK', value_name='Markers (Cumulative)')
        weekly_long = weekly_mut_counts.reset_index().melt(id_vars='WEEK', value_name='Markers (Week)')

        plot_df = pd.merge(cum_long, weekly_long, on=['WEEK', 'AA_MUTATION'])
        plot_df = pd.merge(plot_df, weekly_totals, on='WEEK')
        plot_df['Freq_Weekly'] = (plot_df['Markers (Week)'] / plot_df['Samples (Week)'].replace(0, np.nan)) * 100
        plot_df['Freq_Cum']    = (plot_df['Markers (Cumulative)'] / plot_df['Samples (Cumulative)'].replace(0, np.nan)) * 100
        plot_df = pd.merge(plot_df, hover_info, on='AA_MUTATION', how='left')

        season_plot_data = {}

        for season in seasons:
            s_range = season_ranges.get(season)
            if not s_range:
                continue
            s_start = pd.to_datetime(s_range[0])
            s_end   = pd.to_datetime(s_range[1])

            season_totals = weekly_totals[
                (weekly_totals['WEEK'] >= s_start) & (weekly_totals['WEEK'] <= s_end)
            ].copy()
            season_totals['Samples (Cumulative)'] = season_totals['Samples (Week)'].cumsum()

            season_mut = df_prot[
                (df_prot['WEEK'] >= s_start) & (df_prot['WEEK'] <= s_end)
            ]
            season_week_counts = (
                season_mut.groupby(['WEEK','AA_MUTATION']).size()
                .unstack(fill_value=0)
                .reindex(season_totals['WEEK'], fill_value=0)
            )
            season_cum_counts = season_week_counts.cumsum()

            season_cum_long  = season_cum_counts.reset_index().melt(id_vars='WEEK', value_name='Markers (Cumulative)')
            season_week_long = season_week_counts.reset_index().melt(id_vars='WEEK', value_name='Markers (Week)')

            season_df = pd.merge(season_cum_long, season_week_long, on=['WEEK','AA_MUTATION'])
            season_df = pd.merge(season_df, season_totals, on='WEEK')
            season_df['Freq_Weekly'] = (season_df['Markers (Week)'] / season_df['Samples (Week)'].replace(0, np.nan)) * 100
            season_df['Freq_Cum']    = (season_df['Markers (Cumulative)'] / season_df['Samples (Cumulative)'].replace(0, np.nan)) * 100
            season_df = pd.merge(season_df, hover_info, on='AA_MUTATION', how='left')

            season_plot_data[season] = (season_df, season_totals)

        season_plot_data["All Time"] = (plot_df.copy(), weekly_totals.copy())
        return season_plot_data, weekly_totals

    def build_view_data(season_plot_data, unique_mutations, display_name):
        view_data = {}
        for season_name, (s_plot_df, s_totals) in season_plot_data.items():
            if season_name == "All Time":
                label   = "All Time"
                x_range = all_time_range
            else:
                label   = f"Season {season_name}"
                x_range = season_ranges[season_name]

            s_week_max = s_totals['Samples (Week)'].max()
            s_cum_max  = s_totals['Samples (Cumulative)'].max()
            s_y2 = s_week_max * 1.1 if pd.notnull(s_week_max) and s_week_max > 0 else 10
            s_y4 = s_cum_max  * 1.1 if pd.notnull(s_cum_max)  and s_cum_max  > 0 else 10

            s_totals  = s_totals.copy()
            s_plot_df = s_plot_df.copy()
            s_totals['WEEK']  = pd.to_datetime(s_totals['WEEK'],  errors='coerce')
            s_plot_df['WEEK'] = pd.to_datetime(s_plot_df['WEEK'], errors='coerce')

            visible       = [True, True]
            x_updates     = [s_totals['WEEK'].dt.strftime('%Y-%m-%d').tolist(),
                             s_totals['WEEK'].dt.strftime('%Y-%m-%d').tolist()]
            y_updates     = [s_totals['Samples (Week)'].tolist(),
                             s_totals['Samples (Cumulative)'].tolist()]
            custom_updates = [[], []]

            for mut in unique_mutations:
                mut_df = s_plot_df[s_plot_df['AA_MUTATION'] == mut]
                present = bool(not mut_df.empty and (mut_df['Markers (Week)'].sum() > 0 or mut_df['Markers (Cumulative)'].sum() > 0))
                
                visible.extend([present, present])
                x_updates.extend([
                    mut_df['WEEK'].dt.strftime('%Y-%m-%d').tolist() if present else [],
                    mut_df['WEEK'].dt.strftime('%Y-%m-%d').tolist() if present else []
                ])
                y_updates.extend([
                    mut_df['Freq_Weekly'].tolist() if present else [],
                    mut_df['Freq_Cum'].tolist() if present else []
                ])
                custom_updates.extend([
                    mut_df[['POSITION_REF','EFFECT','FOUND_IN','Markers (Week)','Samples (Week)','Freq_Weekly','AA_MUTATION']].to_numpy().tolist() if present else [],
                    mut_df[['POSITION_REF','EFFECT','FOUND_IN','Markers (Cumulative)','Samples (Cumulative)','Freq_Cum','AA_MUTATION']].to_numpy().tolist() if present else []
                ])

            view_data[season_name] = {
                "data_update": {"visible": visible, "x": x_updates, "y": y_updates, "customdata": custom_updates},
                "layout_update": {
                    "title.text": f"<b>Frequency in Time Report: {display_name} - {label}</b>",
                    "xaxis.range":  x_range,
                    "xaxis2.range": x_range,
                    "yaxis2.range": [0, s_y2],
                    "yaxis4.range": [0, s_y4]
                }
            }
        return view_data

    for segment, proteins in prot_dict.items():
        os.makedirs(f"FrequencyEvolution/{segment}", exist_ok=True)
        
        for protein in proteins:
            df_base_markers = df_markers[df_markers['PROTEIN'] == protein]
            if df_base_markers.empty: continue

            for plot_name in df_base_markers['PLOT_NAME'].unique():
                df_prot_all    = df_base_markers[df_base_markers['PLOT_NAME'] == plot_name]
                df_plot_all_all = df_all[df_all['PLOT_NAME'] == plot_name]

                def extract_mutation_number(mut_string):
                    match = re.search(r'\\d+', str(mut_string))
                    return int(match.group()) if match else float('inf')

                base_wmc = df_prot_all.groupby(['WEEK','AA_MUTATION']).size().unstack(fill_value=0)
                unique_mutations = sorted(base_wmc.columns.tolist(), key=extract_mutation_number)

                base_spd, base_weekly_totals = build_season_plot_data(
                    df_prot_all, df_plot_all_all, unique_mutations
                )

                max_weekly_samples = base_weekly_totals['Samples (Week)'].max()
                max_cum_samples    = base_weekly_totals['Samples (Cumulative)'].max()
                y2_max = max_weekly_samples * 1.1 if pd.notnull(max_weekly_samples) and max_weekly_samples > 0 else 10
                y4_max = max_cum_samples    * 1.1 if pd.notnull(max_cum_samples)    and max_cum_samples    > 0 else 10

                display_name = plot_name.replace('_', ' - ')

                if "${params.protocol}" == "AVIAN":
                    ref = "H5N1"
                elif "${params.protocol}" == "HUMAN":
                    ref = "H1N1"
                else:
                    ref = "Unknown"

                all_view_sets = {}

                for age_val in age_groups:
                    for sex_val in sexes:
                        dm = df_prot_all.copy()
                        da = df_plot_all_all.copy()
                        if age_val != 'All':
                            dm = dm[dm['Age_Group'] == age_val]
                            da = da[da['Age_Group'] == age_val]
                        if sex_val != 'All':
                            dm = dm[dm['Sex'] == sex_val]
                            da = da[da['Sex'] == sex_val]

                        spd, _ = build_season_plot_data(dm, da, unique_mutations)
                        views = build_view_data(spd, unique_mutations, display_name)
                        
                        for season_name, updates in views.items():
                            combo_key = f"{season_name}|||{age_val}|||{sex_val}"
                            all_view_sets[combo_key] = updates

                default_season = seasons[0] if seasons else "All Time"
                default_view = all_view_sets.get(f"{default_season}|||All|||All")
                plot_df_default, weekly_totals_default = base_spd["All Time"]
                weekly_totals_init = base_weekly_totals

                fig = make_subplots(
                    rows=2, cols=1,
                    subplot_titles=("<b>Weekly Frequency</b>", "<b>Cumulative Frequency</b>"),
                    shared_yaxes=True,
                    vertical_spacing=0.12,
                    shared_xaxes=True,
                    specs=[[{"secondary_y": True}], [{"secondary_y": True}]]
                )

                fig.add_trace(
                    go.Bar(
                        x=weekly_totals_init['WEEK'],
                        y=weekly_totals_init['Samples (Week)'],
                        name="Total Samples (n)",
                        marker_color='rgba(180,180,180,0.5)',
                        marker_line_width=0,
                        legend="legend2",
                        hovertemplate="<b>Week:</b> %{x|%V, %Y}<br><b>Total sequenced:</b> %{y} samples<extra></extra>"
                    ), row=1, col=1, secondary_y=True
                )

                fig.add_trace(
                    go.Scatter(
                        x=weekly_totals_init['WEEK'],
                        y=weekly_totals_init['Samples (Cumulative)'],
                        name="Total Samples (n)",
                        marker_color='rgba(180,180,180,0.3)',
                        fillcolor='rgba(180,180,180,0.3)',
                        fill='tozeroy',
                        legend="legend2",
                        hovertemplate="<b>Week:</b> %{x|%V, %Y}<br><b>Total sequenced:</b> %{y} samples<extra></extra>"
                    ), row=2, col=1, secondary_y=True
                )

                for idx, mut in enumerate(unique_mutations):
                    mut_df = plot_df_default[plot_df_default['AA_MUTATION'] == mut]
                    color  = color_list[idx % len(color_list)]

                    fig.add_trace(
                        go.Scatter(
                            x=mut_df['WEEK'], y=mut_df['Freq_Weekly'],
                            name=mut, mode='lines+markers',
                            connectgaps=True,
                            line=dict(color=color),
                            legendgroup=mut,
                            customdata=mut_df[['POSITION_REF','EFFECT','FOUND_IN','Markers (Week)','Samples (Week)','Freq_Weekly','AA_MUTATION']].to_numpy().tolist(),
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
                            customdata=mut_df[['POSITION_REF','EFFECT','FOUND_IN','Markers (Cumulative)','Samples (Cumulative)','Freq_Cum','AA_MUTATION']].to_numpy().tolist(),
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
                    height=900,
                    hovermode='closest',
                    template="plotly_white",
                    barmode='overlay',
                    margin=dict(t=160, b=80, l=80, r=80),
                    legend=dict(
                        title=dict(text="<b>Markers</b>", font=dict(size=14)),
                        x=1.02, y=1.0, xanchor="left", yanchor="top"
                    ),
                    legend2=dict(
                        title=dict(text="<b>Total Samples</b>", font=dict(size=14)),
                        orientation="h",
                        x=0.5, y=-0.1, xanchor="center", yanchor="top"
                    )
                )

                if default_view:
                    initial_data   = default_view['data_update']
                    initial_layout = default_view['layout_update']
                    for idx, trace in enumerate(fig.data):
                        if idx < len(initial_data['visible']):
                            trace.visible = initial_data['visible'][idx]
                            trace.x       = initial_data['x'][idx]
                            trace.y       = initial_data['y'][idx]
                            if initial_data['customdata'][idx]:
                                trace.customdata = initial_data['customdata'][idx]
                    fig.update_layout(initial_layout)
                else:
                    alt_view = all_view_sets.get("All Time|||All|||All")
                    if alt_view:
                        fig.update_layout(alt_view['layout_update'])

                all_view_sets_js = json.dumps(all_view_sets, cls=NumpyEncoder)

                season_options = "".join([f'<option value="{s}">{s}</option>' for s in (seasons + ["All Time"])])
                age_options = "".join([f'<option value="{a}">{a}</option>' for a in age_groups])
                sex_options = "".join([f'<option value="{g}">{g}</option>' for g in sexes])

                graph_html = fig.to_html(
                    full_html=False, include_plotlyjs='cdn', div_id="plotly-graph"
                )

                html_out = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Frequency Evolution: {display_name}</title>
<style>
body {{ font-family: Arial, sans-serif; text-align: center; margin: 0; padding: 0; }}
.sticky-header {{
    position: sticky; top: 0;
    background-color: rgba(255,255,255,0.96);
    padding: 10px 20px 12px 20px; z-index: 1000;
    box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1);
    border-bottom: 1px solid #eaeaea;
}}
.controls-container {{
    display: flex; justify-content: center; align-items: flex-end;
    gap: 30px; margin: 10px auto 0 auto; flex-wrap: wrap;
}}
.control-group {{
    display: flex; flex-direction: column; align-items: center;
}}
.control-group label {{
    font-size: 10px; font-weight: bold; color: #666; margin-bottom: 4px;
}}
.control-group select {{
    padding: 6px; border-radius: 4px; border: 1px solid #ccc;
    background: white; font-size: 14px; min-width: 140px;
}}
.graph-container {{ padding: 10px 20px 20px 20px; }}
</style>
</head>
<body>

<div class="sticky-header">
    <h2 style="margin: 0 0 4px 0;">Frequency in Time Report: {display_name}</h2>
    <div class="controls-container">
        <div class="control-group">
            <label>SEASON</label>
            <select id="seasonSel" onchange="applyFilters()">
                {season_options}
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
    </div>
</div>

<div class="graph-container">
{graph_html}
</div>

<script>
var allDataSets = {all_view_sets_js};

function applyFilters() {{
    var activeSeason = document.getElementById('seasonSel').value;
    var activeAge = document.getElementById('ageSel').value;
    var activeSex = document.getElementById('sexSel').value;
    var comboKey  = activeSeason + "|||" + activeAge + "|||" + activeSex;

    var view = allDataSets[comboKey];
    if (!view) {{ console.warn('No view set for', comboKey); return; }}

    var gd = document.getElementById('plotly-graph');
    if (!gd || !gd.data) return;

    Plotly.update(gd, view.data_update, view.layout_update);
}}
</script>
</body>
</html>'''

                out_path = f"FrequencyEvolution/{segment}/evolution_{plot_name}.html"
                with open(out_path, "w", encoding="utf-8") as fh:
                    fh.write(html_out)

if __name__ == "__main__":
    generate_plots('${excel_file}', '${meta_file}')
    """
}