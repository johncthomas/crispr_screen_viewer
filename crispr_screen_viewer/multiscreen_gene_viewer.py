import pandas as pd
from dash import dash, dcc, html, Input, Output, State
Div = html.Div
import plotly.graph_objs as go
import pathlib, os
from typing import Collection, Union, Dict
#from crispr_screen_viewer.functions_etc import index_of_true, DataSet
from crispr_screen_viewer.shared_components import (
    create_datatable,
    get_lab_val,
    get_reg_stat_selectors,
    get_stat_source_selector,
    colours,
    big_text_style,
    timepoint_labels,
    LOG,
)

from crispr_screen_viewer.functions_etc import (
    doi_to_link,
    datatable_column_dict
)

border_style = {'border': '4px solid #3DD178',
                'border-radius': '14px'}



def initiate(app, data_set, public_version=True) -> Div:
    """Register callbacks to app, generate layout"""

    # if public:
    #     source_display = 'none'
    # else:
    #     source_display = 'inline-block'

    graphid = 'msgv'


    # Graph and table
    graph = dcc.Graph(
        id='msgv-gene-violins',
        figure=go.Figure(),
        style={'height':'800px'}
    )

    table = Div([create_datatable(columns_if_no_df=data_set.comparisons.columns)],
                id='msgv-table-div', className="u-full-width", style={'margin-bottom':'30px'})

    # select genes by name, and comparisons FDR in selected samples
    gene_selector = Div(
        children=[
            html.P('Select genes:', style=big_text_style),
            dcc.Dropdown(id='msgv-gene-selector', placeholder='Select genes', value=[],
                         options=get_lab_val(data_set.genes), multi=True),
            get_stat_source_selector('msgv', 'Analysis:')
        ],

        style={'margin-bottom': '15px'})


    # the graph/table object and the data prefilters
    graph_and_data_selection_div = Div([
        html.H1("Multi-Screen Gene Viewer"),
        html.P("Select your gene(s) of interest. Comparisons that show significant results "
                "(below adjustable FDR) for at least one selected gene are shown below. "
                "Box plots give the overall distribution of scores, and stars show specific genes."),
        Div([
            dcc.Checklist(
                id='msgv-show-boxplots',
                options=[{'label':'Show boxplots', 'value':'msgv-show-boxplots'}],
                value=['msgv-show-boxplots']
            ),
        ], style={'display':'none'}),

        gene_selector,
        # This Tabs will show either graph or datatable
        Div([dcc.Tabs(id='output-tabs', value='graph',
                  children=[
                      dcc.Tab([graph], label='Graph', value='graph'),
                      dcc.Tab([table], label='Table', value='table')
                  ])
             ]),
        # Div([
        #     Div(get_reg_stat_selectors(app, id_prefix='msgv'),
        #         style={'display':source_display, 'width':'170px','vertical-align':'top'}),
        # ])
    ])

    # this is also used for output of one function, so is defined once here
    order_by_categories = ['Mean score', 'Treatment', 'Experiment ID']
    colourable_categories = ['Treatment', 'Cell line', 'Experiment ID', 'Library', 'KO']
    control_bar = Div([
        # Div([
        #     html.P('Use controls below to filter which comparisons are shown. ')
        # ]),
        Div([
            html.Label('Maximum FDR:', htmlFor='msgv-fdr-threshold'),
            dcc.Input(id='msgv-fdr-threshold', type='number', min=0, max=2, step=0.01, value=0.1),
        ], style={ 'display':'inline-block', 'width':'135px'}),
        Div([
            html.Label('Order by:', htmlFor='msgv-order-by'),
            dcc.Dropdown(id='msgv-order-by', value=order_by_categories[0],
                         options=get_lab_val(order_by_categories)),
        ], style={'width':'150px', 'display':'inline-block', 'vertical-align':'top'}),
        Div([
            html.Label('Colour by:', htmlFor='msgv-color-by'),
            dcc.Dropdown(
                id='msgv-color-by', value='Treatment',
                options=get_lab_val(colourable_categories)
            ),
        ], style={'width':'150px', 'display':'inline-block', 'vertical-align':'top'}),
    ])

    # get color map, asssiging colors to the most common values first, so that
    #   common things have different colours.
    def get_colour_map(list_of_things):
        return {thing:colours[i%len(colours)] for i, thing in enumerate(list_of_things)}
    colour_maps = {}
    for color_by in colourable_categories:
        ordered_things = data_set.comparisons.loc[:, color_by].value_counts().index
        cm = get_colour_map(ordered_things)
        colour_maps[color_by] = cm

    # make the dropdowns for filtering
    filter_dropdowns = []
    filter_cols = ['Treatment', 'Experiment ID', 'KO', 'Cell line',
                   'Library', 'Timepoint']

    if not public_version:
        filter_cols.append('Source')

    def filter_options(col):
        vals = data_set.comparisons[col]
        if vals.isna().any():
            LOG.warning(f'Column "{col}" used for data filtering contains NaN, this may be a problem.')
        sorted_vals = sorted(vals.unique())
        if col != 'Timepoint':
            opts = [{'label':v, 'value':v} for v in sorted_vals]
        else:
            opts = [{'label':timepoint_labels[v], 'value':v} for v in sorted_vals]
        return opts

    for col in filter_cols:
        filter_dropdowns.append(
            html.Div([dcc.Dropdown(
                id='msgv-'+col,
                placeholder='Filter by '+col,
                multi=True,
                style={'height':'80px', 'width':'220px'},
                value=[] if col != 'Timepoint' else ['endpoints'],
                options=filter_options(col),
            )], style={'display':'inline-block'})
        )

    # put it all together
    msgv_layout = Div([
        graph_and_data_selection_div,
        html.Br(),
        control_bar,
        html.Br(),
        Div(filter_dropdowns, ),
        #html.Br(),
        #table,
    ])


    # Define callback to update graph
    @app.callback(
        [Output('msgv-gene-violins', 'figure'),
         Output('msgv-table-div', 'children'),
         Output('msgv-order-by', 'options'),],

        [Input('msgv-stat-source-selector', 'value'),

         Input('msgv-gene-selector', 'value'),
         Input('msgv-show-boxplots', 'value'),
         Input('msgv-fdr-threshold', 'value'),

         Input('msgv-order-by', 'value'),
         Input('msgv-color-by', 'value'),
         ] + [Input('msgv-'+cid, 'value') for cid in filter_cols]
    )
    def update_figure(score_type,
                      selected_genes, show_boxplots, fdr_thresh,
                      order_by, color_by, *filters):

        data_tabs = data_set.get_score_fdr(score_type, score_type, )

        # # Identify data sources that are unavailable for selected analyses
        # # hopefully depreciated - will probably take a different tack if it ends up not everything
        # #  can be analysed by everything
        # available_sources = data_set.comparisons.loc[data_tabs['fdr'].columns, 'Source'].unique()
        # missing_data_sources = [d for d in selected_data_sources if d not in available_sources]
        # if missing_data_sources:
        #     missing_data_sources = "(unavailable with current analysis type: "+', '.join(missing_data_sources)+')'
        # else:
        #     missing_data_sources = '(All selections available for analysis type)'

        # get boolean masks for which comparisons to include in the charts
        # first get comparisons filtered by metadata filters
        comparison_mask = pd.Series(True, index=data_set.comparisons.index)
        for filter_id, values in zip(filter_cols, filters):
            if values:
                comparison_mask = comparison_mask & data_set.comparisons[filter_id].isin(values)

        # then by FDR threshold
        fdr_mask = (data_tabs['fdr'].loc[selected_genes] < fdr_thresh).any()
        filtered_scores = data_tabs['score'].loc[selected_genes, (fdr_mask & comparison_mask)]

        # determine order of plots
        if order_by in selected_genes:
            ordered_comps = filtered_scores.loc[order_by].sort_values().index.values
        elif order_by == 'Mean score':
            ordered_comps = filtered_scores.mean().sort_values().index.values
        elif order_by in order_by_categories[1:]:
            # subset the metadata to included comps and then sort by the order_by
            ordered_comps = data_set.comparisons.loc[filtered_scores.columns, order_by].sort_values().index
        else: # this shouldn't happen, though maybe I should have a "whatever" order option
            ordered_comps = filtered_scores.columns.values

        # assemble the figure
        fig = go.Figure()
        # plot the gene scatters
        trace_numbers = [str(i) for i in range(1, len(ordered_comps)+1)]
        for gn in selected_genes:
            fig.add_trace(
                go.Scatter(
                    x=trace_numbers,
                    y=filtered_scores.loc[gn, ordered_comps],
                    mode='markers', name=gn,
                    marker={'size': 15, 'line':{'width':2, 'color':'DarkSlateGrey'}, 'symbol':'hexagram'}),

            )
        # add the boxplot traces if required
        if show_boxplots:
            for trace_i, comp in enumerate(ordered_comps):
                # values that paramatise the box plot
                ys = data_tabs['score'][comp]

                boxlabels = pd.Series(str(trace_i+1), index=ys.index)

                color_group = data_set.comparisons.loc[comp, color_by]

                fig.add_trace(
                    go.Box(x=boxlabels, y=ys, name=color_group, boxpoints=False,
                           line=dict(color=colour_maps[color_by][color_group]))
                )
        # labels
        fig.update_layout(xaxis_title='Boxplot number',
                          yaxis_title=data_set.score_labels[score_type],)

        # DataTable structure:
        # First column trace numbers,
        #   then score/FDR for each gene selected,
        #   then the metadata.

        # Create the stat columns
        selected_fdr = data_tabs['fdr'].loc[filtered_scores.index, ordered_comps]
        filtered_scores = filtered_scores.reindex(columns=ordered_comps)
        filtered_scores.index = filtered_scores.index.map(lambda x: x+' (Effect size)')
        selected_fdr.index = selected_fdr.index.map(lambda x: x+' (FDR)')
        selected_stats = pd.concat([filtered_scores, selected_fdr], sort=False).T

        selected_stats = selected_stats.applymap(lambda n: f"{n:.3}")
        selected_stats.insert(0, 'Boxplot number', trace_numbers)
        cols_oi = ['Experiment ID', 'DOI', 'Treatment', 'Dose', 'KO', 'Growth inhibition %', 'Days grown',
                   'Cell line', 'Library']
        selected_metadata = data_set.comparisons.loc[selected_stats.index, cols_oi]

        data_table_data = pd.concat([selected_stats, selected_metadata], axis=1)

        dtable = create_datatable(data_table_data)

        sort_by_opts = get_lab_val(order_by_categories+ selected_genes)

        return fig, dtable, sort_by_opts

    return msgv_layout
