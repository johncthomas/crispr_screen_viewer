import pandas as pd
from dash import dash, dcc, html, Input, Output, State
Div = html.Div
import plotly.graph_objs as go

import pathlib, os, logging

from typing import Collection, Union, Dict

from crispr_screen_viewer.shared_components import (
    create_datatable,
    get_lab_val,
    get_reg_stat_selectors,
    get_stat_source_selector,
    colours,
    big_text_style,
    timepoint_labels,
    LOG,
    spawn_filter_dropdowns,
)

from crispr_screen_viewer.functions_etc import (
    doi_to_link,
    datatable_column_dict,
    get_metadata_table_columns,
    get_selector_table_filter_keys
)

border_style = {'border': '4px solid #3DD178',
                'border-radius': '14px'}

import dash_bootstrap_components as dbc

PAGE_ID = 'msgv'

def initiate(app, data_set, public=True) -> Div:
    """Register callbacks to app, generate layout"""

    # if public:
    #     source_display = 'none'
    # else:
    #     source_display = 'inline-block'




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
        ],

        style={'margin-bottom': '15px'})


    filter_cols = get_selector_table_filter_keys(public)['comp']

    filter_dropdowns = spawn_filter_dropdowns(
        PAGE_ID,
        'comp',
        filter_cols,
        data_set.comparisons,
        values={'Timepoint':['Matched time points']},
        card_header='Filter results using boxes below.'
    )

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
        ], style={'display':'none'}), # hidden for now, might remove entirely

        gene_selector,
        Div(filter_dropdowns, ),
        # This Tabs will show either graph or datatable
        Div([dcc.Tabs(id='output-tabs', value='graph',
                  children=[
                      dcc.Tab([graph], label='Graph', value='graph', className='data-tab',
                              selected_className='data-tab--selected',),
                      dcc.Tab([table], label='Table', value='table',
                              className='data-tab', selected_className='data-tab--selected',)
                  ])
             ]),
        # Div([
        #     Div(get_reg_stat_selectors(app, id_prefix='msgv'),
        #         style={'display':source_display, 'width':'170px','vertical-align':'top'}),
        # ])
    ], style={'display':'inline-block'})

    ### SIDE PANEL

    stat_source_selectr = get_stat_source_selector('msgv', 'Analysis:', 'Select significance source:')
    fdr_selectr = Div([
        #html.Label('FDR:', htmlFor='msgv-fdr-threshold'),
        dcc.Input(id='msgv-fdr-threshold', type='number', min=0, max=2, step=0.01, value=0.1),
    ], style={ 'display':'block', 'width':'135px'})
    # this is also used for output of one function, so is defined once here
    order_by_categories = ['Mean score', 'Treatment', 'Experiment ID']
    colourable_categories = ['Treatment', 'Cell line', 'Experiment ID', 'Library', 'KO']
    control_bar = Div([
        # Div([
        #     html.P('Use controls below to filter which comparisons are shown. ')
        # ]),

        Div([
            html.Label('Order by:', htmlFor='msgv-order-by'),
            dcc.Dropdown(id='msgv-order-by', value=order_by_categories[0],
                         options=get_lab_val(order_by_categories)),
        ], style={'width':'150px', 'display':'block', 'vertical-align':'top'}),
        Div([
            html.Label('Colour by:', htmlFor='msgv-color-by'),
            dcc.Dropdown(
                id='msgv-color-by', value='Treatment',
                options=get_lab_val(colourable_categories)
            ),
        ], style={'width':'150px', 'display':'block', 'vertical-align':'top'}),
    ])

    side_panel = Div([
        dbc.Card(
            children=[
                dbc.CardHeader('Set maximum FDR:'),
                dbc.CardBody(fdr_selectr),
                stat_source_selectr,
                dbc.CardHeader('Change boxplots:'),
                dbc.CardBody(control_bar)
            ]
        )
    ], style={'display':'inline-block', 'padding-left':'20px'})

    # get color map, asssiging colors to the most common values first, so that
    #   common things have different colours.
    def get_colour_map(list_of_things):
        return {thing:colours[i%len(colours)] for i, thing in enumerate(list_of_things)}
    colour_maps = {}
    for color_by in colourable_categories:
        ordered_things = data_set.comparisons.loc[:, color_by].value_counts().index
        cm = get_colour_map(ordered_things)
        colour_maps[color_by] = cm




    # put it all together
    msgv_layout = Div([
        Div(graph_and_data_selection_div, style={'display':'inline-block'}),
        side_panel,
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
         ] + [Input(f'{PAGE_ID}-comp-filter-{cid}', 'value') for cid in filter_cols]
    )
    def update_figure(score_type,
                      selected_genes, show_boxplots, fdr_thresh,
                      order_by, color_by, *filters):

        data_tabs = data_set.get_score_fdr(score_type, score_type, )

        LOG.debug(f"MSGV: update_figure({score_type}, {selected_genes}, show_boxplots={show_boxplots}, "
                  f"fdr_thres={fdr_thresh}, ...")

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
        LOG.debug(f'update_figure: filters={filters} filter_cols={filter_cols}')
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
        cols_oi = get_metadata_table_columns(public, PAGE_ID)['comp']

        selected_metadata = data_set.comparisons.loc[selected_stats.index, cols_oi]

        data_table_data = pd.concat([selected_stats, selected_metadata], axis=1)

        dtable = create_datatable(data_table_data)

        sort_by_opts = get_lab_val(order_by_categories+ selected_genes)

        return fig, dtable, sort_by_opts

    return msgv_layout

def launch_page(source, port, debug):
    from crispr_screen_viewer.functions_etc import DataSet

    # pycharm debug doesn't like __name__ here for some reason.
    app = dash.Dash('comparison_maker', external_stylesheets=[dbc.themes.BOOTSTRAP])

    if debug:
        LOG.setLevel(logging.DEBUG)

    source_directory = pathlib.Path(source)

    data_set = DataSet(source_directory)

    app.layout = initiate(app, data_set, public=True)

    app.run_server(debug=debug, host='0.0.0.0', port=int(port), )

if __name__ == '__main__':
    from crispr_screen_viewer.functions_etc import get_cmdline_options
    clargs = get_cmdline_options()
    launch_page(*clargs)

