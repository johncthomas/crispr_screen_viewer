#!/usr/bin/env python
import pathlib, os
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash_html_components import Div
from dash.dependencies import Input, Output, State
import pandas as pd
import plotly.graph_objects as go
from dash_table import DataTable
from dash_table.Format import Format, Scheme
#import plotly.express as px
from argparse import ArgumentParser
import typing
from typing import Collection, Union, Dict
from crispr_screen_viewer.util import index_of_true, DataSet
#from dashapp import app as application


VERSION = '1.1.0'

# *updates
# 1.0.3 data_version is an argument
# 1.0.4 stat selector: mixed is now a separate selection
# 1.1 data source selection

#todo filter by suppresor/enhancer
#todo per gene volcano plot
#todo y-axis title
#todo


#from https://sashamaps.net/docs/resources/20-colors/, 99%,
#todo test all these colours
colours = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990',
           '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#333333', ]


def launch_msgv(source_directory:Union[str, os.PathLike], port, debug):
    """A Dash app for showing results from screens for specific screens."""
    print(source_directory)
    data_set = DataSet(source_directory)

    def get_colour_map(list_of_things):
        return {thing:colours[i%len(colours)] for i, thing in enumerate(list_of_things)}

    app = dash.Dash(__name__, external_stylesheets= ['https://codepen.io/chriddyp/pen/bWLwgP.css'])

    # **FIGURE**
    fig = go.Figure()

    # **CONTROLS**
    get_lab_val = lambda arr: [{'label': v, 'value':v} for v in arr]

    # Two scores:
    #   effect size (LFC or NormZ) that is the Y-axis value on the graphs
    #   FDR, that is used for filtering comparisons to be shown
    # Scores can be selected from single analysis type, or mixed.
    # initial radio buttons will list all single types (mageck & drugz at the mo) and mixed
    # selecting mixed pops up new pair of radio buttons allowing selection of source for
    # both effect size and significance values
    selector_style = {'display': 'inline-block',   'border-collapse':'separate','border-spacing':'15px 50px'}
    none_style = {'display':'none'}
    primary_stat_selector = [
        Div([
            html.Label('Analysis type_____', htmlFor='analysis-selector'),
            dcc.RadioItems(
                id='analysis-selector',
                options=[
                    {'label':'DrugZ', 'value':'drugz'},
                    {'label':'MAGeCK',  'value':'mageck'},
                    {'label':'Mixed...', 'value':'mixed'}
                ],
                value='drugz',
            )
        ], style=selector_style),

    ]

    mixed_selectors =  [Div([
        Div([
            html.Label('Effect size_____', htmlFor='score-selector'),
            dcc.RadioItems(
                id='score-selector',
                options=[{'label':'NormZ', 'value':'drz'},
                         {'label':'LFC',  'value':'mag'}],
                value='drz',
            ),
        ], style=selector_style),

        Div([
            html.Label('FDR source', htmlFor='fdr-selector', ),
            dcc.RadioItems(
                id='fdr-selector',
                options=[
                    {'label':'DrugZ', 'value':'drz'},
                    {'label':'MAGeCK',  'value':'mag'}
                ],
                value='drz',
            )
        ], style=selector_style)
    ],  id='mixed-div', style=none_style)]

    data_source_selector = Div([
        dcc.Checklist(
            id='data-source-selector',
            options=get_lab_val(data_set.data_sources),
            value=data_set.data_sources, # set all selected by default
            labelStyle={'display':'inline-block'}
        ),

        html.P([''], id='missing-datasets')
    ])

    # the graph object and things above the graph object
    graph_and_data_selection_div = Div([
        html.H1("Multi-Gene Screen Viewer"),
        Div([
            dcc.Checklist('show-boxplots', options=[{'label':'Show boxplots', 'value':'show-boxplots'}], value=['show-boxplots']),
        ]),

        Div(primary_stat_selector+mixed_selectors),

        data_source_selector,

        Div([dcc.Graph(
            id='gene-violins',
            figure=fig,
            style={'height':'800px'}
        )], ),
    ])

    # this is also used for output of one function, so is defined once here
    order_by_categories = ['Mean score', 'Treatment', 'Experiment ID']

    control_bar = Div([
        Div([
            html.Label('Maximum FDR:', htmlFor='fdr-threshold'),
            dcc.Input('fdr-threshold', type='number', min=0, max=2, step=0.01, value=0.2),
        ], style={'width':'120px', 'display':'inline-block'}),
        Div([
            html.Label('Order by:', htmlFor='order-by'),
            dcc.Dropdown('order-by', value=order_by_categories[0],
                         options=get_lab_val(order_by_categories)),
        ], style={'width':'150px', 'display':'inline-block'}),
        Div([
            html.Label('Colour by:', htmlFor='color-by'),
            dcc.Dropdown(
                'color-by', value='Treatment',
                options=get_lab_val(['Treatment', 'Cell line', 'Experiment ID', 'Library', 'KO'])
            ),
        ]),
    ])

    # select genes by name, and comparisons FDR in selected samples
    gene_selector = Div([
        html.P('Select genes:'),
        dcc.Dropdown('gene-selector', value=[], options=get_lab_val(data_set.genes), multi=True),
    ])

    # ability to filter comparisons based on their metadata
    def get_filter_checklist(column):
        uniques = sorted(data_set.metadata[column].unique())
        get_opts = lambda k: [{'label':k+'  |  ', 'value':k} for k in uniques]
        lab = f"filter-{column.lower().replace(' ', '-')}"
        return dcc.Checklist(
            id=lab,
            options=get_opts(column),
            value=uniques, # set all selected by default
            labelStyle={'display':'inline-block'},
        )
    # Generate individual checklists for metadata columns
    sample_filter = Div([], 'sample-filters')
    for col in ('Treatment', 'Cell line', 'KO'):
        sample_filter.children.extend(
            [html.P(col+':'),
             get_filter_checklist(col)]
        )

    # **DATATABLE**
    def create_datatable(data_df=None):
        #todo make numfmt work
        #numfmt = Format(precision=3, scheme=Scheme.decimal_or_exponent)
        #formatted = Format()
        #numfmt = formatted.precision(3)

        # conditional formatting for the FDR
        thresholds = [0.6, 0.3, 0.1, -0.1]
        fdr_colours = ['#ff0000', '#ff9933', '#ffff00', '#66ff33']
        prev_fdr = 1
        cond_fmts = []
        for fdr, clr in zip(thresholds, fdr_colours):
            fmt = {
                'if':{
                    'filter_query': f'{prev_fdr} <= {{FDR}} < {fdr}',
                    'column_id':'FDR'
                },
                'backgroundColor':clr
            }
            cond_fmts.append(fmt)
            prev_fdr = fdr

        # just a place holder to start with
        if data_df is None:
            return DataTable(
                id='table',
                columns=[{'name':c, 'id':c, } for c in data_set.metadata.columns],
            )

        return DataTable(
            id='table',
            # 'format':get_fmt(c)
            columns=[{'name':c, 'id':c, } for c in data_df.columns],
            data=data_df.to_dict('records'),
            export_format='csv',
            style_data_conditional=cond_fmts,
        )

    table = Div([create_datatable()], id='table-div', className="u-full-width")

    # put it all together
    app.layout = Div([
        graph_and_data_selection_div,
        control_bar,
        html.Br(),
        gene_selector,
        html.Br(),
        sample_filter,
        html.Br(),
        table,
    ])

    @app.callback(
        [Output('score-selector', 'value'),
         Output('fdr-selector', 'value'),
         Output('mixed-div', 'style')],

        [Input('analysis-selector', 'value')],

        [State('score-selector', 'value'),
         State('fdr-selector', 'value')]
    )
    def select_stats_primary(selection, curr_score, curr_fdr):
        if selection == 'mixed':
            return curr_score, curr_fdr, selector_style
        else:
            # this format is kind of vestigial, should use 3-letter abrevs directly
            score, fdr = {'drugz': ('drz', 'drz'),
                          'mageck':('mag', 'mag')}[selection]
            return score, fdr, none_style

    # Define callback to update graph
    @app.callback(
        [Output('gene-violins', 'figure'),
         Output('table-div', 'children'),
         Output('order-by', 'options'),
         Output('missing-datasets', 'children')],

        [Input('score-selector', 'value'),
         Input('fdr-selector', 'value'),
         Input('data-source-selector', 'value'),

         Input('gene-selector', 'value'),
         Input('show-boxplots', 'value'),
         Input('fdr-threshold', 'value'),

         Input('filter-treatment', 'value'),
         Input('filter-cell-line', 'value'),

         Input('order-by', 'value'),
         Input('color-by', 'value')]
    )
    def update_figure(score_source, fdr_source, selected_data_sources,
                      selected_genes, show_boxplots, fdr_thresh,
                      filter_treat, filter_cell,
                      order_by, color_by):

        data_tabs = data_set.get_score_fdr(score_source, fdr_source, selected_data_sources)

        # Identify data sources that are unavailable for selected analyses
        available_sources = data_set.metadata.loc[data_tabs['fdr'].columns, 'Source'].unique()
        missing_data_sources = [d for d in selected_data_sources if d not in available_sources]
        if missing_data_sources:
            missing_data_sources = "Data sources unavailable for selected analyses: "+', '.join(missing_data_sources)
        else:
            missing_data_sources = 'All data sources available'

        # get boolean masks for which comparisons to include in the charts
        # first filter out deselected cell lines and treatments
        comparison_mask = (
                data_set.metadata['Treatment'].isin(filter_treat) &
                data_set.metadata['Cell line'].isin(filter_cell)
        )
        # then by FDR threshold
        fdr_mask = (data_tabs['fdr'].loc[selected_genes] < fdr_thresh).any()
        filtered_scores = data_tabs['score'].loc[selected_genes, (fdr_mask & comparison_mask)]

        # determine order of plots
        if order_by in selected_genes:
            trace_order = filtered_scores.loc[order_by].sort_values().index.values
        elif order_by == 'Mean score':
            trace_order = filtered_scores.mean().sort_values().index.values
        elif order_by in order_by_categories[1:]:
            # subset the metadata to included comps and then sort by the order_by
            trace_order = data_set.metadata.loc[filtered_scores.columns, order_by].sort_values().index
        else: # this shouldn't happen, though maybe I should have a "whatever" order option
            trace_order = filtered_scores.columns.values

        # assemble the figure
        fig = go.Figure()
        # plot the gene scatters
        for gn in selected_genes:
            fig.add_trace(
                go.Scatter(
                    x=trace_order, y=filtered_scores.loc[gn, trace_order],
                    mode='markers', name=gn,
                    marker={'size': 15, 'line':{'width':2, 'color':'DarkSlateGrey'}, 'symbol':'hexagram'})
            )
        # add the boxplot traces if required
        if show_boxplots:
            colour_map = get_colour_map(data_set.metadata.loc[:, color_by].unique())
            for col in trace_order:
                ys = data_tabs['score'][col]
                xs = ys[:]
                xs[:] = col
                color_group = data_set.metadata.loc[col, color_by]
                fig.add_trace(
                    go.Box(x=xs, y=ys, name=color_group, boxpoints=False,
                           line=dict(color=colour_map[color_group]))
                )

        # create the DataTable
        selected_fdr = data_tabs['fdr'].loc[filtered_scores.index, trace_order]
        filtered_scores = filtered_scores.reindex(columns=trace_order)
        filtered_scores.index = filtered_scores.index.map(lambda x: x+' (Effect size)')
        selected_fdr.index = selected_fdr.index.map(lambda x: x+' (FDR)')
        selected_stats = pd.concat([filtered_scores, selected_fdr], sort=False).T

        selected_stats = selected_stats.applymap(lambda n: f"{n:.3}")
        selected_stats.insert(0, 'Comparison', selected_stats.index)
        cols_oi = ['Experiment ID', 'Treatment', 'Dose', 'KO', 'Growth inhibition %', 'Days grown',
                   'Cell line', 'Library']
        selected_metadata = data_set.metadata.loc[selected_stats.index, cols_oi]
        data_table_data = pd.concat([selected_stats, selected_metadata], axis=1)
        dtable = create_datatable(data_table_data)

        sort_by_opts = get_lab_val(order_by_categories+ selected_genes)

        return fig, dtable, sort_by_opts, [missing_data_sources]

    app.run_server(host='0.0.0.0', port=port, debug=debug)

if __name__ == '__main__':
    print('version:', VERSION)
    parser = ArgumentParser(description='MultiGene Screen Viewer')
    parser.add_argument(
        '-p', '--port', metavar='PORT',
        help='Port used to serve the charts'
    )
    parser.add_argument(
        '-d', '--data-version',
        dest='data_version',
        help="Name of the directory within app_data that contains the data from screens."
        # this could potentially take multiple values, with partial data in each...
    )
    parser.add_argument(
        '--debug', action='store_true',
        help='Launch app in debug mode'
    )
    args = parser.parse_args()
    launch_msgv(os.path.join('app_data', args.data_version), args.port, args.debug)
