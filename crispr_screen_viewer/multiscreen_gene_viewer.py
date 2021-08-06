#import pathlib, os
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash_html_components import Div
from dash.dependencies import Input, Output
import pandas as pd
import plotly.graph_objects as go
from dash_table import DataTable
from dash_table.Format import Format, Scheme
#import plotly.express as px
from argparse import ArgumentParser
#from dashapp import app as application


VERSION = '1.0.1'

#todo filter by suppresor/enhancer
#todo per gene volcano plot
#todo y-axis title

# this is outer scope so it doesn't need to be loaded every time
data_date = '20210709'
# `data_selector` will need updating if this changes
data_tables = ['drz_fdr', 'drz_nmz', 'mag_fdr', 'mag_lfc']
APPDATA = {restype:pd.read_csv(f'app_data/{restype}.{data_date}.csv', index_col=0) for restype in data_tables}
GENES = APPDATA[data_tables[0]].index


METADATA = pd.read_csv(f'app_data/metadata.{data_date}.csv', index_col=0)
#METADATA.loc[:, 'KO'] = METADATA.KO.fillna('WT') # this is done earlier now
# Treatment col read as strings, convert to tuples (For hashability), and create printable versions
METADATA.loc[:, 'Treatment_tup'] = METADATA.Treatment.apply(lambda x: tuple(eval(x)))
METADATA.loc[:, 'Treatment'] = METADATA.Treatment_tup.apply(lambda tup: ' & '.join(tup))

#from https://sashamaps.net/docs/resources/20-colors/, 99%,
#todo test all these colours
colours = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990',
           '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#333333', ]


def launch(port, debug):


    def get_colour_map(list_of_things):
        return {thing:colours[i%len(colours)] for i, thing in enumerate(list_of_things)}



    app = dash.Dash(__name__, external_stylesheets= ['https://codepen.io/chriddyp/pen/bWLwgP.css'])

    # **FIGURE**
    fig = go.Figure()



    # **CONTROLS**
    get_lab_val = lambda arr: [{'label': v, 'value':v} for v in arr]

    score_selector = Div([
        Div([
            html.Label('Effect size_____', htmlFor='score-selector'),
            dcc.RadioItems(
                id='score-selector',
                options=[
                    {'label':'NormZ', 'value':'drz_nmz'},
                    {'label':'LFC',  'value':'mag_lfc'}
                ],
                value='drz_nmz'
            )
        ], style={'display': 'inline-block',   'border-collapse':'separate','border-spacing':'15px 50px'}), Div([
            html.Label('FDR source', htmlFor='fdr-selector'),
            dcc.RadioItems(
                id='fdr-selector',
                options=[
                    {'label':'DrugZ', 'value':'drz_fdr'},
                    {'label':'MAGeCK',  'value':'mag_fdr'}
                ],
                value='drz_fdr'
            )
        ], style={'display': 'inline-block',   'border-collapse':'separate','border-spacing':'15px 50px'}),

    ], style={'display': 'inline-block'})


    # the graph object and things above the graph object
    graph = Div([
        html.H1("Cross-screen gene results"),
        Div([
            dcc.Checklist('show-boxplots', options=[{'label':'Show boxplots', 'value':'show-boxplots'}], value=['show-boxplots']),
            score_selector
        ]),
        Div([dcc.Graph(
            id='gene-violins',
            figure=fig,
            style={'height':'800px'}
        )], ),
    ])

    control_bar = Div([
        Div([
            html.Label('Maximum FDR:', htmlFor='fdr-threshold'),
            dcc.Input('fdr-threshold', type='number', min=0, max=1, step=0.01, value=0.2),
        ], style={'width':'120px', 'display':'inline-block'}),
        Div([
            html.Label('Order by:', htmlFor='order-by'),
            dcc.Dropdown('order-by', value='Mean score',
                         options=get_lab_val(['Mean score', 'Treatment', 'Experiment'])),
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
        dcc.Dropdown('gene-selector', value=[], options=get_lab_val(GENES), multi=True),
    ])

    # ability to filter comparisons based on their metadata
    def get_filter_checklist(column):
        if column == 'Treatment':
            uniques = set()
            for unqs in METADATA['Treatment_tup'].unique():
                for u in unqs:
                    uniques.add(u)
            uniques = sorted(list(uniques))
        else:
            try:
                uniques = sorted(METADATA[column].unique())
            except:
                print(column)
                raise

        get_opts = lambda k: [{'label':k+'  |  ', 'value':k} for k in uniques]
        lab = f"filter-{column.lower().replace(' ', '-')}"
        return dcc.Checklist(
            id=lab,
            options=get_opts(column),
            value=uniques, # set all selected by default
            labelStyle={'display':'inline-block'},
        )

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
                #'format':formatted
                columns=[{'name':c, 'id':c, } for c in METADATA.columns],
            )

        # def get_fmt(column):
        #     if '(' in column:
        #         return numfmt
        #     return formatted

        return DataTable(
            id='table',
            # 'format':get_fmt(c)
            columns=[{'name':c, 'id':c, } for c in data_df.columns],
            data=data_df.to_dict('records'),
            export_format='csv',
            style_data_conditional=cond_fmts,
        )

    table = Div([
        create_datatable()
    ], id='table-div', className="u-full-width")

    # put it all together
    app.layout = Div([

        graph,
        control_bar,
        html.Br(),
        gene_selector,
        html.Br(),
        sample_filter,
        html.Br(),
        table,
    ])

    # Define callback to update graph
    @app.callback(
        [Output('gene-violins', 'figure'),
         Output('table-div', 'children'),
         Output('order-by', 'options')],

        [Input('score-selector', 'value'),
         Input('fdr-selector', 'value'),

         Input('gene-selector', 'value'),
         Input('show-boxplots', 'value'),
         Input('fdr-threshold', 'value'),

         Input('filter-treatment', 'value'),
         Input('filter-cell-line', 'value'),

         Input('order-by', 'value'),
         Input('color-by', 'value')]
    )
    def update_figure(score_source, fdr_source,
                      genes, show_boxplots, fdr_thresh,
                      filter_treat, filter_cell,
                      order_by, color_by):

        # get boolean masks for comparisons
        comparison_mask = (
                METADATA['Treatment_tup'].apply(lambda rowtreat: any([t in filter_treat for t in rowtreat])) &
                METADATA['Cell line'].isin(filter_cell)
        )
        fdr_tab = APPDATA[fdr_source]
        score_tab = APPDATA[score_source]
        fdr_mask = (fdr_tab.loc[genes] < fdr_thresh).any()
        selected_data = score_tab.loc[genes, (fdr_mask & comparison_mask)]

        # determine order of plots
        if order_by in genes:
            trace_order = selected_data.loc[order_by].sort_values().index.values
        elif order_by == 'Mean score':
            trace_order = selected_data.mean().sort_values().index.values
        #todo order by exp
        else: # this shouldn't happen, though maybe I should have a "whatever" order option
            trace_order = selected_data.columns.values

        # assemble the figure
        fig = go.Figure()
        # plot the gene scatters
        for gn in genes:
            fig.add_trace(
                go.Scatter(
                    x=trace_order, y=selected_data.loc[gn, trace_order],
                    mode='markers', name=gn,
                    marker={'size': 15, 'line':{'width':2, 'color':'DarkSlateGrey'}, 'symbol':'hexagram'})
            )
        # add the boxplot traces if required
        if show_boxplots:
            colour_map = get_colour_map(METADATA.loc[:, color_by].unique())
            for col in trace_order:
                ys = score_tab[col]
                xs = ys[:]
                xs[:] = col
                color_group = METADATA.loc[col, color_by]
                fig.add_trace(
                    go.Box(x=xs, y=ys, name=color_group, boxpoints=False,
                           line=dict(color=colour_map[color_group]))
                )

        # create the DataTable
        selected_fdr = fdr_tab.loc[selected_data.index, trace_order]
        selected_data = selected_data.reindex(columns=trace_order)
        selected_data.index = selected_data.index.map(lambda x: x+' (Effect size)')
        selected_fdr.index = selected_fdr.index.map(lambda x: x+' (FDR)')
        selected_stats = pd.concat([selected_data, selected_fdr], sort=False).T

        selected_stats = selected_stats.applymap(lambda n: f"{n:.3}")
        selected_stats.insert(0, 'Comparison', selected_stats.index)
        cols_oi = ['Experiment ID', 'Treatment', 'Dose', 'KO', 'Growth inhibition %', 'Days grown',
                   'Cell line', 'Library']
        selected_metadata = METADATA.loc[selected_stats.index, cols_oi]
        data_table_data = pd.concat([selected_stats, selected_metadata], axis=1)
        dtable = create_datatable(data_table_data)

        sort_by_opts = get_lab_val(['Mean score']+genes)

        return fig, dtable, sort_by_opts

    # Run app and display result inline in the notebook
    app.run_server(host='0.0.0.0', port=port, debug=debug)

if __name__ == '__main__':
    print(VERSION)
    parser = ArgumentParser(description='MultiGene Screen Viewer')
    parser.add_argument(
        '-p', '--port', metavar='PORT',
        help='Port used to serve the charts'
    )
    parser.add_argument(
        '--debug', action='store_true',
        help='Launch app in debug mode'
    )
    args = parser.parse_args()
    launch(args.port, args.debug)
