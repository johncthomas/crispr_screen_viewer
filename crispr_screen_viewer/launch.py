#!/usr/bin/env python

from argparse import ArgumentParser
from crispr_screen_viewer import multiscreen_gene_viewer, screen_explorer#, comparison_maker
from functions_etc import DataSet
from crispr_screen_viewer.shared_components import (
    external_stylesheets,
    # get_lab_val,
    # get_reg_stat_selectors,
    # get_annotation_dicts,
    # big_text_style,
    LOG
)
from dash import dash, html, dcc, Input, Output, State
Div = html.Div
import flask

import pathlib, pickle, os


if __name__ == '__main__':

    launcher_parser = ArgumentParser(add_help=False)

    launcher_parser.add_argument(
        '-p', '--port', metavar='PORT',
        help='Port used to serve the app',
        required=True,
    )

    launcher_parser.add_argument(
        '-d', '--data-path',
        dest='data_path',
        help="Name of the directory or pickle file that contains the screens' data.",
        required=True,
    )
    launcher_parser.add_argument(
        '--debug', action='store_true',
        help='Launch app in debug mode'
    )

    launcher_parser.add_argument(
        '--hide-source-selector', action='store_true',
        help="Don't hide the data-source and analysis-type selectors."
             " In the future analysis-type might have its own option."
    )

    server = flask.Flask(__name__)

    #data_set = DataSet(source_directory)

    parser = ArgumentParser(parents=[launcher_parser],
                            description="Dash app for exploring screen data.",
                            add_help=True,)


    args = parser.parse_args()

    if os.path.isfile(args.data_path):
        LOG.info('args.data_path is a file, assuming pickle and loading.')
        with open(args.data_path, 'rb') as f:
            data_set = pickle.load(f)
    else:
        data_set = DataSet(pathlib.Path(args.data_path))

    metadata = data_set.comparisons

    app = dash.Dash(__name__, external_stylesheets=external_stylesheets, server=server)
    server = app.server

    # register the callbacks
    msgv_layout = multiscreen_gene_viewer.init_msgv(app, data_set, hide_data_selectors=args.hide_source_selector)
    se_layout = screen_explorer.init_msgv(app, data_set, hide_data_selectors=args.hide_source_selector)

    landing_page = Div([
        html.H1('DDR CRISPR screens data explorer'), html.Br(),
        html.H3('Search Genes:'), html.Br(),
        html.P("Search gene names to find if they have been significantly enriched/depleted in experiments."),  html.Br(),
        html.H3('Explore screens:'), html.Br(),
        html.P('Find experiments, filtering by things such as treatment used, and view the results of that experiment.'), html.Br(),
    ])

    # links change the url, a callback detects changes to this
    sidebar = Div(
        [dcc.Link('Search genes', href='/gene-explorer'),
         html.Br(),
         dcc.Link('Explore screens', href='/screen-explorer')],
        className='sidenav',
        id='sidebar'
    )

    graphs_div = Div([], id='graph-div', className='graphbox')

    app.layout = Div([
        dcc.Location(id='url', refresh=False),
        sidebar,
        graphs_div,
    ])

    @app.callback(
        Output('graph-div', 'children'),
        [Input('url', 'pathname')],
    )
    def change_div(pathname):
        if pathname == '/gene-explorer':
            return msgv_layout
        elif pathname == '/screen-explorer':
            return se_layout
        elif not pathname or pathname == '/':
            return html.P('Landing page')
        else:
            return html.P('404')

    app.run_server(debug=args.debug, host='0.0.0.0', port=args.port)

