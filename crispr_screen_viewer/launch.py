#!/usr/bin/env python
import sys
from argparse import ArgumentParser
from crispr_screen_viewer import multiscreen_gene_viewer, screen_explorer, comparison_maker
from crispr_screen_viewer.functions_etc import DataSet
from crispr_screen_viewer.shared_components import (
    #external_stylesheets,
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


def load_dataset(paff):
    """If paff is a dir, the dataset is constructed from the files
    within, otherwise it is assumed to be a pickle."""
    if os.path.isfile(paff):
        LOG.info('args.data_path is a file, assuming pickle and loading.')
        with open(paff, 'rb') as f:
            data_set = pickle.load(f)
    else:
        data_set = DataSet(pathlib.Path(paff))
    return data_set


def parse_clargs():
    """Load dataset from command line options, return (dataset, port, debug)"""
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
        '--app-debug', action='store_true',
        help='Launch Dash app in debug mode, tends to cause issues but might be useful.'
    )
    launcher_parser.add_argument(
        '--debug-messages', action='store_true',
        help='Set log level to debug, print messages describing the internal state of the app. '
             'Also hide Werkzeug messages'
    )
    launcher_parser.add_argument(
        '--hide-source-selector', action='store_true',
        help="Don't hide the data-source and analysis-type selectors."
             " In the future analysis-type might have its own option."
    )


    parser = ArgumentParser(parents=[launcher_parser],
                            description="Dash app for exploring screen data.",
                            add_help=True,)

    args = parser.parse_args()

    data_set = load_dataset(args.data_path)

    return data_set, args.port, args.app_debug, args.debug_messages, args.hide_source_selector

def initiate_app(data_set, hide_source_selector=False, ):
    server = flask.Flask(__name__)


    app = dash.Dash(__name__,  server=server,
                    url_base_pathname='/')

    # register the callbacks and get the page layouts
    msgv_layout = multiscreen_gene_viewer.initiate(app, data_set, hide_data_selectors=hide_source_selector)
    se_layout = screen_explorer.initiate(app, data_set, public_version=hide_source_selector)
    cm_layout = comparison_maker.initiate(app, data_set)


    landing_page = Div([
        html.H1('DDR CRISPR screens data explorer'), html.Br(),
        html.H3('Search Genes:'), html.Br(),
        html.P("Search gene names to find if they have been significantly enriched/depleted in experiments."),  html.Br(),
        html.H3('Explore screens:'), html.Br(),
        html.P('Find experiments, filtering by things such as treatment used, and view the results of that experiment.'), html.Br(),
        html.H3('Compare results:'), html.Br(),
        html.P('Compare the results of two treatments against each other.'), html.Br(),
    ])

    # links change the url, a callback detects changes to this
    sidebar = Div(
        [
            dcc.Link('Home', href='/'),
            html.Br(),
            dcc.Link('Search genes', href='/gene-explorer'),
            html.Br(),
            dcc.Link('Explore screens', href='/screen-explorer'),
            html.Br(),
            dcc.Link('Compare results', href='/comparison-explorer')
        ],
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
        elif pathname == '/comparison-explorer':
            return cm_layout
        elif not pathname or pathname == '/':
            return landing_page
        else:
            return html.P('404, page not found')

    return app



if __name__ == '__main__':
    data_set, port, dash_debug, debug_messages, hide_source_selector = parse_clargs()

    import logging
    if debug_messages:
        werklog = logging.getLogger('werkzeug')
        werklog.setLevel(logging.ERROR)
        LOG.setLevel('DEBUG')

    app = intiate_app(data_set, hide_source_selector)
    app.run_server(debug=dash_debug, host='0.0.0.0', port=port)

