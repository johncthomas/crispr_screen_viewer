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
        '--debug', action='store_true',
        help='Launch app in debug mode'
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

    return data_set, args.port, args.debug, args.hide_source_selector

def intiate_app(data_set, hide_source_selector=False, ):
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
    data_set, port, debug, hide_source_selector = parse_clargs()
    LOG.setLevel('DEBUG'); print('LOG level DEBUG')
    import logging
    log = logging.getLogger('werkzeug')
    log.setLevel(logging.ERROR)

    app = intiate_app(data_set, hide_source_selector)
    app.run_server(debug=debug, host='0.0.0.0', port=port)




# This method ended up being weirdly inconsistent, to launch with Gunicorn or
#   whatever, just have a server specific .py that imports initiate_app
#   and creates the server

# def make_falsy_false(var):
#     """Return False if var is a string saying "no" or "false", ignoring
#     capitals. Otherwise returns var.
#     """
#     if type(var) is str:
#         if var.lower() in ('no', 'false'):
#             return False
#     return var
# # if this has been set we get options from the environment
# # otherwise from the command line
# using_env_args = os.getenv('DDRCS', False)
# using_env_args = make_falsy_false(using_env_args)
# if using_env_args:
#     data_set = load_dataset(os.getenv('DDRCS_DATA', None))
#     debug = os.getenv('DDRCS_DEBUG', False)
#     if debug:
#         debug = True
#     port = os.getenv('DDRCS_PORT', 80)
#     # not sure if it needs to be int, but this functions as validation
#     port = int(port)
#
#     # Hide the data selection boxes if it's public
#     private = make_falsy_false(os.getenv('DDRCS_PRIVATE', False))
#     if private:
#         hide_source_selector = False
#     else:
#         hide_source_selector = True
#     print(f"using env args, port = {port}")
# else:
