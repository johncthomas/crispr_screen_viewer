#!/usr/bin/env python
import sys
from argparse import ArgumentParser
from crispr_screen_viewer import multiscreen_gene_viewer, screen_explorer, comparison_maker
from crispr_screen_viewer.functions_etc import (
    doi_to_link,
    load_dataset
)
from crispr_screen_viewer.shared_components import (
    #external_stylesheets,
    # get_lab_val,
    # get_reg_stat_selectors,
    # get_annotation_dicts,
    # big_text_style,
    LOG
)
from dash import dash, html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
Div = html.Div
import flask
from pathlib import Path
import pathlib, pickle, os
import dash_bootstrap_components as dbc

our_colour = '#3db9c2'


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
        help='Launch Dash app in debug mode, may be less functional but messages might be useful.'
    )
    launcher_parser.add_argument(
        '--debug-messages', action='store_true',
        help='Set log level to debug, print messages describing the internal state of the app. '
             'Also hide Werkzeug messages'
    )
    launcher_parser.add_argument(
        '--public-version', action='store_true',
        help="Don't hide the data-source and analysis-type selectors."
             " In the future analysis-type might have its own option."
    )


    parser = ArgumentParser(parents=[launcher_parser],
                            description="Dash app for exploring screen data.",
                            add_help=True,)

    args = parser.parse_args()

    data_set = load_dataset(args.data_path)

    return data_set, args.port, args.app_debug, args.debug_messages, args.public_version

def initiate_app(data_set, public_version=False,):
    server = flask.Flask(__name__)

    app = dash.Dash(__name__,  server=server,
                    url_base_pathname='/',
                    external_stylesheets=[dbc.themes.BOOTSTRAP])



    # register the callbacks and get the page layouts
    msgv_layout = multiscreen_gene_viewer.initiate(app, data_set, public=public_version)
    se_layout = screen_explorer.initiate(app, data_set, public=public_version)
    cm_layout = comparison_maker.initiate(app, data_set, public=public_version)



    landing_page = Div(style={"background-color":our_colour}, children=[
        Div(style={'margin':'auto', 'width':'980px', 'height':'450px'}, children=[
            html.H1(
                'DDRcs – DDR CRISPR screen data explorer',
                className='home-text',
            ),
            html.Br(),
            html.P(
                "A web portal for the exploration of DNA damage reponse (DDR) related CRISPR screens. "
                "DDRcs allows researchers to see results for a few selected genes across all the screens "
                "in the database, explore all the results of a specific screen, and directly compare treatments "
                "to find differential effects.",
                className='home-text',
            ),
            html.Br(),
            # documentation link
            Div(className='centre',
                children=[dbc.Button("View documentation", href='https://docs.google.com/document/d/1RnDvb7NFjlNH52lqPAI5QSRFPWl-xSA58p2XmMM_P9I/edit?usp=sharing',
                       color="light", className="lg-1", size='lg', target="_blank")]),
            html.H2(
                "About the analyses",
                className='home-text',
            ),
            html.P(
                "Screen results presented here come from publically available, or personally provided, "
                "count data analysed using DrugZ and MAGeCK. A pseudocount of 5 was added to all counts "
                "prior to analysis, otherwise they were not adjusted. "
                "Where the experimental design used biological clones, paired analysis "
                "modes were used when appropriate. In the future, results from other analysis methods "
                "will be included, and preprocessing of the data will be investigated.",
                className='home-text',
            )

        ]),
    ])


    # header with links
    header = html.Header(className='myheader', children=[
        html.A(href="/",
            children=html.Img(src="assets/images/DDRCS_LOGO_No_Background.png", alt="SPJ Logo",
                 width='190px')
        ),
        html.Nav(className='navbar', children=[
            html.A(href='/gene-explorer', children="Query Genes"),
            html.A(href='/screen-explorer', children="Explore Screens"),
            html.A(href='/comparison-explorer', children="Compare treatments"),
            html.A(href='/about', children="About"),
        ])
    ])

    # main div that gets updated with content
    contents = Div([], id='graph-div', className='graphbox')

    # logos and external links
    footer = Div(className='myfooter', children=[
        Div(
            html.P("Developed by John C. Thomas. Data processing by John C. Thomas, Vipul Gupta & Simon Lam."),
            className='center-flex',
        ),

        Div(
            className='center-flex',
            children=[
                html.A( rel='noopener nofollow', target='_blank', href='https://www.cam.ac.uk/',
                   children=[html.Img(src="assets/images/univ_cam_logo.jpg",)]
                ),
                html.A(rel='noopener nofollow', target='_blank', href='https://www.gurdon.cam.ac.uk/', children=[
                    html.Img(src="assets/images/gurdon_logo.jpg", )
                ]),
                html.A(rel='noopener nofollow', target='_blank', href='https://www.stevejacksonlab.org/', children=[
                    html.Img(src="assets/images/new_lab_logo.png", )
                ])
            ]
        )
    ])

    app.layout = Div([
        dcc.Location(id='url', refresh=False),
        header,
        contents,
        footer,
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
    data_set, port, dash_debug, debug_messages, public = args = parse_clargs()
    print(args)
    import logging
    if debug_messages:
        print('Debug messages on')
        werklog = logging.getLogger('werkzeug')
        werklog.setLevel(logging.ERROR)
        LOG.setLevel('DEBUG')

    app = initiate_app(data_set, public)
    app.run_server(debug=dash_debug, host='0.0.0.0', port=port)




