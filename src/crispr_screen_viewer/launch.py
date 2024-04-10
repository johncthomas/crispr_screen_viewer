#!/usr/bin/env python
import os.path
from argparse import ArgumentParser
from urllib.parse import urlparse
from importlib import resources

import sqlalchemy

from crispr_screen_viewer import (
    multiscreen_gene_viewer,
    screen_explorer,
    comparison_maker,
    legal
)
from crispr_screen_viewer.functions_etc import set_loguru_level
from crispr_screen_viewer.dataset import (
    DataSet,
    load_dataset
)

from crispr_screen_viewer.shared_components import (
    logger
)
from dash import dash, html, dcc, Input, Output, State
Div = html.Div
import flask
import dash_bootstrap_components as dbc

our_colour = '#3db9c2'

def get_home_page_text():
    landing_page = Div(style={"background-color": our_colour}, children=[
        Div(style={'margin': 'auto', 'width': '980px', },
            children=[
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
                    children=[
                        dbc.Button(
                            "View documentation",
                            href='https://docs.google.com/document/d/1RnDvb7NFjlNH52lqPAI5QSRFPWl-xSA58p2XmMM_P9I/edit?usp=sharing',
                            color="light", className="lg-1", size='lg', target="_blank")
                    ]
                    ),
                html.Br(),
                html.P(
                    [
                        "Please send any comments, feedback, or issues to ",
                        html.A("sl681@cam.ac.uk", className='home-text', href="mailto:sl681@cam.ac.uk"),
                        "."
                    ],
                    className='home-text'
                ),
                html.Br(),
                html.H2(
                    "About the analyses",
                    className='home-text',
                ),
                html.P(
                    "Screen results presented here come from publicly available, or personally provided, "
                    "count data analysed using DrugZ and MAGeCK. A pseudocount of 5 was added to all counts "
                    "prior to analysis, otherwise they were not adjusted. "
                    "Where the experimental design used biological clones, paired analysis "
                    "modes were used when appropriate. In the future, results from other analysis methods "
                    "will be included, and preprocessing of the data will be investigated.",
                    className='home-text',
                ),
                html.Br(),
                html.H2(
                    "Citation",
                    className='home-text'
                ),
                html.P(
                    [
                        "If you use DDRcs in your work, please cite: Awwad, S.W., Serrano-Benitez, A., Thomas, J.C. et al. Revolutionizing DNA repair research and cancer therapy with CRISPR–Cas screens. ",
                        html.I("Nat Rev Mol Cell Biol"),
                        " (2023). ",
                        html.A(
                            "https://doi.org/10.1038/s41580-022-00571-x", className='home-text', href="https://doi.org/10.1038/s41580-022-00571-x"
                        ),
                        (".")
                    ],
                    className='home-text',
                ),
                html.Br(),
            ]
            ),
    ])
    return landing_page


def get_header():
    return html.Header(className='myheader', children=[
        html.A(href='home',
               children=html.Img(src="assets/images/DDRCS_LOGO_No_Background.png", alt="SPJ Logo",
                                 width='190px')
               ),
        html.Nav(className='navbar', children=[
            html.A(href='gene-explorer', children="Query Genes"),
            html.A(href='screen-explorer', children="Explore Screens"),
            html.A(href='comparison-explorer', children="Compare treatments"),
            html.A(href='legal', children="Legal notices"),
            # html.A(href='/about', children="About"),
        ])
    ])


def get_footer():
    return Div(className='myfooter', children=[
        Div(
            html.P("Developed by John C. Thomas & Simon Lam. Data processing by John C. Thomas, Vipul Gupta, Simon Lam, & Tuan-Anh Tran."),
            className='center-flex',
        ),

        Div(
            html.A("https://github.com/johncthomas/crispr_screen_viewer",
                   href="https://github.com/johncthomas/crispr_screen_viewer"),
            className='center-flex',
        ),

        html.Br(),

        Div(
            className='spacey-flex',
            children=[
                html.A(rel='noopener nofollow', target='_blank', href='https://www.cam.ac.uk/',
                    children=[
                           html.Img(src="assets/images/univ_cam_logo.jpg", )
                    ]),
                html.A(rel='noopener nofollow', target='_blank', href='https://www.stevejacksonlab.org/', children=[
                    html.Img(src="assets/images/new_lab_logo.png", )
                ]),
                html.A(rel='noopener nofollow', target='_blank', href='https://www.cruk.cam.ac.uk/', children=[
                    html.Img(src="assets/images/cruk_cambridge_i_pos_cmyk.jpg", )
                ]),
            ]
        ),

    ])


def create_layout(app, data_set:DataSet, public_version=False, url_base='/') \
        -> None:
    """Get and apply layout for all pages of the app, register page changing
    callback."""

    # main div that gets updated with content
    contents = Div([], id='graph-div', className='graphbox')

    # register the callbacks and get the page layouts
    msgv_layout = multiscreen_gene_viewer.initiate(app, data_set, public=public_version)
    se_layout = screen_explorer.initiate(app, data_set, public=public_version)
    cm_layout = comparison_maker.initiate(app, data_set, public=public_version)
    legal_layout = legal.initiate()

    landing_page = get_home_page_text()

    # header with links
    header = get_header()

    # logos and external links
    footer = get_footer()

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
        # I might add extra information to the end of this URL in the future
        #   this gives us just the part we want
        print(pathname)
        pathname = urlparse(pathname).path
        # remove the base pathname if there is one before doing the query
        pathname = pathname.replace(url_base, '')
        if pathname == 'gene-explorer':
            return msgv_layout
        elif pathname == 'screen-explorer':
            return se_layout
        elif pathname == 'comparison-explorer':
            return cm_layout
        elif pathname == 'legal':
            return legal_layout
        elif pathname == 'home':
            return landing_page
        else:
            return landing_page


def from_cli(args):
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
        help="Name of the directory that contains the screens' data.",
        required=True, default=None,
    )
    launcher_parser.add_argument(
        '--app-debug', action='store_true',
        help='Launch Dash app in debug mode, tends to break things, but allows you to look at the'
             ' callback graph and JS debug messages.'
    )
    launcher_parser.add_argument(
        '--debug-messages', action='store_true', default=False,
        help='Set log level to debug – print messages describing the internal state of the app. '
             'Also hide Werkzeug messages'
    )
    launcher_parser.add_argument(
        '--public-version', action='store_true',
        help="Don't hide the data-source and analysis-type selectors."
             " In the future analysis-type might have its own option."
    )
    launcher_parser.add_argument(
        '--url-pathname', default="/",
        help="URL base pathname. Needs to end in a /."
    )
    launcher_parser.add_argument(
        '-u', '--database-url',
        help='A SQL database URL, as understood by SQL Alchemy. '
             '\nSee: https://docs.sqlalchemy.org/en/20/core/engines.html#database-urls.'
             '\nE.g.: sqlite:////home/user/some/dir/databasefile.db',
        default=None, required=False
    )

    parser = ArgumentParser(parents=[launcher_parser],
                            description="Dash app for exploring CRISPR screen data.",
                            add_help=True,)

    args = parser.parse_args(args)

    app = init_app(
        data_path=args.data_path,
        database_url=args.database_url,
        debug_messages=args.debug_messages,
        public_version=args.public_version,
        url_base=args.url_pathname
    )

    app.run_server(debug=args.app_debug, host='0.0.0.0', port=args.port)


def init_app(
        data_path:str=None, database_url=None,
        debug_messages=False,
        public_version=None, url_base='/',
        app_title='CRISPR screen viewer',
        db_echo=False,
):

    if data_path is None:
        data_path = resources.files("crispr_screen_viewer").joinpath("data").__str__()

    if database_url is None:
        # Use absolute path cus
        database_url = f"sqlite:///{os.path.abspath(data_path)}/database.db"
        logger.info(f'Using datbase URL "{database_url}"')

    db_engine = sqlalchemy.create_engine(database_url, echo=db_echo)


    data_set = load_dataset(data_path, db_engine=db_engine)

    import logging
    if debug_messages:
        print('Debug messages on')
        werklog = logging.getLogger('werkzeug')
        werklog.setLevel(logging.ERROR)
        set_loguru_level(logger, 'DEBUG')
    else:
        set_loguru_level(logger, 'INFO')

    server = flask.Flask(__name__)

    app = dash.Dash(__name__,  server=server,
                    url_base_pathname=url_base, )

    create_layout(app, data_set, public_version=public_version, url_base=url_base)

    app.title = app_title

    return app

def get_server(**kwargs):
    """Calls init_app and returns the app.server object, without calling
    app.run_server

    This function exists to be a handle for Gunicorn."""

    return init_app(**kwargs).server

if __name__ == '__main__':
    import sys
    from_cli(sys.argv[1:])




