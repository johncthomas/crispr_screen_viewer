#!/usr/bin/env python
import pathlib, os
#import plotly.express as px
from argparse import ArgumentParser
from crispr_screen_viewer import screen_explorer
#from dashapp import app as application
from crispr_screen_viewer.functions_etc import launcher_parser


if __name__ == '__main__':
    parser = ArgumentParser(
        parents=[launcher_parser],
        description="Dash app for exploring screen data.",
        add_help=True,
    )
    args = parser.parse_args()

    screen_explorer.launch(os.path.join('app_data', args.data_version), args.port, args.debug)


