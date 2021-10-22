#!/usr/bin/env python

from argparse import ArgumentParser
from crispr_screen_viewer import multiscreen_gene_viewer, screen_explorer, comparison_maker

if __name__ == '__main__':
    from argparse import ArgumentParser
    launcher_parser = ArgumentParser(add_help=False)

    launcher_parser.add_argument(
        'VIEWER',
        help="""Name of the viewer to launch. Options are:
        msgv (Multiscreen Gene Viewer), se (Screen Explorer), cm (Comparison Maker)""",
    )

    launcher_parser.add_argument(
        '-p', '--port', metavar='PORT',
        help='Port used to serve the charts',
        required=True,
    )

    launcher_parser.add_argument(
        '-d', '--data-path',
        dest='data_path',
        help="Name of the directory that contains the data from screens.",
        required=True,
    )
    launcher_parser.add_argument(
        '--debug', action='store_true',
        help='Launch app in debug mode'
    )

    parser = ArgumentParser(parents=[launcher_parser],
                        description="Dash app for exploring screen data.",
                        add_help=True,)


    args = parser.parse_args()


    viewer = {'msgv':multiscreen_gene_viewer,
              'cm':comparison_maker,
              'se':screen_explorer}[args.VIEWER]

    viewer.launch(args.data_path, args.port, args.debug)