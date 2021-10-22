#!/usr/bin/env python

from argparse import ArgumentParser
from crispr_screen_viewer import multiscreen_gene_viewer
from crispr_screen_viewer.functions_etc import launcher_parser


if __name__ == '__main__':

    parser = ArgumentParser(parents=[launcher_parser],
                        description="Dash app for exploring screen data.",
                        add_help=True,)

    args = parser.parse_args()

    multiscreen_gene_viewer.launch(args.data_path, args.port, args.debug)