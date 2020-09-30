#!/usr/bin/env python
#todo by default load all available results.

from crispr_screen_viewer import volcano
import crispr_screen_viewer
import pandas as pd
np = pd.np
import os
from pathlib import Path
from crispr_screen_viewer.util import tabulate_mageck, tabulate_score
from argparse import ArgumentParser
import yaml
import os, sys

#todo allow specification of specific csv instead of results dir
# to read multiindex: index_col=0, header=[0,1]


analyses = {
    'jacks_median':{'tabulate':tabulate_score,
             'label':'JACKS',
             'columns':['jacks_score', 'fdr_log10']},
    'mageck':{'tabulate':tabulate_mageck,
              'label':'MAGeCK',
              'columns':['lfc', 'fdr_log10']},
    'jacks':{'tabulate':tabulate_score,
             'label':'JACKS',
             'columns':['jacks_score', 'fdr_log10']},

}

def launch(args):

    port = args.port
    y_filter = args.y_filter
    tabz = {}
    shared_cols = ['LFC or score', 'fdr_log10']

    for expd_fn in args.expd_yaml:
        #expd_fn = args.expd_yaml

        expd = yaml.safe_load(open(expd_fn))

        if args.groupings:
            groupings = pd.read_csv(args.groupings, index_col=0).iloc[:, 0]
        else:
            groupings = None

        expname = expd['exp_name']
        anlname = expd['analysis_name']
        prefix  = expd['file_prefix']
        ctrl_groups = expd['controls'].keys()
        # all tables should have the same columns


        for group in ctrl_groups:
            # load each results set and set the column names

            for analysis in args.analyses.split(','):

                analysis_dict = analyses[analysis]
                tab = analysis_dict['tabulate'](f'{expname}/{anlname}/{analysis}/files/{prefix}.' + group + '.')
                tab = tab.reindex(columns=analyses[analysis]['columns'], level=1)
                tab.columns.set_levels(shared_cols, 1, inplace=True)

                #
                tabz[f"{anlname} - {group} ({analysis_dict['label']})"] = tab.copy()
    print(crispr_screen_viewer.__file__)
    app = volcano.spawn_volcanoes(tabz, shared_cols, filterYLessThan=y_filter,)
    #groupings=groupings
    app.run_server(host='0.0.0.0', port=port)

if __name__ == '__main__':
    # expd_fn, port,
    # fdr filter, analysis types
    parser = ArgumentParser(description='Start serving vocano plots for multiple analysis types.')
    parser.add_argument(
        'expd_yaml', metavar='EXPD_PATH', nargs='+',
        help='Location of expd.yml containing sample replicate and control info.'
    )
    parser.add_argument(
        '-p', '--port', metavar='PORT',
        help='Port used to serve the charts'
    )
    parser.add_argument(
        '--y_filter', metavar='MIN_Y', type=float, default=0,
        help='To be included, a gene must have a Y value (probably -log10(fdr) ) in at least one comparison'
        ' to be included in the chart. Defaults to no filter'
    )
    parser.add_argument(
        '-a', '--analyses', metavar='LIST', default='jacks,jacks_median,mageck',
        help='Comma sep string (no spaces) of analysis types to include, currently supported: '+', '.join(analyses.keys())
    )
    parser.add_argument(
        '--groupings', metavar='PathToCSV', default=None,
        help='An optional CSV that indicates groups that will be plotted as different colours on the chart. '
        'Note that points ommited from the CSV, or with no specified group will not be plotted. '
        '(Header is expected, but ignored).'
    )

    launch(parser.parse_args())