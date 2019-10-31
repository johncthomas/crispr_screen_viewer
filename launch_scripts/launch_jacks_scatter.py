#!/usr/bin/env python
from crispr_screen_viewer import scatter
from crispr_screen_viewer.util import tabulate_score
from argparse import ArgumentParser
import yaml

import sys
import os

# if len(sys.argv) < 3 or sys.argv[1] in ('-h', '--help'):
#     print('launch_jacks_scatter.py EXPD PORT\n'
#           'Assumes the script is ran from the dir above the experiment dir.')
#           #'RESULTS_ROOT is the dir above the one specified by expd["analysis_name"] and is the EXPD dir by default.')
#     exit(0)


def launch(args):

    port = args.port

    _tables = {}

    #todo fix this hacky way of avoiding issues with scatter.valid_comparisons()
    # making a new control dict that has seperate group key for each analysis
    controls = {}

    for expd_fn in args.expd_yaml:
        #expd_fn = args.expd_yaml

        expd = yaml.safe_load(open(expd_fn))

        if not expd['labels']:
            expd['labels'] = list(expd['sample_reps'].keys())

        f = f"{expd['exp_name']}/{expd['analysis_name']}/jacks_median/files/{expd['file_prefix']}."

        for _ctrl_grp in expd['controls'].keys():
            nugroup = expd['analysis_name'] + ' ~ ' +_ctrl_grp
            _tables[nugroup] = tabulate_score(f + _ctrl_grp + '.')

            controls[nugroup] = expd['controls'][_ctrl_grp]

    out_expd = yaml.safe_load(open(args.expd_yaml[0]))
    out_expd['controls'] = controls

    print(_tables.keys())
    # use the first yamlas the expd, currently uses only labels and controls dicts,
    # should probably pass those things separately to make it more flexible.
    app = scatter.spawn_scatter(_tables, 'jacks', out_expd)
    app.run_server(debug=True, host='0.0.0.0',  port=port)


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


    launch(parser.parse_args())