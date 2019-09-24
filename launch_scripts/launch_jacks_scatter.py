#!/usr/bin/env python
from crispr_screen_viewer import scatter
from crispr_screen_viewer.util import tabulate_score

import yaml

import sys
import os

if len(sys.argv) < 3 or sys.argv[1] in ('-h', '--help'):
    print('launch_jacks_scatter.py EXPD PORT\n'
          'Assumes the script is ran from the dir above the experiment dir.')
          #'RESULTS_ROOT is the dir above the one specified by expd["analysis_name"] and is the EXPD dir by default.')
    exit(0)

# if we assume the expd is in the screens dir then we just need the that
expd_fn = sys.argv[1]
# todo: optionally take a second argv giving screen dir
# screens_dir = os.path.dirname(expd_fn)

# for running on ia1
expd = yaml.safe_load(open(expd_fn))
port = int(sys.argv[2])
if len(sys.argv) > 3:
    dist_lim = float(sys.argv[3])
else:
    dist_lim = 0.3


f = f"{expd['exp_name']}/{expd['analysis_name']}/jacks_median/files/{expd['file_prefix']}."
_tables = {}
for _ctrl_grp in expd['controls'].keys():
    _tables[_ctrl_grp] = tabulate_score(f + _ctrl_grp + '.')
# tab = tabulate_score(f)
app = scatter.spawn_scatter(_tables, 'jacks', expd, distance_filter=dist_lim)
server = app.server
app.run_server(debug=True, host='0.0.0.0',  port=port)