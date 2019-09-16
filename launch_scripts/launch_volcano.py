
from crispr_screen_viewer import volcano
import pandas as pd
np = pd.np
import os
from pathlib import Path
from crispr_screen_viewer.util import tabulate_mageck, tabulate_score

import yaml
import os, sys

tabz = {}

expd_fn = sys.argv[1]
port = int(sys.argv[2])

expd = yaml.safe_load(open(expd_fn))


expname = expd['exp_name']
anlname = expd['analysis_name']
prefix  = expd['file_prefix']
ctrl_groups = expd['controls'].keys()

# all tables should have the same columns
shared_cols = ['LFC or score', 'fdr_log10']

for group in ctrl_groups:
    # load each results set and set the column names
    magtab = tabulate_mageck(
        f'{expname}/{anlname}/mageck/files/{prefix}.' + group + '.',
    )
    magtab = magtab.reindex(columns=['lfc', 'fdr_log10'], level=1)
    magtab.columns.set_levels(shared_cols, 1, inplace=True)
    tabz[group+' (MAGeCK)'] = magtab

    jaktab = tabulate_score(
        f'{expname}/{anlname}/jacks_median/files/{prefix}.' + group + '.',
        return_ps=True
    )
    jaktab = jaktab.reindex(columns=['jacks_score', 'fdr_log10'], level=1)
    jaktab.columns.set_levels(shared_cols, 1, inplace=True)
    tabz[group+' (JACKS)'] = jaktab

app = volcano.spawn_volcanoes(tabz, shared_cols)
server = app.server
app.run_server(host='0.0.0.0', port=port)