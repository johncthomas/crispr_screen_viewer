
#!/usr/bin/env python
import logging
import sys
import copy
from typing import List, Dict, Tuple, Collection

import pandas as pd
import numpy as np

from dash.exceptions import PreventUpdate
from dash import (
    dash,
    dcc,
    html,
    dash_table,
    callback_context
)

Div = html.Div
import plotly.graph_objs as go

import pathlib, os
from dash.dependencies import Input, Output, State
import typing

from crispr_screen_viewer.functions_etc import (
    DataSet,
    LOG,
    style_gene_selector_div,
    get_cmdline_options,
    html_small_span,
    launch_page,
)
from crispr_screen_viewer.shared_components import (
    get_lab_val,
    get_gene_dropdown_lab_val,
    get_annotation_dicts,
    register_gene_selection_processor,
    spawn_gene_dropdown,
    spawn_filter_dropdowns,
    select_color, view_color
)

from crispr_screen_viewer.selector_tables import (
    spawn_selector_tables,
    spawn_selector_tabs,
    spawn_treatment_reselector,
    get_selector_table_filter_keys,
    register_exptable_filters_comps,
)


import dash_bootstrap_components as dbc
from dash_bio import Clustergram

PAGE_ID = 'cg'


def initiate(app:dash.Dash, dataset, debug=True):
    score = dataset.get_score_fdr('drz', 'drz')['score']






if __name__ == '__main__':

    #launch_page(*get_cmdline_options())
    launch_page('/Users/johnc.thomas/Dropbox/crispr/DDRcs/app_data/toy2',
                8050, True, 'Clustergram', initiate)