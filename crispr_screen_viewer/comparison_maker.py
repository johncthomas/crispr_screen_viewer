#!/usr/bin/env python
import logging
import sys
import copy
from typing import List, Dict, Tuple, Collection

import pandas as pd
import numpy as np

from dash.exceptions import PreventUpdate
from dash import dash, dcc, html, dash_table, callback_context
Div = html.Div
import plotly.graph_objs as go

import pathlib, os
from dash.dependencies import Input, Output, State
import typing

from crispr_screen_viewer.functions_etc import (
    LOG,
    style_gene_selector_div,
    get_cmdline_options,
    html_small_span,
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

from crispr_screen_viewer.dataset import DataSet

import dash_bootstrap_components as dbc

def initiate(app, data_set:DataSet, public):
    """Source directory should contain the relevant info: metadata.csv,
    screen_analyses and expyaml directories."""

    # this is hardcoded in at least selector_tables.py, filter initialisation
    PAGE_ID = 'cm'

    comparisons = data_set.comparisons
    experiments_metadata = data_set.experiments_metadata
    print(experiments_metadata)

    ###############
    # COMPONENTS

    # ** filter_dropdowns **
    # Used to filterdown the rows of the experiment and comparisons tables.
    filter_keys = get_selector_table_filter_keys(public=public, )

    filter_dropdowns = {tabk:spawn_filter_dropdowns(PAGE_ID, tabk, comptabk, comparisons)
                        for tabk, comptabk in filter_keys.items()}

    selctr_tables = spawn_selector_tables(
        app, data_set, filter_keys, public_version=public, id_prefix=PAGE_ID,
    )

    scatter_chart = Div([dcc.Graph(
        id='graph',
        config={
            # 'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d',
            #                            'autoScale2d', 'resetScale2d'],
            'editable':True,
            'edits':{'annotationPosition':False},
        },
        style={'height': '1000px', 'width':'1000px', 'padding':'0px'},
        figure={'layout':{'clickmode':'event+select', 'dragmode':'select'}}
    )])

    delta_table = dash_table.DataTable(
        id='table',
        columns=[{'name':c, 'id':c} for c in ['Gene', "X normZ", 'Y normZ', "∆NormZ (Y-X)"]],
        sort_action='native',
        style_cell={
            'height': 'auto',
            # all three widths are needed
            'minWidth': '180px', 'width': '180px', 'maxWidth': '180px',
            'whiteSpace': 'normal'
        },
        selected_rows=[],
        #row_selectable='single',
    )

    selected_store = dcc.Store('selected-comps',
                    data={'X':None, 'Y':None})
    previous_row_store = dcc.Store('previous-rows',
                              data=[])
    comp_table_data = dcc.Store('cm-comp-table-data', data=None)

    exptab, comptab = spawn_selector_tabs(
        app,
        PAGE_ID,
        filter_dropdowns,
        selctr_tables,
        exp_tab_text=('Filter treatments shown in the Select Treatments tab to those from specific articles. '
                      'Select Treatments shows all possible by default.'),
        comp_tab_text=('Choose treatments to be displayed on scatter plot using the '
                       'radio buttons in the X/Y columns.'),
    )

    ################
    # LAYOUT
    tabs = dcc.Tabs(
        id=f'{PAGE_ID}-tabs',
        value=f'{PAGE_ID}-comp-tab',
        className='selector-results-tabs',
        children=[
            exptab,
            comptab,
            # results tabs)
            dcc.Tab(
                label='Chart',
                value='chart-tab',
                className='data-tab', selected_className='data-tab--selected',
                children=[
                    Div([scatter_chart]),
                ]
            ),
            dcc.Tab(
                label='Table',
                value='table-tab',
                className='data-tab', selected_className='data-tab--selected',
                children=Div([delta_table], 'table-div', style={'width':'1000px'})
            )
        ]
    )

    layout = Div([
        html.H1('Compare treatments'),
        html.P('Choose two sets of results and display the gene NormZ scores '
               'in a biplot. When comparing results from different experiments, be aware'
               ' that differences between X and Y may result from '
               'relevant biological effects of the treatments; or from differences in '
               'the library performance, cell lines, or technical details of the way '
               'the screens were performed.',
               className='explanatory-paragraph'),
        selected_store,
        previous_row_store,
        comp_table_data,
        Div([
            Div([tabs, ], style={'display':'inline-block',}),
            Div([
                spawn_treatment_reselector(PAGE_ID, is_xy=True)
            ], id=f'{PAGE_ID}-results-control-panel')
        ],  className='tab-content-box'),
        Div(spawn_gene_dropdown(app, PAGE_ID))
    ])


    ################
    # CALLBACKS
    # when an exp and timepoint is selected, update the X/Y options
    def get_xyscores_genes(xk, yk, selected_analysis_type='drz'):
        """return score series from Dataset with unified indexes"""
        score_fdr = data_set.get_score_fdr(selected_analysis_type)

        x, y = [score_fdr['score'][k].dropna() for k in (xk, yk)]
        shared_genes = x.index.intersection(y.index)
        x, y = [xy[shared_genes] for xy in (x,y)]
        return x, y


    # update the keys stored in the data store, and gene selection options
    # when the x/y comparisons are changed.
    @app.callback(
        Output('selected-comps', 'data'),
        Output(f'{PAGE_ID}-gene-dropdown', 'options'),
        Input(f'{PAGE_ID}-x-selector', 'value'),
        Input(f'{PAGE_ID}-y-selector', 'value'),
    )
    def update_selection(xk, yk):
        LOG.debug(f'update_selection({xk}, {yk})')
        selected_comps = {'X':xk, 'Y':yk}
        if (not xk) or (not yk):
            LOG.debug('\tOnly 0-1 comps selected, not making graph')
            return selected_comps, dash.no_update
        x, _ = get_xyscores_genes(xk, yk, 'drz')
        return selected_comps, get_gene_dropdown_lab_val(data_set, x.index)

    # enable selecting genes by interacting with the graph
    register_gene_selection_processor(
        app, 'graph',
        State(f'{PAGE_ID}-gene-dropdown', 'value'),
        Output(f'{PAGE_ID}-gene-dropdown', 'value'),
    )


    # update the chart and table
    @app.callback(
        Output('graph', 'figure'),
        Output('table-div', 'children'),
        #Output('missing-analysis', 'children'),
        Input('selected-comps', 'data'),
        Input(f'{PAGE_ID}-gene-dropdown', 'value'),

    )
    def update_chart_table(selected_comps, selected_genes):
        xk, yk = selected_comps['X'], selected_comps['Y']

        LOG.debug(f'update_chart_table(selected_comps={selected_comps},\n'
                  f'                   selected_genes={selected_genes})')
        LOG.debug(str(callback_context.triggered))
        if (not xk) or (not yk):
            LOG.debug('not updating')
            raise PreventUpdate

        x, y = get_xyscores_genes(xk, yk)
        delta = y-x

        # **FIGURE**
        # get a title
        treatlines = []
        # write a gramatical line for each treatment
        for k in (xk, yk):
            row = comparisons.loc[k]
            if '-KO' not in row['Treatment']:
                if row['KO'] == 'WT':
                    ko = ''
                else:
                    ko = f" {row['KO']}-KO"
            else:
                ko = ''
            treatlines.append(f"<b>{row['Treatment']} in {row['Cell line']}{ko} cells</b>")
        # put them together
        rowx = comparisons.loc[xk]
        rowy = comparisons.loc[yk]
        metadeets_str = "{}, {} library"
        xymet = [metadeets_str.format(row['Citation'], row['Library']) for row in (rowx, rowy)]

        if xymet[0] == xymet[1]:

            xy_metadeets = f'{xymet[0]}<br>{html_small_span(f"(IDs: {xk} & {yk})")}<br>'
        else:

            xy_metadeets = (

                f'X: {xymet[0]}  {html_small_span(f"(ID: {xk})")}<br>'
                f'Y: {xymet[1]}  {html_small_span(f"(ID: {yk})")}<br>'
            )
        title = (f"{treatlines[0]} Vs {treatlines[1]}<br>"+xy_metadeets
                )
        LOG.debug(title)

        # axis labels
        # xy_labs = [f"{comparisons.loc[k, 'Ctrl samp']} ➤ {comparisons.loc[k, 'Treat samp']}"
        #            for k in (xk, yk)]
        xy_labs = xymet

        fig = go.Figure(
            data=go.Scattergl(
                x=x,
                y=y,
                customdata=delta,
                mode='markers',
                text=x.index,
                hovertemplate= (
                        "<b>%{text}</b><br>" +
                        "∆normZ: %{customdata:.2f}<br>" +
                        f"{xy_labs[0]}: "+"%{x:.2f}<br>" +
                        f"{xy_labs[1]}: "+"%{y:.2f}<br>"
                )
            ),
            layout={'clickmode':'event+select',
                    'dragmode':'select'},
        )

        expander = 0.1

        xymax = max(x.max(), y.max())
        xymin = min(x.min(), y.min())
        pad = abs(xymax-xymin)*expander
        xymax += pad
        xymin -= pad
        fig.update_xaxes(range=[xymin, xymax])
        fig.update_yaxes(range=[xymin, xymax])

        fig.update_layout(
            title=title,
            xaxis_title=xy_labs[0],
            yaxis_title=xy_labs[1],
        )

        # annotations
        new_annotations = get_annotation_dicts(x[selected_genes], y[selected_genes], selected_genes)
        for anot in new_annotations:
            fig.add_annotation(
                **anot
            )

        # **TABLE**
        # update the table if the comps have changeded
        if any([trig['prop_id'] == 'selected-comps.data'
                for trig in callback_context.triggered]):

            table_cols = ['Gene', xy_labs[0], xy_labs[1], "∆NormZ (Y-X)"]

            tab = pd.DataFrame(
                dict(zip(table_cols, [x.index, x, y, delta]))
            )
            records = tab.sort_values("∆NormZ (Y-X)").to_dict('records')
            delta_table = dash_table.DataTable(
                id='table',
                columns=[
                    {'name':c, 'id':c, 'type':'numeric', 'format': {'specifier': '.2f'},}
                    for c in table_cols
                ],
                sort_action='native',
                filter_action='native',
                data = records,
                selected_rows=[],
                row_selectable='multi',
                export_format='csv',
            )
            table_output = [delta_table]
        else:
            table_output = dash.no_update

        return fig, table_output

    return layout


if __name__ == '__main__':
    from crispr_screen_viewer.functions_etc import launch_page
    launch_page(*get_cmdline_options(), 'Comparisons', initiate)