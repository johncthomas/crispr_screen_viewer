from dash.dependencies import Input, Output, State
from dash.dash_table import DataTable
from dash.dash_table.Format import Format, Scheme

from dash import dcc, html, callback_context
from crispr_screen_viewer.functions_etc import (
    style_gene_selector_div,
    datatable_column_dict
)
from crispr_screen_viewer.dataset import ANALYSESTYPES
from dash.exceptions import PreventUpdate
Div = html.Div
import pandas as pd

import dash_bootstrap_components as dbc

#from crispr_screen_viewer.functions_etc import DataSet

from typing import Tuple, List, Dict, Callable

from loguru import logger

#from https://sashamaps.net/docs/resources/20-colors/, 99%,

colours = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
           '#42d4f4', '#f032e6', '#fabed4', '#469990',
           '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3',
           '#000075', '#a9a9a9', '#333333', ]

select_color = '#E865DF'
view_color   = '#DFE865'

get_lab_val = lambda arr: [{'label': v, 'value':v} for v in arr]

def get_gene_dropdown_lab_val(data_set, genes):
    return [{'label': data_set.dropdown_gene_label(gn), 'value': gn} for gn in genes]

styles = {'selector':{'display': 'inline-block', 'border-collapse': 'separate',
                      'border-spacing': '15px 50px', 'margin-bottom': '15px', 'width':'150px'},
          'hidden':{'display':'none'}}
big_text_style = {
    'font-family': 'Arial, Helvetica, sans-serif',
    'font-size': '25px',
    'letter-spacing': '-0.4px',
    'word-spacing': '-0.4px',
    'font-weight': '700'}



def get_annotation_dicts(xs,ys,txts, annote_kw=None) -> List[dict]:
    """dicts defining annotations with x/y/text values as given.

    annote_kw used to update/override formatting values of annotations.
    See """
    annotations = []
    if annote_kw is None:
        annote_kw = {}
    for x,y,txt in zip(xs,ys,txts):
        d = dict(
            x=x,
            y=y,
            xref='x',
            yref='y',
            text=txt,
            showarrow=True,
            arrowhead=1,
            arrowwidth=0.4,
            bgcolor='rgba(255,255,255,0.55)',
        )
        # in-place method...
        d.update(annote_kw)
        annotations.append(d)
    return annotations



def get_stat_source_selector(id_prefix, cardheader='Select analysis algorithm') -> dbc.Card:
    """List of single Div with dcc.RadioItems with id
    '{id_prefix}-stat-source-selector'. Options from `options_analyses`"""

    sigsourceid = f'{id_prefix}-stat-source-selector'

    return dbc.Card(
        [
            dbc.CardHeader(
                cardheader
            ),
            dbc.CardBody(
                children=[
                    dcc.RadioItems(
                        id=sigsourceid,
                        options=[{'label':ans.label, 'value':ans.name} for ans in ANALYSESTYPES],
                        value=ANALYSESTYPES.default.name,
                    )
                ]
            ),
        ], style={'width': '170px'}
    )

# **DATATABLE**
def create_datatable(data_df=None, columns_if_no_df=None):

    #numfmt = Format(precision=3, scheme=Scheme.decimal_or_exponent)
    #formatted = Format()
    #numfmt = formatted.precision(3)

    # None of this works and I need to figure out why
    #
    # conditional formatting for the FDR
    thresholds = [0.6, 0.3, 0.1, -0.1]
    fdr_colours = ['#ff0000', '#ff9933', '#ffff00', '#66ff33']
    prev_fdr = 1
    cond_fmts = []
    for fdr, clr in zip(thresholds, fdr_colours):
        fmt = {
            'if':{
                'filter_query': f'{prev_fdr} <= {{FDR}} < {fdr}',
                'column_id':'FDR'
            },
            'backgroundColor':clr
        }
        cond_fmts.append(fmt)
        prev_fdr = fdr

    # just a place holder to start with
    if data_df is None:
        return DataTable(
            id='table',
            columns=[datatable_column_dict(c) for c in columns_if_no_df],
        )

    return DataTable(
        id='table',
        columns=[datatable_column_dict(c) for c in data_df.columns],
        data=data_df.to_dict('records'),
        export_format='csv',
        style_data_conditional=cond_fmts,
    )

def register_gene_selection_processor(app, figure_id, selected_genes_state, output) -> Callable:
    """Registers a function for putting genes selected via interacting with
    the graph into the gene dropdown. Registers a callback with following
    arguments:

    output, Input(figure_id, 'selectedData'), state

    output and state pointing to the same object probably makes the most sense,
    but using Output/State
    """

    # update selected genes by points that are selected on the graph
    # this should only ever add points.
    @app.callback(
        output,
        Input(figure_id, 'selectedData'),
        selected_genes_state,
    )
    def put_selected_genes_into_dropdown(selected_data, dropdown_genes):
        logger.debug(f'Adding genes from {figure_id}.selectedData: {selected_data}')

        if not selected_data:
            raise PreventUpdate

        selected_genes = set()
        for p in selected_data['points']:
            selected_genes.add(p['text'])

        if selected_genes.issubset(dropdown_genes):
            logger.debug('...no new genes, preventing update')
            raise PreventUpdate

        return dropdown_genes+list(selected_genes.difference(dropdown_genes))

    return put_selected_genes_into_dropdown


def spawn_gene_dropdown(app, id_prefix) -> Div:
    """Returns layout for dropdown,
    register_gene_selection_processor(app, fig_id).

    Dropdowns need extra callbacks to communicate with each other. That
    is not handled here."""

    return Div(
        id=f'{id_prefix}-gene-selector-div',
        style=style_gene_selector_div,
        children=[
            html.Label('Label genes:',
                       htmlFor=f'{id_prefix}-gene-dropdown',
                       style=big_text_style),
            dcc.Dropdown(
                id=f'{id_prefix}-gene-dropdown',
                placeholder='Label genes',
                style={'width':'100%', 'height':'100px', },
                value=[],
                options=[],
                clearable=True,
                multi=True,
            )
        ]
    )

def spawn_filter_dropdowns(
        id_prefix, table_str, filter_cols:List[str], comparisons, values:dict=None,
        card_header='Filter table rows by their contents.',
) -> List[dbc.Card]:
    """Get Dropdown objects for layout. No callbacks registered.

    filter_cols: list giving the columns for which filter boxes will be
        spawned.

    values:

    Return list of div with dcc.Dropdowns with
    id=f"{id_prefix}-{table_str}-filter-{col}" for col in filter_cols.
    Options from comparisons[col].

    """

    logger.debug(f"{comparisons=}, {values=}")

    filter_dropdowns = []
    for col in filter_cols:
        idd = f"{id_prefix}-{table_str}-filter-{col}"
        logger.debug(f"Filter dropdown: register ID={idd}")

        try:
            value = values[col]
        except:
            value = []

        drpdwn = dbc.Card(
            [
                dbc.CardBody([
                    dcc.Dropdown(
                        id=idd,
                        placeholder='Filter by '+col,
                        multi=True,
                        style={'min-height':'80px', 'width':'100%'},
                        value=value,
                        options=[{'label':v, 'value':v} for v in sorted(comparisons[col].unique())]
                    ),
                ])
            ],
            className='mt-3'
        )
        filter_dropdowns.append(
            dbc.Tab(
                drpdwn,
                label=col,
                #className='filter-tab'
            )
        )

    filter_tabs = dbc.Card(
        [
            dbc.CardHeader(card_header),
            dbc.CardBody(dbc.Tabs(children=filter_dropdowns, active_tab=None)),
        ]
    )

    return [filter_tabs]