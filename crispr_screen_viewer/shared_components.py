from dash.dependencies import Input, Output, State
import logging
from dash.dash_table import DataTable
from dash.dash_table.Format import Format, Scheme
from crispr_screen_viewer.functions_etc import datatable_column_dict
from dash import dcc, html, callback_context
from functions_etc import (
    cell_text_style,
    cell_number_style,
    LOG,
    timepoint_labels,
    style_gene_selector_div,
)
from dash.exceptions import PreventUpdate
Div = html.Div
import pandas as pd

import dash_bootstrap_components as dbc

#from crispr_screen_viewer.functions_etc import DataSet

from typing import Tuple, List, Dict, Callable

#from https://sashamaps.net/docs/resources/20-colors/, 99%,


colours = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990',
           '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#333333', ]
get_lab_val = lambda arr: [{'label': v, 'value':v} for v in arr]
styles = {'selector':{'display': 'inline-block', 'border-collapse': 'separate',
                      'border-spacing': '15px 50px', 'margin-bottom': '15px', 'width':'150px'},
          'hidden':{'display':'none'}}
big_text_style = {
    'font-family': 'Arial, Helvetica, sans-serif',
    'font-size': '25px',
    'letter-spacing': '-0.4px',
    'word-spacing': '-0.4px',
    'font-weight': '700'}

def get_treatment_label(row, analysis_label='') -> Tuple[str, str]:
    """Pass comparison row (either from data_set.comparisons.loc[compid] or
    from dashtable data), return a pair of strings.

    First string comparison specific, second line library, experiment ID."""
    if '-KO' not in row['Treatment']:
        if row['KO'] == 'WT':
            ko = ''
        else:
            ko = f" {row['KO']}-KO"
    else:
        ko = ''

    if analysis_label:
        analysis_label = f"{analysis_label}, "


    title = (f"Effect of {row['Treatment']} in {row['Cell line']}{ko} cells ({analysis_label}{row['Timepoint']})",
             f"{row['Library']} library, experiment ID {row['Experiment ID']}")

    return title

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



# DEPRECIATED
def get_reg_stat_selectors(app=None, id_prefix='') -> List[Div]:
    """Return radio selectors, for selecting stats to be used in plotting.

    Registers a function that handles the score selection. Needs to have
    output values prefix+"-score-selector" & prefix+"-fdr-selector" captured by the fig
    update callback.

    If app is None, you'll need to reimpliment all that."""

    if id_prefix:
        id_prefix = id_prefix+'-'

    if app is not None:
        @app.callback(
            [Output(id_prefix + 'score-selector', 'value'),
             Output(id_prefix + 'fdr-selector', 'value'),
             Output(id_prefix + 'mixed-div', 'style')],

            [Input(id_prefix + 'analysis-selector', 'value')],

            [State(id_prefix + 'score-selector', 'value'),
             State(id_prefix + 'fdr-selector', 'value')]
        )
        def select_stats_primary(selection, curr_score, curr_fdr):
            if selection == 'mixed':
                return curr_score, curr_fdr, styles['selector']
            else:
                return selection, selection, styles['hidden']

    return [
        Div([
            html.Label('Analysis type:  ', htmlFor='analysis-selector'),
            dcc.RadioItems(
                id=id_prefix + 'analysis-selector',
                options=[
                    {'label':'DrugZ', 'value':'drz'},
                    {'label':'MAGeCK',  'value':'mag'},
                    #{'label':'Mixed...', 'value':'mixed'}
                ],
                value='drz',
                labelStyle={'display': 'inline-block'},
            )
        ], style={**styles['selector'], **{'width':170}}),

        # This Div is hidden unless "Mixed" is chosen from Div above.
        #   Currently inaccessible
        Div([
            Div([
                html.Label('Effect size_____', htmlFor=id_prefix+'score-selector'),
                dcc.RadioItems(
                    id=id_prefix + 'score-selector',
                    options=[{'label':'NormZ', 'value':'drz'},
                             {'label':'LFC',  'value':'mag'}],
                    value='drz',
                ),
            ], style=styles['selector']),

            Div([
                html.Label('FDR source', htmlFor=id_prefix+'fdr-selector', ),
                dcc.RadioItems(
                    id=id_prefix + 'fdr-selector',
                    options=[
                        {'label':'DrugZ', 'value':'drz'},
                        {'label':'MAGeCK',  'value':'mag'}
                    ],
                    value='drz',
                )], style=styles['selector'])
        ],  id=id_prefix + 'mixed-div', style=styles['hidden'])
    ]

# def get_data_source_selector(data_set, id_prefix='') -> List[dcc.Checklist]:
#     """A Div with a checklist that will be populated data_sources and a
#     paragraph for reporting missing datasets
#
#     IDs: data-source-selector, missing-datasets"""
#
#     if id_prefix:
#         id_prefix = id_prefix+'-'
#
#     return [
#         html.Label('Select data sources:', htmlFor=id_prefix+'data-source-selector'),
#         dcc.Checklist(
#             id=id_prefix+'data-source-selector',
#             options=get_lab_val(data_set.data_sources),
#             value=data_set.data_sources, # set all selected by default
#             labelStyle={'display':'inline-block'}
#         ),
#         html.P([''], id=id_prefix+'missing-datasets'),
#     ]

options_analyses = [
    {'label':'DrugZ', 'value':'drz'},
    {'label':'MAGeCK', 'value':'mag'}
]


def get_stat_source_selector(id_prefix, label, cardheader='Select analysis algorithm') -> Div:
    """List of single Div with dcc.RadioItems with id
    '{id_prefix}-stat-source-selector'. Options from `options_analyses`"""

    sigsourceid = f'{id_prefix}-stat-source-selector'

    return Div(dbc.Card(
        [
            dbc.CardHeader(
                cardheader
            ),
            dbc.CardBody(
                children=[
                    html.Label(label, htmlFor=sigsourceid),
                    dcc.RadioItems(
                        id=sigsourceid,
                        options=options_analyses,
                        value='drz',
                    )
                ]
            ),
        ]
    ))

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
        LOG.debug(f'Adding genes from {figure_id}.selectedData: {selected_data}')

        if not selected_data:
            raise PreventUpdate

        selected_genes = set()
        for p in selected_data['points']:
            selected_genes.add(p['text'])

        if selected_genes.issubset(dropdown_genes):
            LOG.debug('...no new genes, preventing update')
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
    filter_dropdowns = []
    for col in filter_cols:
        idd = f"{id_prefix}-{table_str}-filter-{col}"
        LOG.debug(f"Filter dropdown: register ID={idd}")

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