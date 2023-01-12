#!/usr/bin/env python


import flask
from typing import List, Dict, Tuple, Collection
import typing

import pandas as pd
import numpy as np

from dash.exceptions import PreventUpdate
from dash import dash, dcc, html, dash_table, callback_context, callback
Div = html.Div
import plotly.graph_objs as go

import pathlib, os
from dash.dependencies import Input, Output, State


from crispr_screen_viewer.functions_etc import (
    LOG,
    get_cmdline_options,
    html_small_span,
)
from crispr_screen_viewer.shared_components import (
    get_lab_val,
    get_annotation_dicts,
    register_gene_selection_processor,
    spawn_gene_dropdown,

)


import dash_bootstrap_components as dbc

# will try to make generic and reusable as possible, but options are always dependent on the data
# Data is a dictionary of DF.
# * Each key a comp
# * Each table containing the results of the comparison analysis. LFC, adj-p, gene symbol.
# because, later on, we'll want to have duplicate rows in the pure DE proteomics, we will not use meaningful indicies.

# logger.add(sys.stderr, level='DEBUG',
#                format='<level>{message}</level> | '
#                       '<level>{level: <8}</level> | '
#                       '<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan>  '
#                       '<green>({time:YYYY-MM-DD HH:mm:ss})</green>'
#            )

def unify_indicies(tabx, taby):
    """return tables with aligned rows and reset Index."""

    # Indexes are reset so that I don't accidently
    #  use them for something important

    shared = tabx.index.intersection(taby.index)
    return [tab.loc[shared].reset_index(drop=False) for tab in (tabx, taby)]


def get_customdata(tabs):
    shared = tabs[0].loc[:, ['Gene', 'Accessions', 'Description', 'GeneSymbol']]
    chunks = [shared]
    # the stat columns prefixed with 'X_' or 'Y_'
    for xy, tab in ('X', tabs[0]), ('Y', tabs[1]):
        tab = tab.loc[:, 'log2FC	AvgIntensity	adj.P.Val	B	log10FDR'.split('\t')].copy()
        tab.columns = [f"{xy}_{col}" for col in tab.columns]
        chunks.append(tab)

    return pd.concat(chunks, axis=1)


def empty_figure(width=1200, height=1200):

    return go.Figure(
        layout=dict(
            width=width,
            height=height,
            clickmode='event+select',
            dragmode='select',
        )
    )


def load_data(dir):
    DATA = {}
    prevdir = os.getcwd()
    os.chdir(dir)
    for fn in os.listdir():
        if fn.endswith('.csv'):
            DATA[fn[:-4]] = pd.read_csv(fn, index_col=0)
    os.chdir(prevdir)
    return DATA


def spawn_comp_selector(comps) -> List[Div]:
    """list of Divs of the Dropdowns with Labels.

    IDs (x|y)-selector, used as Input for selecting data"""

    output_selectors = []

    for xy in 'xy':
        label = {'x':'X treatment',
                 'y':'Y treatment',}[xy]
        id=f'{xy}-selector'
        xy_selector = dcc.Dropdown(
            id=id,
            placeholder=f"Select {label}",
            className='comp-selector',
            options=get_lab_val(comps),
            #value={'x':'MAEA_KO_vs_WT', 'y':'MAEA_KO_vs_MAEA_M396R'}[xy] #for debug
            value=None,
        )
        label = html.Label(f'{label}:',htmlFor=id,style= {'display':'block'})
        block = Div([label, xy_selector,])
        output_selectors.append(block)

    return output_selectors


def add_annotations(fig, x, y, labels, selected):
    """call fig.add_annotations, and return the row mask.

    x, y, labels should share index

    rows annotated are labels.isin(selected)
    """
    m = labels.isin(selected)
    new_annotations = get_annotation_dicts(
        x[m], y[m], labels[m])
    for anot in new_annotations:
        fig.add_annotation(
            **anot
        )
    return m

def spawn_scatter(
        results, labelcol='GeneSymbol',) -> dcc.Graph:

    @callback(
        Output('scatter-graph', 'figure'),
        Output('prot-gene-dropdown', 'options'),
        Output('selected-table', 'data'),
        Input('x-selector', 'value'),
        Input('y-selector', 'value'),
        Input('prot-gene-dropdown', 'value'),
        Input('stat-selector', 'value'),
    )
    def update_graph_table(xk, yk, selected_genes, selected_stat):

        LOG.debug(f"{xk} {yk}")
        LOG.debug(f"selected_genes: {selected_genes}")

        if (xk is None) or (yk is None):
            raise PreventUpdate

        tabs = unify_indicies(results[xk], results[yk])

        # *FIGURE*
        fig = empty_figure()
        customdat = get_customdata(tabs)

        def cdstr(col, fmt=':.3'):
            i = customdat.columns.get_loc(col)
            return f"%{{customdata[{i}]{fmt}}}"

        x = tabs[0][selected_stat]
        y = tabs[1][selected_stat]

        fig.add_trace(
            go.Scattergl(
                x=x, y=y,
                text=tabs[0][labelcol],
                mode='markers',
                customdata=customdat,
                hovertemplate=(
                    f"<b>{cdstr('GeneSymbol', '')}</b><br>"
                    f"X LFC: {cdstr('X_log2FC')}<br>"
                    f"X FDR: {cdstr('X_adj.P.Val')}<br>"
                    "<br>"
                    f"Y LFC: {cdstr('Y_log2FC')}<br>"
                    f"Y FDR: {cdstr('Y_adj.P.Val')}<br>"
                    "<extra></extra>"
                )
            )
        )
        #todo data labels specificity not good
        axistitle_fmt = lambda k: f"{k.replace('_vs_', ' vs ')} ({selected_stat})"
        fig.update_layout(
            dict(xaxis_title=axistitle_fmt(xk),
                 yaxis_title=axistitle_fmt(yk))
        )

        # *GENE DROPDOWN OPTIONS*
        # options for the gene dropdown. Rows will be selected by mask, not index
        opts = [{'label':v, 'value':v} for v in sorted(tabs[0][labelcol].dropna().unique())]

        LOG.debug(opts[:5])

        # *ANNOTATIONS*
        # rows for annotation selected by having `labelcol` row value in selected_genes
        m = add_annotations(fig, x, y, tabs[0][labelcol], selected_genes)

        # *DATATABLE DATA*
        # construct a new df and export records
        if m.sum():
            newdat = [tabs[0].loc[m, ['Gene', 'GeneSymbol']]]
            for tab in tabs:
                newdat.append(tab.loc[m, ['log2FC', 'adj.P.Val']])
            newtab = pd.concat(newdat, axis=1, )
            newtab.columns = ['Protein', 'Gene', 'X LFC', 'X FDR', 'Y LFC', 'Y FDR']
            newdat = newtab.to_dict('records')
        else:
            newdat = dash.no_update

        return fig, opts, newdat

    graph = dcc.Graph(
        id='scatter-graph',
        config={
            'editable':True,
            'edits':{'annotationPosition':False},
        },
        figure=empty_figure()
    )

    return graph


def spawn_volcano(results, xy, dim=(600, 800), labelcol='GeneSymbol'):
    graph = dcc.Graph(
        id=f'{xy}-volcano',
        config={
            'editable':True,
            'edits':{'annotationPosition':False},
        },
        figure=empty_figure(*dim)
    )

    @callback(
        Output(f'{xy}-volcano', 'figure'),
        Input(f'{xy}-selector', 'value'),
        Input('prot-gene-dropdown', 'value'),
    )
    def update_volc(comp, selected_genes):
        tab = results[comp]
        lfc = tab['log2FC']
        fdr = tab['adj.P.Val'].apply(lambda p: -np.log10(p))
        fig = empty_figure(*dim)

        #todo hovertext shows protein.
        fig.add_trace(
            go.Scattergl(
                x=lfc, y=fdr,
                text=tab[labelcol],
                mode='markers',
            )
        )

        add_annotations(fig, lfc, fdr, tab[labelcol], selected_genes)

        fig.update_layout(
            dict(xaxis_title='LFC',
                 yaxis_title='-log10(FDR)',
                 title=comp.replace('_vs_', ' vs '))
        )

        return fig

    return graph


def initiate(app, results:Dict[str, pd.DataFrame], port):

    scatter_plot = spawn_scatter(results)
    comp_dropdowns = spawn_comp_selector(results.keys())
    gene_selector = spawn_gene_dropdown(app, 'prot')

    table = dash_table.DataTable(
        id='selected-table',
        columns=[
            dict(name=v, id=v) for v in ('Protein', 'Gene', 'X LFC', 'X FDR', 'Y LFC', 'Y FDR')
        ],
        export_format='csv',
        fill_width=False,
        sort_action='native',
    )

    stats={'LFC':'log2FC', 'FDR':'log10FDR'}
    stat_selector = dcc.RadioItems(
        id='stat-selector',
        options=[{'label': k, 'value': v} for k, v in stats.items()],
        style={'display':'inline-block'},
        value=stats['FDR'],
    )


    controls_row = dbc.Row([
        dbc.Card([
            dbc.CardBody(comp_dropdowns, ),
        ], style={'width': '400px'}),
        dbc.Card([
            dbc.CardHeader('Plotted statistic'),
            dbc.CardBody(stat_selector)
        ], style={'width': '100px'})
    ])

    volcanos = Div([
        spawn_volcano(results, xy,) for xy in 'xy'
    ], style={'display':'flex'})

    layout = Div([
        controls_row,
        Div(gene_selector),
        Div([
            Div([scatter_plot], style={'display':'inline-block'}),
            Div([table], style={'padding-top':'75px'}),
         ], style={'display':'flex'}),
        volcanos,
    ])

    # selections on the graph propogate to the dropdown
    register_gene_selection_processor(
        app,
        'scatter-graph',
        State('prot-gene-dropdown', 'value'),
        Output('prot-gene-dropdown', 'value'),
    )

    app.layout = layout
    return app


if __name__ == '__main__':
    import sys
    argv = sys.argv
    if (len(argv) < 2) or (argv[1] in ('-h', '--help')):
        print('proteomics_biplots.py DATA PORT\n\tDATA can be a path to a .pickle or dir with csv.')
        exit(1)

    if argv[1].endswith('.pickle'):
        import pickle

        with open(argv[1], 'rb') as f:
            data = pickle.load(f)
    elif os.path.isdir(argv[1]):
        data = load_data(argv[1])
    else:
        print(f"{argv[1]} is not a .pickle or directory.")
        exit(1)
    try:
        port = int(argv[2])
    except:
        print(f"{argv[2]} (port) needs to be a number.")
        exit(1)

    server = flask.Flask(__name__)

    app = dash.Dash(
        'proteomics_bioplots',
        server=server,
        external_stylesheets=[dbc.themes.BOOTSTRAP])
    app = initiate(app, data, port)
    app.run(host='0.0.0.0', port=port)
