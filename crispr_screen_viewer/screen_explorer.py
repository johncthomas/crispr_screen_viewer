#!/usr/bin/env python

import logging
import copy
import pandas as pd
import numpy as np

from dash import dash, dcc, html, Input, Output, State, dash_table, callback_context

from dash.exceptions import PreventUpdate

import dash_bootstrap_components as dbc

import plotly.graph_objs as go

import pathlib, os
from dash.dependencies import Input, Output, State
from typing import Collection, Union, Dict, List, Callable

from crispr_screen_viewer.functions_etc import (
    datatable_column_dict,
    get_cmdline_options,
    get_table_title_text,
    get_treatment_label
)

from crispr_screen_viewer.shared_components import (
    get_gene_dropdown_lab_val,
    get_annotation_dicts,
    LOG,
    get_stat_source_selector,
    register_gene_selection_processor,
    spawn_gene_dropdown,
)

from crispr_screen_viewer.selector_tables import (
    spawn_selector_tables,
    spawn_selector_tabs,
    spawn_treatment_reselector,
    get_selector_table_filter_keys,
)
from crispr_screen_viewer.shared_components import spawn_filter_dropdowns
from crispr_screen_viewer.dataset import DataSet
Div = html.Div


idprfx_res_table = 'gene-results-table'
# String used to distinguish reused components across pages. Component IDs used by dash
#   preceeded by this string
PAGE_ID = 'se'

def spawn_volcano_graph(app, fig_id=f'{PAGE_ID}-volcano'):
    """Return layout containing the plotly Graph object and stat/gene
    selectors. Register volcano chart callbacks.

    Callbacks:
        render_volcano: produces volcano figure data, to plot selected data.
            Gets data (and is triggered by changes to graph-data.data)
        add_volcano_selection: deals with gene selections coming from users
            interacting with the graph
        calls register_gene_selection_processor(app, fig_id)
    """
    volcano = dcc.Graph(
        id=fig_id,
        config={
            'editable':True,
            'edits':{'annotationPosition':False},
        },
        style={'width': '1200px',
               'height':'800px',
               'padding':'0px'},
        # # note, figure gets replaced in update_volcano callback so add
        # # any changes there.
        figure={'layout':{'clickmode':'event+select',
                          'dragmode':'select'}, }
    )

    volcano_layout = [
        Div([volcano]),
    ]

    ####
    # render volc should trigger only when a comparison has been selected
    # another call back should deal with selection, look at comparison maker
    # Some more of that should be shared, like the title generator.
    @app.callback(
        Output(f'{PAGE_ID}-volcano', 'figure'),

        Input(f'{PAGE_ID}-gene-dropdown', 'value'),
        Input(f'{PAGE_ID}-graph-data', 'data'),

        State(f'{PAGE_ID}-comp-table', 'data'),
        State(f'{PAGE_ID}-comp-table', 'selected_rows'),
    )
    def render_volcano(selected_genes, xy_genes,
                       table_data, selected_row):

        LOG.debug(f"Rendering volcano")

        if not selected_row:
            raise PreventUpdate

        # get values and create the figure
        x, fdr, genes = [xy_genes[k] for k in ('score', 'fdr', 'genes')]
        # dcc.Store converts this to a list...
        x, fdr = [pd.Series(xy, index=genes) for xy in (x,fdr)]
        y = fdr.apply(lambda _x: -np.log10(_x))

        fig = go.Figure(
            data=go.Scattergl(
                x=x.values,
                y=y.values,
                mode='markers',
                customdata=fdr,
                text=genes,
                hovertemplate= (
                        "<b>%{text}</b><br>" +
                        "LFC: %{x:.2f}<br>" +
                        "FDR: %{customdata:.2e}" +
                        "<extra></extra>"
                    ),
                ),
            layout={'clickmode':'event+select',
                    'dragmode':'select'},

        )
        LOG.debug(f"{(len(x), len(y))}")

        # add some titles
        # consider using data_set.comparisons.loc[<Input('selected-comp', 'data')>]
        row = table_data[selected_row[0]]
        line1,line2 = get_treatment_label(row)
        title = f"<b>{line1}</b><br>{line2}"
        fig.update_layout(
            title=title,
            xaxis_title='LFC',
            yaxis_title='-Log10(FDR)',

        )

        # Add annotations for the selected genes
        new_annotations = get_annotation_dicts(x[selected_genes], y[selected_genes], selected_genes)
        for anot in new_annotations:
            fig.add_annotation(
                **anot
            )

        LOG.debug(str(fig))
        LOG.debug('Finished generating volcano Figure')
        return fig


    return volcano_layout



def initiate(app, data_set:DataSet, public=False) -> Div:
    """Source directory should contain the relevant info: metadata.csv,
    screen_analyses and expyaml directories."""

    comparisons = data_set.comparisons
    # try:
    #     comparisons = comparisons.drop('Available analyses', 1)
    # except:
    #     pass


    # **GRAPHS****
    volcano_layout = spawn_volcano_graph(app, f'{PAGE_ID}-volcano')

    # # ***TABLES***
    # # Timepoint, renaming values to something more readable
    # for val, lab in timepoint_labels.items():
    #     m = comparisons.Timepoint.str.startswith(val)
    #     comparisons.loc[m, 'Time point group'] = lab

    # The Experiments and Comparisons tables keys 'exp' and 'comp'
    # User selection of row in exp table switches to comp tab, selection in comp tab
    #   updates plot and datatable.

    filter_keys = get_selector_table_filter_keys(public)
    selctr_tables = spawn_selector_tables(
        app, data_set, filter_keys, public, PAGE_ID)

    datatable = dash_table.DataTable(
        id=idprfx_res_table+'-table',
        columns=[datatable_column_dict(x) for x in ('Effect size', 'FDR', 'Selected')],
        sort_action='native',
        sort_mode='multi',
        filter_action='native',
        # not sure what this does but I put it in the other table for some reason
        css=[{"selector": "p", "rule": "margin: 0"}],
    )

    results_table_layout = [
        Div([], id=idprfx_res_table+'-title',),
        Div(datatable, style={'padding-right':'20px'})
    ]

    filter_dropdowns = {tabk:spawn_filter_dropdowns(PAGE_ID, tabk, filtercols, comparisons)
                        for tabk, filtercols in filter_keys.items()}

    exptab, comptab = spawn_selector_tabs(
        app, PAGE_ID, filter_dropdowns, selctr_tables,
        ('Choose an experiment below to filter the options in "Select Treatment" table. '
         'Go straight to Select Treatment to see all options.'),
        ('Select a specific treatment using the table below. Click on tabs to '
         'the right to see results for a selected treatment'),
    )



    # ############## #
    # ****LAYOUT**** #
    # ############## #


    tabs = dcc.Tabs(
        id=f'{PAGE_ID}-tabs', value=f'{PAGE_ID}-comp-tab', children=[
        exptab,
        comptab,

        # volcano
        dcc.Tab(
            label='Volcano plot', value='volcano-chart-tab',
            className='data-tab', selected_className='data-tab--selected',
            children=volcano_layout
        ),
        dcc.Tab(
            label='Results Table', value=idprfx_res_table+'-tab',
            className='data-tab', selected_className='data-tab--selected',
            children=results_table_layout
        ),
    ])


    se_layout = Div([
        Div([
            html.H1("Explore screen results"),
            html.P("View results from a single treatment. Select data using the blue tabs, and view results in the green tabs."),
            Div([
                # inline block so the results control panel goes to the side of it
                Div([tabs,], style={'display':'inline-block',}),
                Div([
                    spawn_treatment_reselector(PAGE_ID, is_xy=False),
                    get_stat_source_selector(PAGE_ID, 'Analysis method:'),
                ], id=f'{PAGE_ID}-results-control-panel')
            ],  className='tab-content-box'),

            spawn_gene_dropdown(app, PAGE_ID),



            dcc.Store(id=f'{PAGE_ID}-selected-exp', storage_type='session'),
            dcc.Store(id=f'{PAGE_ID}-selected-comp', storage_type='session'),
            dcc.Store(id=f'{PAGE_ID}-graph-data', storage_type='session'),
            dcc.Store(id=f'{PAGE_ID}-selected-genes', storage_type='session'),
        ]),
    ], style={'display':'inline-block'})


    # ***CALLBACKS***

    # Note on gene selection dropdown callbacks:
    #   I want them all to be updated together (currently 2, but am building for an arbitrary number)
    #   but I only want the currently selected graph to be updated. So, the call down for the current
    #   tab takes Input from the current tab's gene dropdown and Outputs to the next tabs dropdown.
    #   The final one outputs to a data storage object, which triggers a callback checking if all
    #   of the gene dropdowns have the same values. If they do, raise PreventUpdate, if not, pass
    #   selections to the first tab's gene dropdown.

    register_gene_selection_processor(
        app, f"{PAGE_ID}-volcano",
        State(f'{PAGE_ID}-gene-dropdown', 'value'),
        Output(f'{PAGE_ID}-gene-dropdown', 'value')
    )

    @app.callback(
        Output(f'{PAGE_ID}-graph-data', 'data'),
        Output(f'{PAGE_ID}-gene-dropdown', 'options'),

        Input(f'{PAGE_ID}-comp-selector', 'value' ),
        Input(f'{PAGE_ID}-stat-source-selector', 'value'),
    )
    def store_selected_data(compid, sig_source, ):
        """Output data for the charts.

        Returns:
        {'x':score, 'y':fdr, 'genes':"""

        args_for_printing = {k:v for k, v in zip(
            'selected_row, sig_source, table_stat_source'.split(', '),
            [compid, sig_source, ]
        )}
        LOG.debug(f'CALLBACK: update_volcano_data with {args_for_printing}')

        if not compid:
            raise PreventUpdate

        # get x, y and genes values
        # todo LFC as it's own table!
        score_fdr = data_set.get_score_fdr('mag', sig_source)
        score, fdr = [score_fdr[k][compid].dropna() for k in ('score', 'fdr')]

        # some genes may get filtered out
        unified_index = score.index.intersection(fdr.index)
        score,fdr = [xy.reindex(unified_index) for xy in (score,fdr)]

        volcano_data = {'score': score, 'fdr': fdr, 'genes': score.index}
        gene_options = get_gene_dropdown_lab_val(data_set, score.index)

        LOG.debug(f'End of update_volcano_data with:')
        LOG.debug('     datatable:  '+'\n'.join([f"{k}={volcano_data[k].head()}" for k in ('score', 'fdr')]))

        return (
            volcano_data,
            gene_options,
        )

    # update the contents of the gene results tab
    @app.callback(
        Output(idprfx_res_table+'-table', 'columns'),
        Output(idprfx_res_table+'-table', 'data'),
        Output(idprfx_res_table+'-title', 'children'),

        Input(f'{PAGE_ID}-comp-selector', 'value' ),
        Input(f'{PAGE_ID}-gene-dropdown', 'value'),
        Input(f'{PAGE_ID}-stat-source-selector', 'value'),
        Input(f'{PAGE_ID}-tabs', 'value'),

    )
    def update_results_table_data(selected_comp, selected_genes, stat_source,
                                  selected_tab, ):
        LOG.debug(f'CALLBACK: update_results_table_data({selected_comp}, {stat_source})')
        if not selected_comp:
            raise PreventUpdate
        if selected_tab != idprfx_res_table+'-tab':
            raise PreventUpdate

        # labels for the analysis type
        ans_lab = data_set.analysis_labels[stat_source]
        score_lab = data_set.score_labels[stat_source]

        # data for the table
        dat = data_set.get_score_fdr(stat_source, stat_source)
        lfc = data_set.get_score_fdr('mag', 'mag')['score'][selected_comp]
        score = dat['score'][selected_comp]
        fdr = dat['fdr'][selected_comp]

        # index is genes
        index = lfc.index.union(score.index).union(fdr.index)
        is_selected = index.map(lambda g: ['❌','Selected'][g in selected_genes])

        # Build the table
        results_tab = pd.DataFrame(
            {
                'Log2(FC)':lfc,
                score_lab:score,
                'FDR':fdr,
                'Gene selected':is_selected,
            },
            index=index
        )
        no_stats = results_tab[[score_lab, 'FDR']].isna().any(axis=1)
        results_tab = results_tab.loc[~no_stats]
        results_tab.index.name = 'Gene'

        # Create the Output objects
        results_data = results_tab.reset_index().to_dict('records')

        columns = [datatable_column_dict(x) for x in results_data[0].keys()]


        treatment_para = get_table_title_text(
            comparisons.loc[selected_comp],
            ans_lab
        )

        return (columns, results_data, treatment_para, )

    @app.callback(
        Output(f'{PAGE_ID}-comp-selector', 'value' ),
        Output(f'{PAGE_ID}-comp-selector', 'options'),
        Input(f'{PAGE_ID}-comp-table', 'selected_rows'),
        State(f'{PAGE_ID}-comp-table', 'data'),
        State(f'{PAGE_ID}-comp-selector', 'options'),
    )
    def select_comp(selected_row, table_data, opt):

        if not selected_row:
            raise PreventUpdate

        # get the comparison ID, select the relevant data from dataset
        compid = table_data[selected_row[0]]['Comparison ID']
        LOG.debug(f'selected comp: {compid}')
        if compid not in [d['value'] for d in opt]:
            opt.append(
                {'value':compid, 'label':compid}
            )

        return compid, opt

    return se_layout




if __name__ == '__main__':
    from crispr_screen_viewer.functions_etc import launch_page
    launch_page(*get_cmdline_options(), 'ScreensExplorer', initiate)

