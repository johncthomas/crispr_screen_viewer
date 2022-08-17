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
    DataSet,
    LOG,
    style_gene_selector_div,
    get_cmdline_options,
)
from crispr_screen_viewer.shared_components import (
    get_lab_val,
    get_annotation_dicts,
    register_gene_selection_processor,
    spawn_gene_dropdown,
)

from crispr_screen_viewer.selector_tables import (
    spawn_selector_tables,
    spawn_selector_tabs,
    spawn_treatment_reselector,
    get_selector_table_filter_keys,
    register_exptable_filters_comps,

)
from shared_components import spawn_filter_dropdowns

import dash_bootstrap_components as dbc

def initiate(app, data_set, public):
    """Source directory should contain the relevant info: metadata.csv,
    screen_analyses and expyaml directories."""
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

    def get_xy_choice_panel(xy:str) -> List[Div]:
        """Return button and output Divs.
        Button: id="cm-choose-{xy.lower()}-text"
        Para:   id="cm-chosen-{xy.lower()}-text"
        """
        xy = xy.lower()
        XY = xy.upper()
        xybutton = html.Button(f'Choose {XY} treatment', id=f'{PAGE_ID}-choose-{xy}', n_clicks=0)
        return [
            Div([xybutton]),
            Div(html.P(
                id=f"{PAGE_ID}-chosen-{xy}-text",
                children=f'No {XY} treatment selected. Choose from list and press "Choose" button'
            ))
        ]

    ## IF we add another analysis method that's appropriate for this
    #    then this can be uncommented
    # analysis_selector = [
    #     html.Label('Analysis type', htmlFor='analysis-selector'),
    #     dcc.RadioItems(
    #         id='analysis-selector',
    #         options=[{'value':a, 'label':data_set.analysis_labels[a]} for a in data_set.available_analyses],
    #         value='drz', # set all selected by default
    #         labelStyle={'display':'inline-block'}
    #     ),
    #     html.P([''], id='missing-analysis'),
    # ]

    scatter_chart = Div([dcc.Graph(
        id='graph',
        config={
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d',
                                       'autoScale2d', 'resetScale2d'],
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
                    data={'xk':None, 'yk':None})
    previous_row_store = dcc.Store('previous-rows',
                              data=[])


    exptab, comptab = spawn_selector_tabs(
        app,
        PAGE_ID,
        filter_dropdowns,
        selctr_tables,
        exp_tab_text=('Choose an experiment from the table below. This will '
                       'filter the available treatments in the "Select Treatments" tab. '
                       'Optionally, go straight to "Select Treatments" to see all possible '
                       'samples.'),
        comp_tab_text=('Choose treatments by selecting a row and then pressing the '
                   '"Choose X/Y Treatment" buttons. Once an X and Y treatment '
                   'has been selected, the comparison will appear in the "Chart" '
                   'and "Table" tabs.'),
        comp_choice_panel=[Div(get_xy_choice_panel('x')),
                           Div(get_xy_choice_panel('y'))],
    )

    ################
    # LAYOUT
    tabs = dcc.Tabs(

        id=f'{PAGE_ID}-tabs',
        value=f'{PAGE_ID}-exp-tab',
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
        # Div(
        #     style={'display':'inline-block'},
        #     children=[
        #         Div(selectors, style={'display':'inline-block'}),
        # ]),
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

    # When 'cm-choose-x/y' pressed (n_clicks), updated ('x/y-selector', 'value'&'options')
    #   with the comp chosen via ('cm-comp-table', 'selected_rows')

    # for updating what the chosen treatments are
    @app.callback(
        # dropdown selections/options and text below button
        Output(f"{PAGE_ID}-x-selector", 'value'),
        Output(f"{PAGE_ID}-x-selector", 'options'),
        Output(f'{PAGE_ID}-chosen-x-text', 'children'),

        Output(f"{PAGE_ID}-y-selector", 'value'),
        Output(f"{PAGE_ID}-y-selector", 'options'),
        Output(f'{PAGE_ID}-chosen-y-text', 'children'),

        Input(f'{PAGE_ID}-choose-x', 'n_clicks'),
        Input(f'{PAGE_ID}-choose-y', 'n_clicks'),

        State(f'{PAGE_ID}-comp-table', 'selected_rows'),
        State(f'{PAGE_ID}-comp-table', 'data'),
        State(f'{PAGE_ID}-x-selector', 'options'),
        #State(f'{PAGE_ID}-y-selector', 'options'),

    )
    def select_treat_for_cm(
            xbutton, ybutton, selected_row, table_data,
            opt
    ):
        if not selected_row:
            PreventUpdate

        try:
            triggered = callback_context.triggered_id
        except AttributeError:
            # v<2.4
            triggered = callback_context.triggered[0]['prop_id'].split('.')[0]
        if not triggered:
            raise PreventUpdate

        xy = triggered.split('-')[-1]
        compid = table_data[selected_row[0]]['Comparison ID']
        LOG.debug(f'Choosing {compid} for axis {xy}')

        if compid not in [d['value'] for d in opt]:
            opt.append(
                {'value':compid, 'label':compid}
            )

        if xy == 'x':
            return [compid, opt, compid,
                    dash.no_update, opt, dash.no_update]
        else:
            return [dash.no_update, opt , dash.no_update,
                    compid, opt, compid]

    def get_xyscores_genes(xk, yk, selected_analysis_type='drz'):
        """return score series from Dataset with unified indexes"""
        score_fdr = data_set.get_score_fdr(selected_analysis_type)

        x, y = [score_fdr['score'][k].dropna() for k in (xk, yk)]
        shared_genes = x.index.intersection(y.index)
        x, y = [xy[shared_genes] for xy in (x,y)]
        return x, y


    # update the keys stored in the data store, and gene selection options
    @app.callback(
        Output('selected-comps', 'data'),
        Output(f'{PAGE_ID}-gene-dropdown', 'options'),
        Input(f'{PAGE_ID}-x-selector', 'value'),
        Input(f'{PAGE_ID}-y-selector', 'value'),
    )
    def update_selection(xk, yk):
        LOG.debug(f'update_selection({xk}, {yk})')
        if (not xk) or (not yk):
            LOG.debug('not updating')
            raise PreventUpdate
        x, _ = get_xyscores_genes(xk, yk, 'drz')
        return {'xk':xk, 'yk':yk}, get_lab_val(x.index)

    # enable selecting genes by interacting with the graph
    register_gene_selection_processor(
        app, 'graph',
        State(f'{PAGE_ID}-gene-dropdown', 'value'),
        Output(f'{PAGE_ID}-gene-dropdown', 'value'),
    )


    # @app.callback(
    #     Output('comparison-selectors-div', 'style'),
    #     Output(f'{PAGE_ID}-gene-selector-div', 'style'),
    #     Input(f'{PAGE_ID}-tabs', 'value'),
    # )
    # def hide_comp_gene_selectors(selected_tab):
    #     """Set comparison/gene selector dropdown container's
    #     style={display: none}
    #     when selection tabs are active, show them when chart or table tabs
    #     are active."""
    #     if not selected_tab:
    #         raise PreventUpdate
    #     # any of these should be fine, if we add things later
    #     LOG.debug(f"Selected tab = {selected_tab}")
    #     show_treat_select = False
    #     for n in ('chart', 'graph', 'table'):
    #         if n in selected_tab:
    #             show_treat_select = True
    #     if show_treat_select:
    #         return (style_comparisons_card, style_gene_selector_div)
    #     else:
    #         return (style_hidden, style_hidden)

    # update the chart and table
    @app.callback(
        Output('graph', 'figure'),
        Output('table-div', 'children'),
        #Output('missing-analysis', 'children'),
        Input('selected-comps', 'data'),
        Input(f'{PAGE_ID}-gene-dropdown', 'value'),

    )
    def update_chart_table(selected_comps, selected_genes):
        xk, yk = selected_comps['xk'], selected_comps['yk']

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
        metadeets_str = "Experiment ID {}, {} library"
        xymet = [metadeets_str.format(row['Experiment ID'], row['Library']) for row in (rowx, rowy)]

        if xymet[0] == xymet[1]:
            xy_metadeets = xymet[0]
        else:
            xy_metadeets = f"X: {xymet[0]}<br>Y: {xymet[1]}"
        title = (f"{treatlines[0]} Vs {treatlines[1]}<br>"+xy_metadeets
                )
        LOG.debug(title)

        # axis labels
        xy_labs = [f"{comparisons.loc[k, 'Ctrl samp']} ➤ {comparisons.loc[k, 'Treat samp']}"
                   for k in (xk, yk)]

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
                # style_cell={ # doesn't work...
                #     # all three widths are needed
                #     'minWidth': '150px', 'width': '150px', 'maxWidth':'180px',
                #     'whiteSpace': 'normal'
                # },
                selected_rows=[],
                row_selectable='multi',
                export_format='csv',
            )
            table_output = [delta_table]
        else:
            table_output = dash.no_update

        return fig, table_output

    return layout




def launch_page(source, port, debug):

    # pycharm debug doesn't like __name__ here for some reason.
    app = dash.Dash('comparison_maker', external_stylesheets=[dbc.themes.BOOTSTRAP])

    if debug:
        LOG.setLevel(logging.DEBUG)

    source_directory = pathlib.Path(source)

    data_set = DataSet(source_directory)

    app.layout = initiate(app, data_set, public=True)

    app.run_server(debug=debug, host='0.0.0.0', port=int(port), )

if __name__ == '__main__':

    launch_page(*get_cmdline_options())



    # # create strings that will be used in the experiment selector dropdown,
    # #   gives some info about the experiment
    # experiment_options = []
    # for exp_id in sorted(comparisons['Experiment ID'].unique()):
    #     comps_for_exp = comparisons.loc[comparisons['Experiment ID'] == exp_id]
    #
    #     expd = experiments_metadata.loc[exp_id]
    #
    #     if 'Treatments' in expd:
    #         treats = expd['Treatments']
    #     else:
    #         # get some of the treatments at least...
    #         treats = sorted(comps_for_exp.Treatment.dropna().apply(lambda t: t.split(' ➤ ')[-1]).unique())
    #         treats = ", ".join([t for t in treats if t and t not in ('DMSO', 'None')])
    #
    #     maxlen = 60
    #     if len(treats) > maxlen:
    #         treats = treats[:maxlen-3]+'...'
    #
    #     # get a nice string
    #     #exp_metadata = data_set.experiments_metadata.loc[exp_id]
    #     exp_metadata = experiments_metadata.loc[exp_id]
    #     exp_metadata = exp_metadata.dropna()
    #     #who = exp_metadata['Investigator']
    #     date_k = exp_metadata.index.map(lambda s: s.startswith('Date'))
    #     date_k = exp_metadata.index[date_k]
    #     when = exp_metadata[date_k].values[0]
    #     exp_str = f"{exp_id}, {when} – {treats}"
    #
    #     experiment_options.append(
    #         {'label':exp_str, 'value':exp_id}
    #     )

    ##
    ## DONT DELETE
    ## these are from when comparisons were selected by dropdowns. Might bring em back to
    ## quickly flick between treatments
    # when experiment is selected, update the group options
    # @app.callback(
    #     Output('timepoint-selector', 'options'),
    #     Output('timepoint-selector', 'disabled'),
    #     Input('experiment-selector', 'value')
    # )
    # def update_timepoint_options(exp_id):
    #     if not exp_id:
    #         return [], True
    #     timepoints = comparisons.loc[comparisons['Experiment ID'] == exp_id, 'Timepoint'].unique()
    #     return [{'value':tp, 'label':tp} for tp in timepoints], False
    #
    # # when an exp and timepoint is selected, update the X/Y options
    # @app.callback(
    #     Output('xy-selectors', 'children'),
    #     Input('experiment-selector', 'value'),
    #     Input('timepoint-selector', 'value'),
    # )
    # def update_xy_selectors(exp_id, timepoint):
    #     if (not exp_id) or (not timepoint):
    #         raise PreventUpdate
    #
    #     comps = comparisons.loc[
    #         (comparisons.Timepoint == timepoint)
    #         & (comparisons['Experiment ID'] == exp_id)
    #     ]
    #
    #     comp_options = []
    #     for idx, row in comps.iterrows():
    #         # format the KO part of the cell description, emptry string if none
    #         if '-KO' not in row['Treatment']:
    #             if row['KO'] == 'WT':
    #                 ko = ''
    #             else:
    #                 ko = f" {row['KO']}-KO"
    #         else:
    #             ko = ''
    #         # nice strings
    #         treatstr = f"{row['Ctrl samp']} ➤ {row['Treat samp']}"
    #         if row['Timepoint'] == 'endpoints':
    #             comp_description = f"{treatstr} – Effect of {row['Treatment']} in {row['Cell line']}{ko} cells"
    #         else:
    #             comp_description = f"{treatstr} – {row['Treatment']} in {row['Cell line']}{ko} cells"
    #         comp_options.append(
    #             {'value':idx, 'label':comp_description}
    #         )
    #
    #     return get_xy_selectors(comp_options)

# group_dropdown = dcc.Dropdown(
#     id='timepoint-selector',
#     placeholder="Select timepoint group",
#     disabled=True,
#     style={'width':'300px', 'display':'inline-block'},
#     value=[],
#     options=[]
# )

    # # # get selections from the table
    ## # DONT DELETE This might come back later - for now prioritising selection via the graph.
    # @app.callback(
    #     Output(f'{PAGE_ID}-gene-dropdown', 'value'),
    #     Output('previous-rows', 'data'),
    #     Input('table', 'selected_rows'),
    #     State('previous-rows', 'data'),
    #     State(f'{PAGE_ID}-gene-dropdown', 'value'),
    #     State('table', 'data'),
    # )
    # def select_gene_with_table(selected_rows, prev_rows, selected_genes, table_data):
    #     """Whenever a row is selected or deselected from the table
    #     add or remove that gene from the selected genes. Additions
    #     made to or from the dropdown do not effect the table"""
    #     if selected_genes is None:
    #         selected_genes = []
    #     print(selected_genes)
    #     LOG.debug('select_gene_with_table()')
    #     selected_rows = set(selected_rows)
    #     row_selected_genes = set([table_data[i]['Gene'] for i in selected_rows])
    #     prev_rows = set(prev_rows)
    #     new_gene = row_selected_genes.difference(prev_rows)
    #     # Writing this like it's possible to do multiple operations at once
    #     output_selection = selected_genes+list(new_gene)
    #     removed_genes = prev_rows.difference(row_selected_genes)
    #     output_selection = [g for g in output_selection if g not in removed_genes]
    #     if not new_gene and not removed_genes:
    #         LOG.debug('No update to selected genes')
    #         raise PreventUpdate
    #     LOG.debug('select_gene_by_table -> '+f"{output_selection}, {list(row_selected_genes)}")
    #     return output_selection, list(row_selected_genes)

    # exp_dropdown = dcc.Dropdown(
    #     id='experiment-selector',
    #     placeholder="Select experiment",
    #     style={'width':'700px', 'display':'inline-block'},
    #     value=[],
    #     options=experiment_options
    # )
# @app.callback(
#     Output('xy-selectors', 'children'),
#     Input('experiment-selector', 'value'),
#     Input('timepoint-selector', 'value'),
# )
# def update_xy_selectors(exp_id, timepoint):
#     if (not exp_id) or (not timepoint):
#         raise PreventUpdate
#
#     comps = comparisons.loc[
#         (comparisons.Timepoint == timepoint)
#         & (comparisons['Experiment ID'] == exp_id)
#     ]
#
#     comp_options = []
#     for idx, row in comps.iterrows():
#         # format the KO part of the cell description, emptry string if none
#         if '-KO' not in row['Treatment']:
#             if row['KO'] == 'WT':
#                 ko = ''
#             else:
#                 ko = f" {row['KO']}-KO"
#         else:
#             ko = ''
#         # nice strings
#         treatstr = f"{row['Ctrl samp']} ➤ {row['Treat samp']}"
#         if row['Timepoint'] == 'endpoints':
#             comp_description = f"{treatstr} – Effect of {row['Treatment']} in {row['Cell line']}{ko} cells"
#         else:
#             comp_description = f"{treatstr} – {row['Treatment']} in {row['Cell line']}{ko} cells"
#         comp_options.append(
#             {'value':idx, 'label':comp_description}
#         )
#
#     return get_xy_selectors(comp_options)