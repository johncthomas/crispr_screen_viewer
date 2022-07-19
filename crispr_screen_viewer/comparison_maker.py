#!/usr/bin/env python
import logging
import sys

import pandas as pd
import numpy as np

from dash.exceptions import PreventUpdate
from dash import dash, dcc, html, dash_table, callback_context
Div = html.Div
import plotly.graph_objs as go

import pathlib, os
from dash.dependencies import Input, Output, State
import typing
from crispr_screen_viewer.functions_etc import DataSet
from crispr_screen_viewer.shared_components import (
    get_lab_val,
    get_annotation_dicts,
    LOG,
    timepoint_labels,
)


def initiate(app, data_set):
    """Source directory should contain the relevant info: metadata.csv,
    screen_analyses and expyaml directories."""


    comparisons = data_set.comparisons
    experiments_metadata = data_set.experiments_metadata
    print(experiments_metadata)

    ###############
    # COMPONENTS

    # create strings that will be used in the experiment selector dropdown,
    #   gives some info about the experiment
    experiment_options = []
    for exp_id in sorted(comparisons['Experiment ID'].unique()):
        comps_for_exp = comparisons.loc[comparisons['Experiment ID'] == exp_id]

        expd = experiments_metadata.loc[exp_id]

        if 'Treatments' in expd:
            treats = expd['Treatments']
        else:
            # get some of the treatments at least...
            treats = sorted(comps_for_exp.Treatment.dropna().apply(lambda t: t.split(' ➤ ')[-1]).unique())
            treats = ", ".join([t for t in treats if t and t != 'DMSO'])

        maxlen = 60
        if len(treats) > maxlen:
            treats = treats[:maxlen-3]+'...'

        # get a nice string
        #exp_metadata = data_set.experiments_metadata.loc[exp_id]
        exp_metadata = experiments_metadata.loc[exp_id]
        #who = exp_metadata['Investigator']
        date_k = exp_metadata.index.map(lambda s: s.startswith('Date'))
        date_k = exp_metadata.index[date_k]
        when = exp_metadata[date_k]
        exp_str = f"{exp_id}, {when} – {treats}"

        experiment_options.append(
            {'label':exp_str, 'value':exp_id}
        )

    exp_dropdown = dcc.Dropdown(
        id='experiment-selector',
        placeholder="Select experiment",
        style={'width':'700px', 'display':'inline-block'},
        value=[],
        options=experiment_options
    )

    group_dropdown = dcc.Dropdown(
        id='timepoint-selector',
        placeholder="Select timepoint group",
        disabled=True,
        style={'width':'300px', 'display':'inline-block'},
        value=[],
        options=[]
    )

    comp_selector_style = {'width':'800px', }
    def get_xy_selectors(options):
        if not options:
            disabled = True
        else:
            disabled = False

        x_selector = dcc.Dropdown(
            id='x-selector',
            placeholder="Select X treatment",
            disabled=disabled,
            style=comp_selector_style,
            value=[],
            options=options,
        )
        y_selector = dcc.Dropdown(
            id='y-selector',
            placeholder="Select Y treatment",
            disabled=disabled,
            style=comp_selector_style,
            value=[],
            options=options,
        )
        return [
            html.Label('3. X comparison:',htmlFor='x-selector'),
            x_selector,
            html.Label('4. Y comparison:',htmlFor='y-selector'),
            y_selector,
        ]

    selectors = [
        Div([
            Div([
                html.Label('1. Experimental dataset:',htmlFor='experiment-selector',),
                exp_dropdown,
            ], style={'display':'inline-block'}),
            Div([
                html.Label('2. Time point group:',htmlFor='timepoint-selector',),
                group_dropdown
            ], style={'display':'inline-block'}),
        ]),
        Div(get_xy_selectors([]), 'xy-selectors'),
    ]

    gene_dropdown = [
        Div([
            html.Label('5 (optional). Highlight genes :',htmlFor='gene-selector'),
            dcc.Dropdown(
                id='gene-selector',
                placeholder='Label genes',
                style={'width':'1000px', 'height':'100px', },
                value=[],
                options=[],
                clearable=True,
                multi=True,
            )
        ])
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
            'editable':True
        },
        style={'height': '1000px', 'width':'1000px', 'padding':'0px'},
        figure={'layout':{'clickmode':'event+select'}}
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
        row_selectable='single',
    )



    selected_store = dcc.Store('selected-comps',
                    data={'xk':None, 'yk':None})
    previous_rows = dcc.Store('previous-rows',
                              data=[])

    ################
    # LAYOUT
    tabs = dcc.Tabs(
        value='chart-tab',
        style={'width':1050, 'display':'inline-block'},
        children=[
            # comparison selector
            dcc.Tab(
                label='Chart',
                value='chart-tab',
                children=[
                    Div([scatter_chart])
                ]
            ),
            dcc.Tab(
                label='Table',
                value='table-tab',
                children=Div([delta_table], 'table-div', style={'width':'1000px'})
            )
        ]
    )

    layout = Div([
        html.H1('Compare treatments'),
        html.P('Choose two sets of results and plot the gene NormZ scores against each other. '
               'The two results must be from the same experiment and time point group so that the '
               'comparison is likely to be meaningful. This can be used to, for example, observe '
               'synthetic effects between two treatments.'),
        selected_store,
        previous_rows,
        Div(
            style={'display':'inline-block'},
            children=[
                Div(selectors, style={'display':'inline-block'}),
        ]),
        Div([tabs]),
        Div(gene_dropdown)
    ])


    ################
    # CALLBACKS
    # when experiment is selected, update the group options
    @app.callback(
        Output('timepoint-selector', 'options'),
        Output('timepoint-selector', 'disabled'),
        Input('experiment-selector', 'value')
    )
    def update_timepoint_options(exp_id):
        if not exp_id:
            return [], True
        timepoints = comparisons.loc[comparisons['Experiment ID'] == exp_id, 'Timepoint'].unique()
        return [{'value':tp, 'label':timepoint_labels[tp]} for tp in timepoints], False

    # when an exp and timepoint is selected, update the X/Y options
    @app.callback(
        Output('xy-selectors', 'children'),
        Input('experiment-selector', 'value'),
        Input('timepoint-selector', 'value'),
    )
    def update_xy_selectors(exp_id, timepoint):
        if (not exp_id) or (not timepoint):
            raise PreventUpdate

        comps = comparisons.loc[
            (comparisons.Timepoint == timepoint)
            & (comparisons['Experiment ID'] == exp_id)
        ]

        comp_options = []
        for idx, row in comps.iterrows():
            # format the KO part of the cell description, emptry string if none
            if '-KO' not in row['Treatment']:
                if row['KO'] == 'WT':
                    ko = ''
                else:
                    ko = f" {row['KO']}-KO"
            else:
                ko = ''
            # nice strings
            treatstr = f"{row['Ctrl samp']} ➤ {row['Treat samp']}"
            if row['Timepoint'] == 'endpoints':
                comp_description = f"{treatstr} – Effect of {row['Treatment']} in {row['Cell line']}{ko} cells"
            else:
                comp_description = f"{treatstr} – {row['Treatment']} in {row['Cell line']}{ko} cells"
            comp_options.append(
                {'value':idx, 'label':comp_description}
            )

        return get_xy_selectors(comp_options)


    def get_xyscores_genes(xk, yk, selected_analysis_type='drz'):
        """return score series from Dataset with unified indexes"""
        score_fdr = data_set.get_score_fdr(selected_analysis_type)

        # if selected_analysis_type not in comparisons.loc[xk, 'Available analyses']:
        #     LOG.debug(f'Unavailable analysis type {str(selected_analysis_type)}, comp {xk}')
        #     return dash.no_update, [f'Analysis type, {analysis_label}, not available for this experiment']

        x, y = [score_fdr['score'][k].dropna() for k in (xk, yk)]
        shared_genes = x.index.intersection(y.index)
        x, y = [xy[shared_genes] for xy in (x,y)]
        return x, y


    # update the keys stored in the data store, and gene selection options
    @app.callback(
        Output('selected-comps', 'data'),
        Output('gene-selector', 'options'),
        Input('x-selector', 'value'),
        Input('y-selector', 'value'),
    )
    def update_selection(xk, yk):
        LOG.debug(f'update_selection({xk}, {yk})')
        if (not xk) or (not yk):
            LOG.debug('not updating')
            raise PreventUpdate
        x, _ = get_xyscores_genes(xk, yk, 'drz')
        return {'xk':xk, 'yk':yk}, get_lab_val(x.index)

    # get selections from the table
    @app.callback(
        Output('gene-selector', 'value'),
        Output('previous-rows', 'data'),
        Input('table', 'selected_rows'),
        State('previous-rows', 'data'),
        State('gene-selector', 'value'),
        State('table', 'data'),
    )
    def select_gene_with_table(selected_rows, prev_rows, selected_genes, table_data):
        """Whenever a row is selected or deselected from the table
        add or remove that gene from the selected genes. Additions
        made to or from the dropdown do not effect the table"""
        LOG.debug('select_gene_with_table()')
        selected_rows = set(selected_rows)
        row_selected_genes = set([table_data[i]['Gene'] for i in selected_rows])
        prev_rows = set(prev_rows)
        new_gene = row_selected_genes.difference(prev_rows)
        # Writing this like it's possible to do multiple operations at once
        output_selection = selected_genes+list(new_gene)
        removed_genes = prev_rows.difference(row_selected_genes)
        output_selection = [g for g in output_selection if g not in removed_genes]
        if not new_gene and not removed_genes:
            LOG.debug('No update to selected genes')
            raise PreventUpdate
        LOG.debug('select_gene_by_table -> '+f"{output_selection}, {list(row_selected_genes)}")
        return output_selection, list(row_selected_genes)


    # update the chart and table
    @app.callback(
        Output('graph', 'figure'),
        Output('table-div', 'children'),
        #Output('missing-analysis', 'children'),
        Input('selected-comps', 'data'),
        Input('gene-selector', 'value'),

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
        title = (f"{treatlines[0]} Vs {treatlines[1]}<br>{row['Library']}" +
                 f" library, experiment ID {row['Experiment ID']}")
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
            )
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


if __name__ == '__main__':
    def main():
        args = sys.argv
        if (len(args) == 1) or (args[1] in ('-h', '--help')):
            print('usage: comparison_maker.py source_dir port [debug]\n    Any value in the debug position means True.')
        source = sys.argv[1]
        port = sys.argv[2]
        if len(sys.argv) > 3:
            debug = True
        else:
            debug = False

        app = dash.Dash(__name__)

        if debug:
            LOG.setLevel(logging.DEBUG)

        source_directory = pathlib.Path(source)

        data_set = DataSet(source_directory)

        initiate(app, data_set)
        app.run_server(debug=debug, host='0.0.0.0', port=int(port))

    main()