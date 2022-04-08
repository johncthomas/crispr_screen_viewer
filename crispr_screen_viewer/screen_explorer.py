#!
import logging

import pandas as pd
import numpy as np

from dash import dash, dcc, html, Input, Output, State, dash_table

from dash.exceptions import PreventUpdate

import plotly.graph_objs as go

import pathlib, os
from dash.dependencies import Input, Output, State
from typing import Collection, Union, Dict
from crispr_screen_viewer.functions_etc import DataSet
from crispr_screen_viewer.shared_components import (
    get_lab_val,
    get_reg_stat_selectors,
    get_annotation_dicts,
    big_text_style,
    LOG
)

from crispr_screen_viewer.functions_etc import (
    DataSet,
    parse_expid
)

import copy

Div = html.Div

def init_msgv(app, data_set:DataSet, public_version=False) -> Div:
    """Source directory should contain the relevant info: metadata.csv,
    screen_analyses and expyaml directories."""

    comparisons = data_set.comparisons


    if public_version:
        source_display = 'none'
    else:
        source_display = 'block'

    # ****Componenets****
    volcano = dcc.Graph(
        id='volcano0',
        config={
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d',
                                       'autoScale2d', 'resetScale2d'],
            'editable':True
        },
        style={'height': '800px', 'padding':'0px'},
        figure={'layout':{'clickmode':'event+select'}}
    )

    # make the dropdowns for filtering
    filter_dropdowns = []
    filter_cols = ['Treatment', 'Experiment ID', 'Cell line', 'Library', 'Source']
    for col in filter_cols:
        filter_dropdowns.append(
            html.Div(
                children=[dcc.Dropdown(
                    id='se-'+col,
                    placeholder='Filter by '+col,
                    multi=True,
                    style={'height':'80px', 'width':'250px'},
                    value=[],
                    options=[{'label':v, 'value':v} for v in sorted(comparisons[col].unique())]
                )],
                style={'display':'inline-block'})
        )

    # construct data table with selected columns from the metadata
    if not public_version:
        metadata_columns = ['Comparison ID', 'Experiment ID', 'Treatment', 'Dose', 'Growth inhibition %',
                           'Days grown', 'Cell line', 'KO',  'Library', 'Source',
                           'Available analyses']
        table_data = comparisons.to_dict('records')
    else:
        metadata_columns = [ 'Treatment', 'KO', 'Cell line', 'Library', 'Dose', 'DOI', 'Citation']
        table_data = comparisons

        # get the experiment data for these comparisons
        for k in ('DOI', 'Citation'):
            table_data.loc[:, k] = table_data['Experiment ID'].apply(
                lambda exp: data_set.experiments_metadata.loc[exp, k]
            )

        # had weird issue of new columns not getting updated data without doing
        # this here.
        table_data = table_data.to_dict('records')

    table_of_comparisons = dash_table.DataTable(
        id='comparisons-table',
        columns=[{'name':c, 'id':c} for c in metadata_columns],
        data=table_data,
        sort_action='native',
        sort_mode="multi",
        selected_rows=[],
        row_selectable='single',
    )

    gene_selector = [
        html.Label('Select genes:',
                   htmlFor='se-gene-dropdown',
                   style=big_text_style),
        dcc.Dropdown(
            id='se-gene-dropdown',
            placeholder='Select genes by name',
            multi=True,
            #style={'height':'100px', },
            value=[],
            options=[]),
    ]

    # Tabs, one for selecting a comparison, one for viewing the graph
    #   maybe another for a datatable of the selected comps
    tabs = dcc.Tabs([
        # comparison selector
        dcc.Tab([
            html.P([(
                'Sele'+'ct an experiment from the table below. ' # broke up to stop pycharm treating it like SQL...
                'Use dropdowns to filter rows to find comparisons '
                'with specific treatments, etc.'
            )], style={'margin-top': '15px'}),
            Div(filter_dropdowns, style={'margin-bottom': '15px', }),
            Div([table_of_comparisons])
        ], label='Select experiment/comparison', value='comparison-selector'),
        # graph
        dcc.Tab([
            Div([volcano]),
            Div(get_reg_stat_selectors(app, id_prefix='se'), style={'display':source_display}),
            Div(html.P('', id='se-missing-analyses', style={"background-color":"#e60000", 'color':'white'}) ),
            Div(gene_selector)
        ], label='Graph', value='graph'),
    ], value='comparison-selector')

    # ****LAYOUT**** to be returned
    se_layout = Div([
        html.H1("Screens explorer"),
        tabs,
        Div([html.P(id='debug')], ),
        dcc.Store('volcano-data', 'session',
                  data={'x':[], 'y':[], 'genes':[]})
    ])


    # ***CALLBACKS***
    # callbacks. Filter datatable, select comparison
    @app.callback(
        Output('comparisons-table', 'data'),
        Output('comparisons-table', 'selected_rows'),
        # Filters from metadata columns, additional filtrs go after
        [Input('se-'+cid, 'value') for cid in filter_cols],
        State('comparisons-table', 'selected_rows'),
        State('comparisons-table', 'data'),
        # additional states will require change to zip line.
    )
    def filter_datatable(*filters):
        LOG.debug(f'CALLBACK: filter_datable()')
        selected_row = filters[-2]
        table_data   = filters[-1]

        # filter the table
        filters = filters[:-2]
        filtered_table = comparisons.copy()
        # go through each of the filters, apply them to the table_of_comparisons
        for col, filtr in zip(filter_cols, filters):
            if not filtr:
                continue
            # select the rows that contain filtered values
            filtered_table = filtered_table.loc[
                filtered_table[col].apply(lambda val: val in filtr)
            ]

        # If the selected row still exists in the table, maintain that selection
        #   otherwise deselect the row so we don't show some random graph
        #   (happily, deselecting prevents update so we keep the same graph).
        # Get the pre-filtered compID
        if selected_row:
            correct_compid = table_data[selected_row[0]]['Comparison ID']
            # new, if it's still in the table
            if correct_compid in filtered_table.index:
                # it's a list, cus selected_rows in datatable could be multiple
                new_selected_row = [filtered_table.index.get_loc(correct_compid)]
            else:
                new_selected_row = []
        else:
            new_selected_row = []


        return (filtered_table.to_dict('records'),
                new_selected_row)


    # render volcano for the selected comparison.
    @app.callback(
        Output('volcano-data', 'data'),
        Output('se-gene-dropdown', 'options'),
        # used for printing error message when missing score/fdr types chosen
        Output('se-missing-analyses', 'children'),

        Input('comparisons-table', 'selected_rows'),
        Input('se-score-selector', 'value'),
        Input('se-fdr-selector', 'value'),

        State('comparisons-table', 'data'),
    )
    def select_experiment_stats(selected_row, score_type, fdr_type,  table_data):

        args_for_printing = {k:v for k, v in zip(
            'selected_row, score_type, fdr_type,  table_data'.split(', '),
            [selected_row, score_type, fdr_type,  type(table_data)]
        )}
        LOG.debug(f'CALLBACK: select_experiment_stats with {args_for_printing}')
        if not selected_row:
            raise PreventUpdate

        # get the comparison ID, select the relevant data from dataset
        compid = table_data[selected_row[0]]['Comparison ID']

        available_analyses = data_set.comparisons.loc[compid, 'Available analyses']
        analyses_are_available = all(
            [ ans_type in available_analyses
             for ans_type in (score_type, fdr_type) ]
        )
        # if an analysis is missing, prevent the other things from trying to update
        if not analyses_are_available:
            availablility_text = ('  Analyses type not available for this experiment.\n'
                                  f"Available type(s): {[data_set.analysis_labels[k] for k in available_analyses]}\n"
                                  "Graph has not updated.")
            volcano_data = dash.no_update
            gene_options = dash.no_update

        else:
            availablility_text = ''

            # get x, y and genes values
            score_fdr = data_set.get_score_fdr(score_type, fdr_type)
            score, fdr = [score_fdr[k][compid].dropna() for k in ('score', 'fdr')]

            # # not supporting mixing score types
            # # combining analyses may result in different gene lists, so use interection
            # if score_type != fdr_type:
            #     unified_index = score.index.intersection(fdr.index)
            #     score,fdr = [xy.reindex(unified_index) for xy in (score,fdr)]


            volcano_data = {'x': score, 'y': fdr, 'genes': score.index}
            gene_options = get_lab_val(score.index)

        LOG.debug(f'End of select_experiment_stats with:')
        LOG.debug('     datatable:  '+'\n'.join([f"{k}={volcano_data[k].head()}" for k in ('x', 'y')]))

        return (
            volcano_data,
            gene_options,
            availablility_text,
        )


    @app.callback(
        Output('volcano0', 'figure'),

        Input('se-gene-dropdown', 'value'),
        Input('volcano-data', 'data'),

        State('comparisons-table', 'data'),
        State('comparisons-table', 'selected_rows'),
        State('se-score-selector', 'value'),
    )
    def render_volcano(selected_genes, xy_genes,
                       table_data, selected_row, score_type):
        LOG.debug("CALLBACK: render_volcano")

        # debug device
        def count_upper():
            n = 0
            while True:
                n+=1
                yield n
        counter = count_upper()
        LOG.debug(next(counter))


        if not selected_row:
            raise PreventUpdate

        # get values and create the figure
        x, fdr, genes = [xy_genes[k] for k in ('x', 'y', 'genes')]
        # dcc.Store converts this to a list...
        x, fdr = [pd.Series(xy, index=genes) for xy in (x,fdr)]
        y = fdr.apply(lambda _x: -np.log10(_x))

        fig = go.Figure(
            data=go.Scattergl(
                x=x,
                y=y,
                mode='markers',
                customdata=fdr,
                text=genes,
                hovertemplate= ("<b>%{text}</b><br>" +
                                f"{data_set.score_labels[score_type]}"+": %{x:.2f}<br>" +
                                "FDR: %{customdata:.2e}")
            ),
        )
        LOG.debug(f"{(len(x), len(y))}")
        LOG.debug(next(counter))
        # add some titles
        row = table_data[selected_row[0]]
        if '-KO' not in row['Treatment']:
            if row['KO'] == 'WT':
                ko = ''
            else:
                ko = f" {row['KO']}-KO"
        else:
             ko = ''
        title = (f"<b>Effect of {row['Treatment']} in {row['Cell line']}{ko} cells<br></b>"
                 f"{row['Library']} library, experiment ID {row['Experiment ID']}")

        fig.update_layout(
            title=title,
            xaxis_title=data_set.score_labels[score_type],
            yaxis_title='-Log10(FDR)',
        )
        LOG.debug(next(counter))
        # Add annotations for the selected genes
        new_annotations = get_annotation_dicts(x[selected_genes], y[selected_genes], selected_genes)
        for anot in new_annotations:
            fig.add_annotation(
                **anot
            )
        LOG.debug(next(counter))
        LOG.debug(str(fig))
        return fig

    return se_layout



# if __name__ == '__main__':
#     # expd_fn, port,
#     # fdr filter, analysis types
#     parser = argparse.ArgumentParser(description='')
#     parser.add_argument(
#         '-d', '--data-version',
#         dest='data_version',
#         required=True,
#         help="Name of the directory within app_data that contains the data from screens."
#     )
#     parser.add_argument(
#         '-p', '--port',
#         required=True,
#         help='Port used to serve the charts'
#     )
#     parser.add_argument(
#         '--debug', action='store_true',
#         help='Launch app in debug mode'
#     )
#     args = parser.parse_args()
#     launch_volcano(os.path.join('app_data', args.data_version), args.port, args.debug)
