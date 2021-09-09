import pandas as pd
import numpy as np

import dash
import dash_table
from dash.exceptions import PreventUpdate

import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import pathlib, os
from dash.dependencies import Input, Output, State
from typing import Collection, Union, Dict
from crispr_screen_viewer.functions_etc import DataSet
from crispr_screen_viewer.shared_components import (external_stylesheets,
                               get_lab_val,
                               get_reg_stat_selectors,
                               get_annotation_dicts)

#todo deal with missing analysis type for selected comparison

def launch(source_directory:Union[str, os.PathLike], port, debug=False):
    """Source directory should contain the relevant info: metadata.csv,
    screen_analyses and expyaml directories."""
    source_directory = pathlib.Path(source_directory)

    data_set = DataSet(source_directory)
    metadata = data_set.metadata

    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


    # ****Componenets****
    graph = dcc.Graph(
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
    filter_cols = ['Treatment', 'Experiment ID', 'Cell line', 'Library', 'Source'] # list of ids (lowercase, underscored), for populate table callback
    for col in filter_cols:
        filter_dropdowns.append(
            html.Div([dcc.Dropdown(
                id=col,
                placeholder='Filter by '+col,
                multi=True,
                style={'height':'80px', 'width':'250px'},
                value=[],
                options=[{'label':v, 'value':v} for v in sorted(metadata[col].unique())]
            )], style={'display':'inline-block'})
        )

    metadata_columns = ['Comparison ID', 'Experiment ID', 'Treatment', 'Dose', 'Growth inhibition %',
                       'Days grown', 'Cell line', 'KO',  'Library', 'Source',
                       'Available analyses']
    table_of_comparisons = dash_table.DataTable(
        id='comparisons-table',
        columns=[{'name':c, 'id':c} for c in metadata_columns],
        data=metadata.to_dict('records'),
        sort_action='native',
        sort_mode="multi",
        selected_rows=[],
        row_selectable='single',
    )

    gene_selector = dcc.Dropdown(
        id='gene-dropdown',
        placeholder='Select genes by name',
        multi=True,
        #style={'height':'100px', },
        value=[],
        options=[]
    )

    # ****LAYOUT****
    Div = html.Div
    app.layout = Div([
        get_reg_stat_selectors(app),
        Div(html.P('', id='missing-analyses', style={"background-color":"#e60000", 'color':'white'}) ),
        #get_data_source_selector(data_set.data_sources),
        Div(graph),
        Div(gene_selector),
        Div(filter_dropdowns),
        Div([table_of_comparisons]),
        Div([html.P(id='debug')], ),
        dcc.Store('xy-genes', 'session',
                  data={'x':[], 'y':[], 'genes':[]})
    ])


    # ***CALLBACKS***
    # callbacks. Filter datatable, select comparison
    @app.callback(
        Output('comparisons-table', 'data'),
        Output('comparisons-table', 'selected_rows'),
        # Filters from metadata columns, additional filtrs go after
        [Input(cid, 'value') for cid in filter_cols],
        State('comparisons-table', 'selected_rows'),
        State('comparisons-table', 'data'),
        # additional states will require change to zip line.
    )
    def filter_datatable(*filters):
        selected_row = filters[-2]
        table_data   = filters[-1]

        # filter the table
        filters = filters[:-2]
        filtered_table = metadata.copy()
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
        #   (ideally, we'd always maintain selection, but that requires some
        #   rejiggering)
        # get the pre-filtered compID
        if selected_row:
            correct_compid = table_data[selected_row[0]]['Comparison ID']
            # new, if it's still in the table
            if correct_compid in filtered_table.index:
                # it's a list, cus data table could multiple
                new_selected_row = [filtered_table.index.get_loc(correct_compid)]
            else:
                new_selected_row = []
        else:
            new_selected_row = []


        return (filtered_table.to_dict('records'),
                new_selected_row)


    @app.callback(
        Output('xy-genes', 'data'),
        Output('gene-dropdown', 'options'),
        # used for printing error message when missing score/fdr types chosen
        Output('missing-analyses', 'children'),

        Input('comparisons-table', 'selected_rows'),
        Input('score-selector', 'value'),
        Input('fdr-selector', 'value'),

        State('comparisons-table', 'data'),
    )
    def select_experiment_stats(selected_row, score_type, fdr_type,  table_data):

        if not selected_row:
            raise PreventUpdate

        # get the comparison ID, select the relevant data from dataset
        compid = table_data[selected_row[0]]['Comparison ID']

        available_analyses = data_set.metadata.loc[compid, 'Available analyses']
        analyses_are_available = all(
            [ ans_type in available_analyses
             for ans_type in (score_type, fdr_type) ]
        )
        # if an analysis is missing, prevent the other things from trying to update
        if not analyses_are_available:
            availablility_text = ('  Analyses type not available for this experiment.\n'
                                  f"Available type(s): {set([data_set.analysis_labels[k] for k in available_analyses])}\n"
                                  "Graph has not updated.")
            data = dash.no_update
            gene_options = dash.no_update

        else:
            availablility_text = ''

            # get x, y and genes values
            score_fdr = data_set.get_score_fdr(score_type, fdr_type)
            x, y = [score_fdr[k][compid].dropna() for k in ('score', 'fdr')]

            # combining analyses may result in different gene lists, so use interection
            if score_type != fdr_type:
                unified_index = x.index.intersection(y.index)
                x,y = [xy.reindex(unified_index) for xy in (x,y)]

            data = {'x': x, 'y': y, 'genes': x.index}
            gene_options = get_lab_val(x.index)

        return (
            data,
            gene_options,
            availablility_text,
        )


    @app.callback(
        Output('volcano0', 'figure'),

        Input('gene-dropdown', 'value'),
        Input('xy-genes', 'data'),

        State('comparisons-table', 'data'),
        State('comparisons-table', 'selected_rows'),
        State('score-selector', 'value'),
    )
    def render_volcano(selected_genes, xy_genes,
                       table_data, selected_row, score_type):

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

        # Add annotations for the selected genes
        new_annotations = get_annotation_dicts(x[selected_genes], y[selected_genes], selected_genes)
        for anot in new_annotations:
            fig.add_annotation(
                **anot
            )

        return fig

    app.run_server(debug=debug, host='0.0.0.0', port=port)


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
