
import dash
import dash_table
from dash.exceptions import PreventUpdate
from dash_table.Format import Format, Scheme
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
from typing import Union, Dict
import pandas as pd
from typing import Iterable


from crispr_screen_viewer.util import tabulate_score, load_mageck_tables


import os, logging

"""Functionality is selecting genes by box select on the graph or by inputting names into
the select genes by name box

Todo:
    Markdown giving instructions
    Formatting options
    Stop being able to select Gene as a sample
    Clean up. This has the legacy of multiple failed attempts at reactivity littering up the place."""



def get_annotation_dicts(df, xkey, ykey, annote_kw=None):
    """A list of dict defining annotations from specified columns
    of a dataframe."""
    annotations = []
    if annote_kw is None:
        annote_kw = {}
    for txt, (x, y) in df.loc[:, [xkey, ykey]].iterrows():
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

def get_formatted_col(hdr):
    """Returns dict(id=hdr, name='Nice Header')
    Nice headers are from a static dictionary"""
    # for table headers and axis labels
    nice_headers = dict(
        fdr='FDR',
        fdr_log10='Log10(FDR)',
        lfc='Log2(Fold Change)',
        score='JACKS score',
        jacks_score='JACKS score'
    )

    if hdr not in nice_headers.keys():
        return {'name':hdr, 'id':hdr}
    else:
        return {
            'id': hdr,
            'name': nice_headers[hdr]
        }

#
# def get_options_thing(lst):
#     """pass some kind of iterable, get [{'label':v, 'value':v}, ...]"""
#     return

LOG = logging.getLogger('volcano')
LOG.setLevel(logging.DEBUG)
# no index_col, is converted using .to_dict('rows') which doesn't keep the index.


def union_gt_Y(tables:Dict[str,pd.DataFrame],
               alpha:float,
               ykey='fdr_log10'):
    """Give dict of DF, get bool series, True where ykey exceeds alpha in at
    least one of the DF."""
    # get the starting mask of Falses. iloc[:, 0] returns a series even with multindex apparently
    tab = list(tables.values())[0]
    yfilter = tab.iloc[:, 0].apply(lambda x: False)

    # get the union of y > filter across all tables
    for k, tab in tables.items():
        # if 'jacks' in k.lower():
        #     continue
        nextfilter = (tab.loc[:, (slice(None, None), ykey)] > alpha).any(1)
        yfilter = yfilter | nextfilter
        #print(yfilter.sum(), tab.shape[0])

    return yfilter

#todo: what this actually needs is a pair of sample selectors, rather than a single list of all analyses.
def spawn_volcanoes(tables:Union[pd.DataFrame, Dict[str, pd.DataFrame]],
                    xy_keys=('lfc', 'fdr_log10'), filterYLessThan:float=None,
                    groupings:pd.Series=None):
    """A DataFrame containing all available data determined in the filtering step
    multiindexed by (exp, stat), or a dictionary of such DF.
    Keys of the dict will define the first level of sample selection

    filterYLessThan: remove genes that never go above a  given -log10(fdr)
        to limit the number of points rendered. Ignores samples with 'jacks' in the
        name.
    """
    WIDTH = "80%"
    XKEY, YKEY = xy_keys
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    # individual analyses accessable through the dict
    if type(tables) is pd.DataFrame:
        tables = {'-':tables}

    for tab in tables.values():
        # give each subtable the gene columns
        for samp in tab.columns.levels[0]:
            tab.loc[:, (samp, 'gene')] = tab.index

    # filter by fdr, will remove all items that never appear
    if filterYLessThan:
        yfilter = union_gt_Y(tables, filterYLessThan, xy_keys[1])
        for k, tab in tables.items():
            tables[k] = tab.loc[yfilter]




    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    available_analyses = list(tables.keys())
    gene_list = list(tables[available_analyses[0]].index)

    #################
    class COMPONENTS:
        pass
    COMPONENTS()
    ################

    graph_config = {'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d',
                                                   'resetScale2d'],
                        'editable':True}

    analysis_dropdown = dcc.Dropdown(
        id='analysis-dropdown0',
        options=[{'label':v,'value':v} for v in available_analyses],
        value=None if len(available_analyses) > 1 else available_analyses[0],
        placeholder='Select an experiment or control group'
    )

    sample_dropdown = dcc.Dropdown(
        id='samples-dropdown0',
        value=None,
        placeholder='Select a sample',
        options=[]
    )

    graph = dcc.Graph(
        id='volcano0',
        config=graph_config,
        style={'height': '800px', 'padding':'0px'},
        figure={'layout':{'clickmode':'event+select'}}
    )

    empty_datatable = dash_table.DataTable(
        id='table0',
        columns=[{'id':'Selected','name':'Selected'}]+[get_formatted_col(hdr) for hdr in
                                                       list(list(tables.values())[0].columns.levels[1]) ],
        data=[],
        sort_action='native',
        sort_mode="multi",
    )

    def get_gene_dropdown(value):
        return dcc.Dropdown(
            id='gene-dropdown',
            placeholder='Select genes by name',
            multi=True,
            style={'height':'100px', },
            value=value,
            options=[{'label':v, 'value':v} for v in gene_list]

        )

    markdown = dcc.Markdown(
        """Select control sample group in the top dropdown, then a treatment in the second dropdown.  
        To select genes type the name in the box below, or choose box select in the top right of the
        chart.  \nAfter box-selecting it switches back to zoom-mode, if you accidently zoom in, 
        double click on the chart to zoom out again."""
    )

    #################
    # Layout
    app.layout = html.Div(
        [#
            html.Div([analysis_dropdown, sample_dropdown, markdown, graph],
            style={'width':WIDTH, 'display':'inline-block'}),
            html.Div(id='genes-box', children=[get_gene_dropdown([])],), #style={'width':'250px'}),
            html.Div(id='selected-genes', style={'display': 'none'},
                     children=[]),
            html.Div(children=[empty_datatable], style={'width':'90%', 'display':'inline-block', 'float':'center'})
        ]
    )


    ###################
    # callbacks
    class CALLBACKS:
        pass

    # nonlocal list of selected genes
    previous_points = []

    #Whenever selection or sample changes update the chart and table
    @app.callback(
        # outputting the selected text on all volcs
        [Output('volcano0', 'figure'),
         Output('table0', 'data'),
         Output('selected-genes', 'children'),
         Output('genes-box', 'children')],

        # input will come from a single volc each time
        [Input('volcano0', 'selectedData'),
         Input('samples-dropdown0', 'value'),],
        [State('analysis-dropdown0', 'value')]

    )
    def plot_all_points_and_some_text(selectedData, sampleName, analysisName):
        """If a selection changes in one chart, replot each plot with only the text
        markers in the current selection"""
        #print(analysisName, sampleName)
        #print('plot callback...')
        nonlocal previous_points

        # when the graph is first initialised...
        if sampleName is None or analysisName is None:
            raise PreventUpdate

        # maintain selected genes when swapping between samples
        if selectedData is None:
            previous_points = []
            genes = []
        else:
            #print(selectedData)
            #((does genes and current selectn have to be different things))
            genes = [d['text'] for d in previous_points]
            new_select = selectedData['points']
            if previous_points != new_select:
                previous_points = new_select
                genes = [d['text'] for d in previous_points]



        # update table data

        sample_df = tables[analysisName][sampleName]
        data_df = sample_df.copy()
        #format these stats
        for k in ['lfc', 'fdr_log10', 'score', 'jacks_score']:
            if k not in sample_df.columns:
                continue
            data_df.loc[:, k] = data_df[k].apply(lambda s: f'{s:.4}')
        for k in ['fdr', 'p']:
            if k not in sample_df.columns:
                continue
            data_df.loc[:, k] = data_df[k].apply(lambda s: f'{float(s):.2e}' if len(str(s)) > 4 else str(s))
        data_df.loc[:, 'Selected'] = 'No'
        data_df.loc[genes, 'Selected'] = '!Yes'
        cols = list(data_df.columns)
        # put is selected first
        data_df = data_df.reindex(columns=[cols[-1]] + cols[:-1])

        # Annotations are part of the layout
        layout = go.Layout(annotations = get_annotation_dicts(sample_df.loc[genes], XKEY, YKEY))
        print(sample_df.columns)
        _graph = go.Scatter(x=sample_df[XKEY], y=sample_df[YKEY],text=sample_df.index, mode='markers')

        new_genedropdown = get_gene_dropdown(genes)

        # return figure props and table data
        return [
            # the figure values
            {'data':[_graph],
             'layout':layout,},

            # table data
            data_df.to_dict('rows'),

            # dropdown values to the hidden div
            genes,

            # new dropdown
            [new_genedropdown]
        ]

    # define available samples based on the selected analysis
    @app.callback(
        Output('samples-dropdown0', 'options'),
        [Input('analysis-dropdown0', 'value')]
    )
    def populate_sample_dropdown_cb(analysisName):
        if analysisName is None:
            raise PreventUpdate
        samps = tables[analysisName].columns.levels[0]
        #print(list(tables[analysisName].index))
        return [dict(label=v, value=v) for v in samps]



    # update the selected genes from the dropdown or from the on graph selection.
    @app.callback(
        Output('volcano0', 'selectedData'),
        [Input('gene-dropdown', "value")],
        [State('volcano0', 'selectedData'), State('samples-dropdown0', 'value'),
         State('analysis-dropdown0', 'value'),
         State('volcano0', 'figure'), State('selected-genes', 'children')]
    )
    def gene_selection_cb(dropdown_values, selectedData, selectedSample, selectedAnalysis,
                          volcData, selected_genes, ):
        """"""
        #print('dropdown callback')
        #print(new_gene)
        if selectedData is None:
            selectedData = {'points':[], 'range':{'x':[0,0], 'y':[0,0]}}
        #currently_selected_genes = [p['text'] for p in selectedData['points']]
        #new_selections = [gn for gn in dropdown_values if gn not in selected_genes]
        # if dropdown_values is None or len(new_selections) == 0:
        #     raise PreventUpdate

        if sorted(dropdown_values) == sorted(selected_genes):
            raise PreventUpdate

        table = tables[selectedAnalysis]

        new_selected_points = []
        for gn in dropdown_values:
            dati = volcData['data'][0]['text'].index(gn)
            xy = [table.loc[gn, (selectedSample, k)] for k in (XKEY, YKEY)]
            point = dict(
                curveNumber=0,
                pointNumber=dati,
                pointIndex=dati,
                x=xy[0],
                y=xy[1],
                text=gn,
            )
            new_selected_points.append(point)
        selectedData['points'] = new_selected_points
        return selectedData

    # app.run_server(debug=True, port=8001)
    # return app.server
    return app

if __name__ == '__main__':

    import yaml
    #_expd = yaml.safe_load(open('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/david_756-7^ddrV2/dav756-7.yaml'))
    #f = '/Users/johnc.thomas/Dropbox/crispr/screens_analysis/david_756-7^ddrV2/dav_756-7/take1-firstruns/jacks_median/files/dav_756-7.'
    _expd = yaml.safe_load(open('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/ramsay_759+/ram_759-80.3.repmap.yaml'))

    f = "/Users/johnc.thomas/Dropbox/crispr/pkg/dash_charts/ram_759-80/take3/mageck/files/ram_759-80."
    _tables = load_mageck_tables(f, list(_expd['controls'].keys())+['EXTRA'])



    _tables = {k+" (MAGECK)":tab for k, tab in _tables.items()}
    jf = '/Users/johnc.thomas/Dropbox/crispr/pkg/dash_charts/ram_759-80/take3/jacks_median/files/ram_759-80.'
    for _ctrl_grp in _expd['controls'].keys():
        _tables[_ctrl_grp + "(JACKS)"] = tabulate_score(jf + _ctrl_grp + '.')

    app = spawn_volcanoes(_tables)
    server = app.server
    app.run_server(debug=True)



# #paff = '/Users/johnc.thomas/Dropbox/crispr/screens_analysis/david_756-7^ddrV2/dav_756-7/take1-firstruns/'
# tabz = {}
# for group in ('D3', 'NT', 'pretreat'):
#     magtab = tabulate_mageck(
#         'data/dav_756-7/take1-firstruns/mageck/files/dav_756-7.'+group+'.'
#     )
#     shared_cols = ['LFC or score', 'fdr', 'fdr_log10']
#     magtab.columns = shared_cols
#     tabz[group+' (MAGeCK)'] = magtab
#     jaktab = tabulate_score(
#         'data/dav_756-7/take1-firstruns/jacks/files/dav_756-7.'+group+'.'
#     )
#     jaktab.columns = shared_cols
#     tabz[group+' (JACKS)'] = jaktab
#
# app = spawn_volcanoes(tabz)
# server = app.server
# app.run_server(debug=True)