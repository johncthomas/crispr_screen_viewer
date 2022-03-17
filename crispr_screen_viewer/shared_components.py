from dash.dependencies import Input, Output, State
import logging
from dash.dash_table import DataTable
from dash.dash_table.Format import Format, Scheme
from dash import dcc, html
Div = html.Div
from .functions_etc import DataSet

logging.basicConfig()
LOG = logging.getLogger('screen_viewers')

#from https://sashamaps.net/docs/resources/20-colors/, 99%,
from typing import List

colours = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990',
           '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#333333', ]
get_lab_val = lambda arr: [{'label': v, 'value':v} for v in arr]
styles = {'selector':{'display': 'inline-block', 'border-collapse': 'separate',
                      'border-spacing': '15px 50px', 'margin-bottom': '15px'},
          'hidden':{'display':'none'}}
big_text_style = {
    'font-family': 'Arial, Helvetica, sans-serif',
    'font-size': '25px',
    'letter-spacing': '-0.4px',
    'word-spacing': '-0.4px',
    'font-weight': '700'}

def get_annotation_dicts(xs,ys,txts, annote_kw=None):
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

def get_data_source_selector(data_set, id_prefix='') -> List[dcc.Checklist]:
    """A Div with a checklist that will be populated data_sources and a
    paragraph for reporting missing datasets

    IDs: data-source-selector, missing-datasets"""

    if id_prefix:
        id_prefix = id_prefix+'-'

    return [
        html.Label('Select data sources:', htmlFor=id_prefix+'data-source-selector'),
        dcc.Checklist(
            id=id_prefix+'data-source-selector',
            options=get_lab_val(data_set.data_sources),
            value=data_set.data_sources, # set all selected by default
            labelStyle={'display':'inline-block'}
        ),
        html.P([''], id=id_prefix+'missing-datasets'),
    ]



# **DATATABLE**
def create_datatable(data_df=None, columns_if_no_df=None):
    #todo make numfmt work
    # to do that, it should be set per column
    #numfmt = Format(precision=3, scheme=Scheme.decimal_or_exponent)
    #formatted = Format()
    #numfmt = formatted.precision(3)

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
            columns=[{'name':c, 'id':c, } for c in columns_if_no_df],
        )

    return DataTable(
        id='table',
        # 'format':get_fmt(c)
        columns=[{'name':c, 'id':c, } for c in data_df.columns],
        data=data_df.to_dict('records'),
        export_format='csv',
        style_data_conditional=cond_fmts,
    )



external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']