from dash_table import DataTable
from dash_table.Format import Format, Scheme
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
from dash_html_components import Div

#from https://sashamaps.net/docs/resources/20-colors/, 99%,
#todo test all these colours
colours = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990',
           '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#333333', ]
get_lab_val = lambda arr: [{'label': v, 'value':v} for v in arr]
styles = {'selector':{'display': 'inline-block', 'border-collapse': 'separate', 'border-spacing': '15px 50px'},
          'hidden':{'display':'none'}}


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


def get_reg_stat_selectors(app=None):
    """Return radio selectors, for selecting stats to be used in plotting.

    Registers a function that handles the score selection. Needs to have
    output values "score-selector" & "fdr-selector" captured by the fig
    update callback.

    If app is None, you'll need to reimpliment all that."""

    if app is not None:
        @app.callback(
            [Output('score-selector', 'value'),
             Output('fdr-selector', 'value'),
             Output('mixed-div', 'style')],

            [Input('analysis-selector', 'value')],

            [State('score-selector', 'value'),
             State('fdr-selector', 'value')]
        )
        def select_stats_primary(selection, curr_score, curr_fdr):
            if selection == 'mixed':
                return curr_score, curr_fdr, styles['selector']
            else:
                return selection, selection, styles['hidden']

    return Div([
        Div([
            html.Label('Analysis type_____', htmlFor='analysis-selector'),
            dcc.RadioItems(
                id='analysis-selector',
                options=[
                    {'label':'DrugZ', 'value':'drz'},
                    {'label':'MAGeCK',  'value':'mag'},
                    {'label':'Mixed...', 'value':'mixed'}
                ],
                value='drz',
            )
        ], style=styles['selector']),

        # This Div is hidden unless "Mixed" is chosen from Div above.
        Div([
            Div([
                html.Label('Effect size_____', htmlFor='score-selector'),
                dcc.RadioItems(
                    id='score-selector',
                    options=[{'label':'NormZ', 'value':'drz'},
                             {'label':'LFC',  'value':'mag'}],
                    value='drz',
                ),
            ], style=styles['selector']),

            Div([
                html.Label('FDR source', htmlFor='fdr-selector', ),
                dcc.RadioItems(
                    id='fdr-selector',
                    options=[
                        {'label':'DrugZ', 'value':'drz'},
                        {'label':'MAGeCK',  'value':'mag'}
                    ],
                    value='drz',
                )], style=styles['selector'])
        ],  id='mixed-div', style=styles['hidden'])
    ])

def get_data_source_selector(data_sources) -> Div:
    """A Div with a checklist that will be populated data_sources and a
    paragraph for reporting missing datasets

    IDs: data-source-selector, missing-datasets"""
    return Div([
        dcc.Checklist(
            id='data-source-selector',
            options=get_lab_val(data_sources),
            value=data_sources, # set all selected by default
            labelStyle={'display':'inline-block'}
        ),

        html.P([''], id='missing-datasets')
    ])



# **DATATABLE**
def create_datatable(data_df=None, columns_if_no_df=None):
    #todo make numfmt work
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