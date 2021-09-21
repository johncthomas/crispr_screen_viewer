import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
Div = html.Div
import plotly.graph_objs as go
import pathlib, os
from dash.dependencies import Input, Output, State
from typing import Collection, Union, Dict
from crispr_screen_viewer.functions_etc import index_of_true, DataSet
from crispr_screen_viewer.shared_components import (
    create_datatable,
    get_data_source_selector,
    get_lab_val,
    get_reg_stat_selectors,
    colours,
    big_text_style,
)

border_style = {'border': '4px solid #3DD178',
                'border-radius': '14px'}

#css
"""#funBorder {
border: 4px solid #3DD178;
border-radius: 14px;
}"""


# *updates
# 1.0.3 data_version is an argument
# 1.0.4 stat selector: mixed is now a separate selection
# 1.1 data source selection
msgv_version = '1.1.0'
def launch(source_directory:Union[str, os.PathLike], port, debug):
    """A Dash app for showing results from screens for specific screens."""
    print(source_directory)
    data_set = DataSet(source_directory)

    def get_colour_map(list_of_things):
        return {thing:colours[i%len(colours)] for i, thing in enumerate(list_of_things)}

    app = dash.Dash(__name__, external_stylesheets= ['https://codepen.io/chriddyp/pen/bWLwgP.css'])

    # Graph and table
    graph = dcc.Graph(
        id='gene-violins',
        figure=go.Figure(),
        style={'height':'800px'}
    )

    table = Div([create_datatable(columns_if_no_df=data_set.metadata.columns)],
                id='table-div', className="u-full-width", style={'margin-bottom':'30px'})

    # select genes by name, and comparisons FDR in selected samples
    gene_selector = Div([
        html.P('Select genes:', style=big_text_style),
        dcc.Dropdown('gene-selector', placeholder='Select genes', value=[],
                     options=get_lab_val(data_set.genes), multi=True),

    ], style={'margin-bottom': '15px'})

    # the graph/table object and the data prefilters
    graph_and_data_selection_div = Div([
        html.H1("Multi-Screen Gene Viewer"),
        html.P("Select your gene(s) of interest. Comparisons that show significant results "
                "(below adjustable FDR) for at least one selected gene are shown below. "
                "Box plots give the overall distribution of scores, and stars show specific genes."),
        Div([
            dcc.Checklist(
                'show-boxplots',
                options=[{'label':'Show boxplots', 'value':'show-boxplots'}],
                value=['show-boxplots']
            ),
        ], style={'display':'none'}),

        gene_selector,
        # This Tabs will show either graph or datatable
        Div([dcc.Tabs(id='output-tabs', value='graph',
                  children=[
                      dcc.Tab([graph], label='Graph', value='graph'),
                      dcc.Tab([table], label='Table', value='table')
                  ])
             ]),
        Div([
            Div(get_reg_stat_selectors(app),
                style={'display':'inline-block', 'width':'170px','vertical-align':'top'}),
            Div(get_data_source_selector(data_set.data_sources), style={'display':'inline-block'})
        ])
    ])

    # this is also used for output of one function, so is defined once here
    order_by_categories = ['Mean score', 'Treatment', 'Experiment ID']

    control_bar = Div([
        # Div([
        #     html.P('Use controls below to filter which comparisons are shown. ')
        # ]),
        Div([
            html.Label('Maximum FDR:', htmlFor='fdr-threshold'),
            dcc.Input('fdr-threshold', type='number', min=0, max=2, step=0.01, value=0.2),
        ], style={ 'display':'inline-block', 'width':'135px'}),
        Div([
            html.Label('Order by:', htmlFor='order-by'),
            dcc.Dropdown('order-by', value=order_by_categories[0],
                         options=get_lab_val(order_by_categories)),
        ], style={'width':'150px', 'display':'inline-block', 'vertical-align':'top'}),
        Div([
            html.Label('Colour by:', htmlFor='color-by'),
            dcc.Dropdown(
                'color-by', value='Treatment',
                options=get_lab_val(['Treatment', 'Cell line', 'Experiment ID', 'Library', 'KO'])
            ),
        ], style={'width':'150px', 'display':'inline-block', 'vertical-align':'top'}),
    ])



    # make the dropdowns for filtering
    filter_dropdowns = []
    filter_cols = ['Treatment', 'Experiment ID', 'KO', 'Cell line', 'Library', 'Source']
    for col in filter_cols:
        filter_dropdowns.append(
            html.Div([dcc.Dropdown(
                id=col,
                placeholder='Filter by '+col,
                multi=True,
                style={'height':'80px', 'width':'220px'},
                value=[],
                options=[{'label':v, 'value':v} for v in sorted(data_set.metadata[col].unique())]
            )], style={'display':'inline-block'})
        )

    # put it all together
    app.layout = Div([
        graph_and_data_selection_div,
        html.Br(),
        control_bar,
        html.Br(),
        Div(filter_dropdowns, ),
        #html.Br(),
        #table,
    ])


    # Define callback to update graph
    @app.callback(
        [Output('gene-violins', 'figure'),
         Output('table-div', 'children'),
         Output('order-by', 'options'),
         Output('missing-datasets', 'children')],

        [Input('score-selector', 'value'),
         Input('fdr-selector', 'value'),
         Input('data-source-selector', 'value'),

         Input('gene-selector', 'value'),
         Input('show-boxplots', 'value'),
         Input('fdr-threshold', 'value'),

         Input('order-by', 'value'),
         Input('color-by', 'value'),
         ] + [Input(cid, 'value') for cid in filter_cols]
    )
    def update_figure(score_type, fdr_type, selected_data_sources,
                      selected_genes, show_boxplots, fdr_thresh,
                      order_by, color_by, *filters):

        data_tabs = data_set.get_score_fdr(score_type, fdr_type, selected_data_sources)

        # Identify data sources that are unavailable for selected analyses
        available_sources = data_set.metadata.loc[data_tabs['fdr'].columns, 'Source'].unique()
        missing_data_sources = [d for d in selected_data_sources if d not in available_sources]
        if missing_data_sources:
            missing_data_sources = "(unavailable with current analysis type: "+', '.join(missing_data_sources)+')'
        else:
            missing_data_sources = '(All selections available for analysis type)'

        # get boolean masks for which comparisons to include in the charts
        # first get comparisons filtered by metadata filters
        comparison_mask = pd.Series(True, index=data_set.metadata.index)
        for filter_id, values in zip(filter_cols, filters):
            if values:
                comparison_mask = comparison_mask & data_set.metadata[filter_id].isin(values)

        # then by FDR threshold
        fdr_mask = (data_tabs['fdr'].loc[selected_genes] < fdr_thresh).any()
        filtered_scores = data_tabs['score'].loc[selected_genes, (fdr_mask & comparison_mask)]

        # determine order of plots
        if order_by in selected_genes:
            trace_order = filtered_scores.loc[order_by].sort_values().index.values
        elif order_by == 'Mean score':
            trace_order = filtered_scores.mean().sort_values().index.values
        elif order_by in order_by_categories[1:]:
            # subset the metadata to included comps and then sort by the order_by
            trace_order = data_set.metadata.loc[filtered_scores.columns, order_by].sort_values().index
        else: # this shouldn't happen, though maybe I should have a "whatever" order option
            trace_order = filtered_scores.columns.values

        # assemble the figure
        fig = go.Figure()
        # plot the gene scatters
        for gn in selected_genes:
            fig.add_trace(
                go.Scatter(
                    x=trace_order, y=filtered_scores.loc[gn, trace_order],
                    mode='markers', name=gn,
                    marker={'size': 15, 'line':{'width':2, 'color':'DarkSlateGrey'}, 'symbol':'hexagram'}),

            )
        # add the boxplot traces if required
        if show_boxplots:
            colour_map = get_colour_map(data_set.metadata.loc[:, color_by].unique())
            for col in trace_order:
                ys = data_tabs['score'][col]
                xs = ys[:]
                xs[:] = col
                color_group = data_set.metadata.loc[col, color_by]
                fig.add_trace(
                    go.Box(x=xs, y=ys, name=color_group, boxpoints=False,
                           line=dict(color=colour_map[color_group]))
                )
        # labels
        fig.update_yaxes(title_text=data_set.score_labels[score_type])

        # create the DataTable
        selected_fdr = data_tabs['fdr'].loc[filtered_scores.index, trace_order]
        filtered_scores = filtered_scores.reindex(columns=trace_order)
        filtered_scores.index = filtered_scores.index.map(lambda x: x+' (Effect size)')
        selected_fdr.index = selected_fdr.index.map(lambda x: x+' (FDR)')
        selected_stats = pd.concat([filtered_scores, selected_fdr], sort=False).T

        selected_stats = selected_stats.applymap(lambda n: f"{n:.3}")
        selected_stats.insert(0, 'Comparison', selected_stats.index)
        cols_oi = ['Experiment ID', 'Treatment', 'Dose', 'KO', 'Growth inhibition %', 'Days grown',
                   'Cell line', 'Library']
        selected_metadata = data_set.metadata.loc[selected_stats.index, cols_oi]
        data_table_data = pd.concat([selected_stats, selected_metadata], axis=1)
        dtable = create_datatable(data_table_data)

        sort_by_opts = get_lab_val(order_by_categories+ selected_genes)

        return fig, dtable, sort_by_opts, [missing_data_sources]

    app.run_server(host='0.0.0.0', port=port, debug=debug)