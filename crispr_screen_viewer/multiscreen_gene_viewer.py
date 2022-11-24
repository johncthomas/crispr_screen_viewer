import itertools
import inspect
import typing

import pandas as pd
from dash import dash, dcc, html, Input, Output, State
from dash.exceptions import PreventUpdate
Div = html.Div
import plotly.graph_objs as go

import pathlib, os, logging

from typing import Collection, Union, Dict

from crispr_screen_viewer.shared_components import (
    create_datatable,
    get_lab_val,
    get_gene_dropdown_lab_val,
    #get_reg_stat_selectors,
    get_stat_source_selector,
    colours,
    big_text_style,
    timepoint_labels,
    LOG,
    spawn_filter_dropdowns,
)

from crispr_screen_viewer.functions_etc import (
    doi_to_link,
    datatable_column_dict,
    get_metadata_table_columns,
    get_selector_table_filter_keys,
    getfuncstr
)

border_style = {'border': '4px solid #3DD178',
                'border-radius': '14px'}

import dash_bootstrap_components as dbc

PAGE_ID = 'msgv'

# todo clustergram
# Maybe give it it's own page,  use the selector-tables to select comps, and gene box
# But for now it's going to be a tab inside of MSGV

# todo data selection
# Currently MSGV is a single callback function, needs to be split
#  - Select the data and outputs to a datastore
#  - Update the boxplots/clustergram
#    * 2 functions, each has input from datastore and selected tab
#    * Output is the below Tabs Div

# todo the rendering of the graphs
# 2 cols, one gene vs comp, other heatmap of pearsons between comps using the selected genes
# Where some comps don't have some genes, we do one row of graphs excluding comps, one excluding genes

# todo controls
# The boxplot specific controls aren't needed and should probably be excluded
# Lock selected graphs needs to be a thing

# todo Make MSGV work like the others with the comparison selection tables
#   with select by FDR as an additional option.
# Obviously this is a bigger deal and further down the pipe

# Crucially, lots of callbacks sharing inputs is *fine*. You just can't have multiple
#   callbacks that share outputs.
# So the structure here is selected data (comps and genes, tables?) are selected and
#   put into a store, this store is the trigger for updating the tabs (tabs don't also
#   don't update unless selected).

# datatab on every page? yes


def initiate(app, data_set, public=True) -> Div:
    """Register callbacks to app, generate layout"""

    # if public:
    #     source_display = 'none'
    # else:
    #     source_display = 'inline-block'

    def spawn_boxplot_graph():

        graph = dcc.Graph(
            id='msgv-gene-boxplots',
            figure=go.Figure(),
            style={'height':'800px', 'width':'1500px'}
        )

        @app.callback(
            Output('msgv-gene-boxplots', 'figure'),
            Output('msgv-order-by', 'options'),
            Output('ordered-comp-store', 'data'),

            Input('msgv-stat-source-selector', 'value'),
            Input('comp-store', 'data'),
            Input('msgv-tabs', 'value'),
            Input('msgv-order-by', 'value'),
            Input('msgv-color-by', 'value'),

            State('msgv-gene-selector', 'value'),
            State('msgv-fdr-threshold', 'value'),
        )
        def update_boxplot_figure(
                score_type,
                comps_dict:Dict[str, typing.List[str]],
                selected_tab,
                order_by,
                color_by,

                selected_genes,
                fdr_thresh
        ):
            if selected_tab != 'msgv-boxplot-tab':
                raise PreventUpdate



            comps = comps_dict['selected_comps']
            LOG.debug(f"{getfuncstr()}: {comps}")
            data_tabs = data_set.get_score_fdr(score_type, score_type, )
            filtered_scores = data_tabs['score'].loc[selected_genes, comps]

            # *assemble the figure*
            fig = go.Figure()

            if (len(comps) == 0) and selected_genes:
                fig.add_annotation(x=4, y=4,
                                   text=f"None of the selected genes have FDR <= {fdr_thresh}",
                                   showarrow=False, )

            # determine order of plots
            if order_by in selected_genes:
                ordered_comps = filtered_scores.loc[order_by].sort_values().index.values
            elif order_by == 'Mean score':
                ordered_comps = filtered_scores.mean().sort_values().index.values
            elif order_by in order_by_categories[1:]:
                # subset the metadata to included comps and then sort by the order_by
                ordered_comps = data_set.comparisons.loc[filtered_scores.columns, order_by].sort_values().index
            else:  # this shouldn't happen, though maybe I should have a "whatever" order option
                ordered_comps = filtered_scores.columns.values

            # x tick labels
            trace_numbers = pd.Series([str(i) for i in range(1, len(ordered_comps) + 1)])

            if len(ordered_comps):
                citations = data_set.comparisons.loc[ordered_comps, 'Citation']
                treats = data_set.comparisons.loc[ordered_comps, 'Treatment'].fillna('')
                x_tick_labels = trace_numbers.values + '. ' + citations.values + ', ' + treats.values
            else:
                x_tick_labels = trace_numbers.values

            # x_tick_labels = pd.Series(x_tick_labels, index=ordered_comps)

            for gn in selected_genes:
                fdrs = data_tabs['fdr'].loc[gn, ordered_comps]
                mrkrs: pd.Series = fdrs <= fdr_thresh
                mrkrs[mrkrs == True] = 'diamond'
                mrkrs[mrkrs == False] = 'square'

                fig.add_trace(
                    go.Scatter(
                        x=x_tick_labels,
                        y=filtered_scores.loc[gn, ordered_comps],
                        mode='markers', name=gn,
                        marker_symbol=mrkrs.values,
                        marker={'size': 15, 'line': {'width': 2, 'color': 'DarkSlateGrey'}},
                        customdata=fdrs.apply(lambda n: f'{float(f"{n:.3g}"):g}'),
                        hovertemplate=f"{gn}" + "<br>Score: %{y}<br>FDR: %{customdata}<extra></extra>"
                    ),

                )

            # add the boxplot traces
            included = set()
            # add a boxplot trace for each comparison
            for trace_i, comp in enumerate(ordered_comps):
                # these values define the boxplot
                ys = data_tabs['score'][comp]
                fdr = data_tabs['fdr'][comp]
                # This gives the X value for each y value used to create the boxplot, which
                #   is apparently required? I guess this is approximating Tidy formated data?
                boxlabels = pd.Series(x_tick_labels[trace_i], index=ys.index)

                # Get the value by which the box will be coloured
                colorable_value = data_set.comparisons.loc[comp, color_by]

                # key-word args for the box
                boxkw = dict(
                    x=boxlabels, y=ys, name=colorable_value, boxpoints='outliers',
                    legendgroup=colorable_value,
                    customdata=fdr,
                    text=ys.index,
                    line=dict(color=box_colour_maps[color_by][colorable_value]),

                    hovertemplate=(
                            "<b>%{text}</b><br>" +
                            "Score: %{y:.2f}<br>" +
                            "FDR:   %{customdata:.2f}"
                    )
                )
                # include each treatment/whatever in the legend only once.
                if colorable_value in included:
                    boxkw['showlegend'] = False
                included.add(colorable_value)

                fig.add_trace(
                    go.Box(**boxkw)
                )

            # labels
            fig.update_layout(xaxis_title='Boxplot number',
                              yaxis_title=data_set.score_labels[score_type], )
            sort_by_opts = get_lab_val(order_by_categories + selected_genes)
            return fig, sort_by_opts, list(ordered_comps)

        return graph


    def spawn_datatable():
        table = Div(
            children=[
                create_datatable(columns_if_no_df=data_set.comparisons.columns)
            ],
            id=f'{PAGE_ID}-table-div',
            className="u-full-width",
            style={'margin-bottom': '30px'}
        )

        @app.callback(
            Output('msgv-table-div', 'children'),
            Input('msgv-stat-source-selector', 'value'),
            Input('ordered-comp-store', 'data'),
            State('msgv-gene-selector', 'value'),
            State('msgv-tabs', 'value'),
        )
        def update_datatable(score_type, ordered_comps, selected_genes, selected_tab):

            LOG.debug(f"{getfuncstr()}:\n\tordered_comps={ordered_comps}")

            data_tabs = data_set.get_score_fdr(score_type, score_type, )
            filtered_scores = data_tabs['score'].loc[selected_genes, ordered_comps]

            # DataTable structure:
            # First column trace numbers,
            #   then score/FDR for each gene selected,
            #   then the metadata.

            # Create the stat columns

            selected_fdr = data_tabs['fdr'].loc[filtered_scores.index, ordered_comps]
            filtered_scores = filtered_scores.reindex(columns=ordered_comps)
            filtered_scores.index = filtered_scores.index.map(lambda x: x + ' (Effect size)')
            selected_fdr.index = selected_fdr.index.map(lambda x: x + ' (FDR)')
            selected_stats = pd.concat([filtered_scores, selected_fdr], sort=False).T

            selected_stats = selected_stats.applymap(lambda n: f"{n:.3}")
            selected_stats.insert(0, 'Boxplot number', list(range(len(ordered_comps))))
            cols_oi = get_metadata_table_columns(public, PAGE_ID)['comp']

            selected_metadata = data_set.comparisons.loc[selected_stats.index, cols_oi]

            data_table_data = pd.concat([selected_stats, selected_metadata], axis=1)
            LOG.debug(f"{getfuncstr()} datatable data head:\n{str(data_table_data.head())}")
            return create_datatable(data_table_data)

        return table

    table = spawn_datatable()


    # select genes by name, and comparisons FDR in selected samples
    gene_selector = Div(
        children=[
            html.P('Select genes:', style=big_text_style),
            dcc.Dropdown(id='msgv-gene-selector', placeholder='Select genes', value=[],
                         options=get_gene_dropdown_lab_val(data_set, data_set.genes),
                         multi=True),
        ],
        style={'margin-bottom': '15px'})

    filter_cols = get_selector_table_filter_keys(public)['comp']

    filter_dropdowns = spawn_filter_dropdowns(
        PAGE_ID,
        'comp',
        filter_cols,
        data_set.comparisons,
        values={'Timepoint':['Matched time points']},
        card_header='Filter samples by categories below:'
    )

    text_header = Div([
        html.H1("Multi-Screen Gene Viewer"),
        html.P("Select your gene(s) of interest. Comparisons that show significant results "
               "(below adjustable FDR) for at least one selected gene are shown below. "
               "Box plots give the overall distribution of scores, and markers show specific genes. "
               "Diamonds indicate significant genes, and squares non-significant genes."),
    ])


    ### CONTROL PANEL for the plot

    stat_source_selectr = get_stat_source_selector('msgv', 'Analysis:', 'Significance source:')

    fdr_selectr = dbc.Card([
        dbc.CardHeader('Maximum FDR:'),
        #html.Label('FDR:', htmlFor='msgv-fdr-threshold'),
        dbc.CardBody([
            dcc.Input(id='msgv-fdr-threshold', type='number', min=0, max=2, step=0.01, value=0.1),
        ])
    ], style={'max-width': '170px'})

    # this is also used for output of one function, so is defined once here
    order_by_categories = ['Mean score', 'Treatment', 'Citation']
    colourable_categories = ['Treatment', 'Cell line', 'Citation', 'Library', 'KO']

    control_order_by = dbc.Card([
        dbc.CardHeader('Order by:'),
        dbc.CardBody([
            #html.Label('Order by:', htmlFor='msgv-order-by'),
            dcc.Dropdown(id='msgv-order-by', value=order_by_categories[0],
                         options=get_lab_val(order_by_categories),style={'width':'150px'}),
        ]),
    ], style={'max-width': '200px'})

    control_colour_by = dbc.Card([
        dbc.CardHeader('Colour by:'),
        dbc.CardBody([
            dcc.Dropdown(
                id='msgv-color-by', value='Treatment',
                options=get_lab_val(colourable_categories),
             style = {'width': '150px'}
            ),
        ]),
    ], style={'max-width': '200px'})


    control_panel = dbc.CardGroup(
        [
            fdr_selectr,
            control_order_by,
            control_colour_by,
            stat_source_selectr
        ],
        style={'max-width':'1100px'}
    )

    def spawn_comp_store() -> dcc.Store:
        """Spawn a dcc.Store that acts as a trigger for all output tabs.:
            id='comp-store'
            data={'selected_comps':List[str]}

        Register callback for storing mask of filtered comparisons.
        """

        store = dcc.Store(
            'comp-store',
            data={'selected_comps':[]}
        )
        @app.callback(
            Output('comp-store', 'data'),
            [Input('msgv-stat-source-selector', 'value'),
             Input('msgv-gene-selector', 'value'),
             Input('msgv-fdr-threshold', 'value'), ]
            + [Input(f'msgv-comp-filter-{cid}', 'value') for cid in filter_cols],

        )
        def filter_store_comps(
                score_type,
                selected_genes,
                fdr_thresh,
                *filters
        ):
            """Identify comparisons retained by filters, and containing significant
            genes."""
            if not selected_genes:
                raise PreventUpdate

            LOG.debug(f'{inspect.currentframe().f_code.co_name}\n'
                      f"\tgenes={selected_genes}, fdr_thresh={fdr_thresh}\n",
                      f"\tfilters={filters}")

            data_tabs = data_set.get_score_fdr(score_type, score_type, )

            filter_mask = pd.Series(True, index=data_set.comparisons.index)
            for filter_id, values in zip(filter_cols, filters):
                if values:
                    filter_mask = filter_mask & data_set.comparisons[filter_id].isin(values)

            fdr_mask = (data_tabs['fdr'].loc[selected_genes] <= fdr_thresh).any()
            comp_mask = (fdr_mask & filter_mask)
            comps = comp_mask[comp_mask].index
            return {'selected_comps': comps}

        return store

    # User selecting tabs changes tabs.value to individual tab value
    # NOTE: changing any Tab(value="") will require updating callbacks
    tabs = dcc.Tabs(
        id=f"msgv-tabs",
        value=f"msgv-boxplot-tab",
        className='selector-results-tabs',
        children=[
            dcc.Tab(
                label='Boxplots',
                value=f'msgv-boxplot-tab',
                className='data-tab', selected_className='data-tab--selected',
                children=[
                    spawn_boxplot_graph()
                ]
            ),

        ]
    )


    ### FINAL LAYOUT
    msgv_layout = Div([
        text_header,
        gene_selector,
        Div(filter_dropdowns, ),
        control_panel,
        tabs,
        table,
        spawn_comp_store(),
        dcc.Store('ordered-comp-store', data=[]) #used by datatable to order rows
    ], style={'display':'inline-block'})



    # get color map, asssiging colors to the most common values first, so that
    #   common things have different colours.
    def get_colour_map(list_of_things):
        return {thing:colours[i%len(colours)] for i, thing in enumerate(list_of_things)}
    box_colour_maps = {}
    for color_by in colourable_categories:
        ordered_things = data_set.comparisons.loc[:, color_by].value_counts().index
        cm = get_colour_map(ordered_things)
        box_colour_maps[color_by] = cm

    #todo create Tabs, id=f"{PAGE_ID}_tabs"
    #callbacks:

    # Callback to update graph
    # @app.callback(
    #     [Output('msgv-gene-boxplots', 'figure'),
    #      Output('msgv-table-div', 'children'),
    #      Output('msgv-order-by', 'options'),],
    #
    #     [Input('msgv-stat-source-selector', 'value'),
    #
    #      Input('msgv-gene-selector', 'value'),
    #      #Input('msgv-show-boxplots', 'value'),
    #      Input('msgv-fdr-threshold', 'value'),
    #
    #      Input('msgv-order-by', 'value'),
    #      Input('msgv-color-by', 'value'),
    #      ] + [Input(f'{PAGE_ID}-comp-filter-{cid}', 'value') for cid in filter_cols]
    # )
    # def update_figure(score_type,
    #                   selected_genes,
    #                   fdr_thresh,
    #                   order_by, color_by, *filters):
    #
    #     if not selected_genes:
    #         raise PreventUpdate
    #
    #     data_tabs = data_set.get_score_fdr(score_type, score_type, )
    #
    #     LOG.debug(f"MSGV: update_figure({score_type}, {selected_genes}, "
    #               f"fdr_thres={fdr_thresh}, ...")
    #
    #     # get boolean masks for which comparisons to include in the charts
    #     # Filter comparisons by metadata filters
    #     comparison_mask = pd.Series(True, index=data_set.comparisons.index)
    #     LOG.debug(f'update_figure: filters={filters} filter_cols={filter_cols}')
    #     for filter_id, values in zip(filter_cols, filters):
    #         if values:
    #             comparison_mask = comparison_mask & data_set.comparisons[filter_id].isin(values)
    #
    #     # Filter comparisons by FDR threshold
    #     fdr_mask = (data_tabs['fdr'].loc[selected_genes] <= fdr_thresh).any()
    #     filtered_scores = data_tabs['score'].loc[selected_genes, (fdr_mask & comparison_mask)]
    #
    #     # *assemble the figure*
    #     fig = go.Figure()
    #
    #     if (sum((fdr_mask & comparison_mask)) == 0) and selected_genes:
    #         fig.add_annotation(x=4, y=4,
    #                            text=f"None of the selected genes have FDR <= {fdr_thresh}",
    #                            showarrow=False,)
    #
    #     # determine order of plots
    #     if order_by in selected_genes:
    #         ordered_comps = filtered_scores.loc[order_by].sort_values().index.values
    #     elif order_by == 'Mean score':
    #         ordered_comps = filtered_scores.mean().sort_values().index.values
    #     elif order_by in order_by_categories[1:]:
    #         # subset the metadata to included comps and then sort by the order_by
    #         ordered_comps = data_set.comparisons.loc[filtered_scores.columns, order_by].sort_values().index
    #     else: # this shouldn't happen, though maybe I should have a "whatever" order option
    #         ordered_comps = filtered_scores.columns.values
    #
    #     # *plot the gene scatters*
    #
    #     # x tick labels
    #     trace_numbers = pd.Series([str(i) for i in range(1, len(ordered_comps)+1)])
    #
    #     if len(ordered_comps):
    #         citations = data_set.comparisons.loc[ordered_comps, 'Citation']
    #         treats = data_set.comparisons.loc[ordered_comps, 'Treatment'].fillna('')
    #         x_tick_labels = trace_numbers.values + '. ' +citations.values + ', ' + treats.values
    #     else:
    #         x_tick_labels = trace_numbers.values
    #
    #     #x_tick_labels = pd.Series(x_tick_labels, index=ordered_comps)
    #
    #     for gn in selected_genes:
    #         fdrs = data_tabs['fdr'].loc[gn, ordered_comps]
    #         mrkrs:pd.Series = fdrs <= fdr_thresh
    #         mrkrs[mrkrs == True] = 'diamond'
    #         mrkrs[mrkrs == False] = 'square'
    #
    #         fig.add_trace(
    #             go.Scatter(
    #                 x=x_tick_labels,
    #                 y=filtered_scores.loc[gn, ordered_comps],
    #                 mode='markers', name=gn,
    #                 marker_symbol=mrkrs.values,
    #                 marker={'size': 15, 'line':{'width':2, 'color':'DarkSlateGrey'}},
    #                 customdata=fdrs.apply(lambda n: f'{float(f"{n:.3g}"):g}'),
    #                 hovertemplate=f"{gn}"+"<br>Score: %{y}<br>FDR: %{customdata}<extra></extra>"
    #             ),
    #
    #         )
    #
    #     # add the boxplot traces
    #     included = set()
    #     # add a boxplot trace for each comparison
    #     for trace_i, comp in enumerate(ordered_comps):
    #         # these values define the boxplot
    #         ys = data_tabs['score'][comp]
    #         fdr = data_tabs['fdr'][comp]
    #         # This gives the X value for each y value used to create the boxplot, which
    #         #   is apparently required? I guess this is approximating Tidy formated data?
    #         boxlabels = pd.Series(x_tick_labels[trace_i], index=ys.index)
    #
    #         # Get the value by which the box will be coloured
    #         colorable_value = data_set.comparisons.loc[comp, color_by]
    #
    #         # key-word args for the box
    #         boxkw = dict(
    #             x=boxlabels, y=ys, name=colorable_value, boxpoints='outliers',
    #             legendgroup=colorable_value,
    #             customdata=fdr,
    #             text=ys.index,
    #             line=dict(color=box_colour_maps[color_by][colorable_value]),
    #
    #             hovertemplate = (
    #                     "<b>%{text}</b><br>" +
    #                     "Score: %{y:.2f}<br>" +
    #                     "FDR:   %{customdata:.2f}"
    #             )
    #         )
    #         # include each treatment/whatever in the legend only once.
    #         if colorable_value in included:
    #             boxkw['showlegend'] = False
    #         included.add(colorable_value)
    #
    #         fig.add_trace(
    #             go.Box(**boxkw)
    #         )
    #
    #     # labels
    #     fig.update_layout(xaxis_title='Boxplot number',
    #                       yaxis_title=data_set.score_labels[score_type],)
    #
    #     # DataTable structure:
    #     # First column trace numbers,
    #     #   then score/FDR for each gene selected,
    #     #   then the metadata.
    #
    #     # Create the stat columns
    #     selected_fdr = data_tabs['fdr'].loc[filtered_scores.index, ordered_comps]
    #     filtered_scores = filtered_scores.reindex(columns=ordered_comps)
    #     filtered_scores.index = filtered_scores.index.map(lambda x: x+' (Effect size)')
    #     selected_fdr.index = selected_fdr.index.map(lambda x: x+' (FDR)')
    #     selected_stats = pd.concat([filtered_scores, selected_fdr], sort=False).T
    #
    #     selected_stats = selected_stats.applymap(lambda n: f"{n:.3}")
    #     selected_stats.insert(0, 'Boxplot number', trace_numbers)
    #     cols_oi = get_metadata_table_columns(public, PAGE_ID)['comp']
    #
    #     selected_metadata = data_set.comparisons.loc[selected_stats.index, cols_oi]
    #
    #     data_table_data = pd.concat([selected_stats, selected_metadata], axis=1)
    #
    #     dtable = create_datatable(data_table_data)
    #
    #     sort_by_opts = get_lab_val(order_by_categories+ selected_genes)
    #
    #     return fig, dtable, sort_by_opts

    return msgv_layout



if __name__ == '__main__':
    from crispr_screen_viewer.functions_etc import get_cmdline_options
    from crispr_screen_viewer.functions_etc import launch_page
    launch_page(*get_cmdline_options(), name='MSGV', initiate=initiate)


