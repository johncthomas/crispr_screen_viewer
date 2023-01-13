#!/usr/bin/env python

import itertools
import inspect
import typing

import pandas as pd
from dash import dash, dcc, html, Input, Output, State
import dash_bio
from dash.exceptions import PreventUpdate
Div = html.Div
import plotly.graph_objs as go

import pathlib, os, logging

from typing import Collection,  Union, Dict

from crispr_screen_viewer.shared_components import (
    create_datatable,
    get_lab_val,
    get_gene_dropdown_lab_val,
    #get_reg_stat_selectors,
    get_stat_source_selector,
    colours,
    big_text_style,
    LOG,
    spawn_filter_dropdowns,
)

from crispr_screen_viewer.functions_etc import (
    get_metadata_table_columns,
    get_selector_table_filter_keys,
    getfuncstr,
    index_of_true,
    bicolour_cmap,
)


border_style = {'border': '4px solid #3DD178',
                'border-radius': '14px'}

import dash_bootstrap_components as dbc
import plotly.express as px

PAGE_ID = 'msgv'


# Maybe clustergram it's own page,  using the selector-tables to select comps, and gene box

# todo Make MSGV work like the others with the comparison selection tables
#   with select by FDR as an additional option.
# Obviously this is a bigger deal and further down the pipe

# Crucially, lots of callbacks sharing inputs is *fine*. You just can't have multiple
#   callbacks that share outputs.
# So the structure here is selected data (comps and genes, tables?) are selected and
#   put into a store, this store is the trigger for updating the tabs (tabs don't also
#   don't update unless selected).



def initiate(app, data_set, public=True) -> Div:
    """Register callbacks to app, generate layout"""

    order_by_categories = ['Mean score', 'Treatment', 'Citation']
    colourable_categories = ['Treatment', 'Cell line', 'Citation', 'Library', 'KO']

    def get_numbered_tick_labels(comps, just_numbers=False):
        """list of strings giving number, citation and treatment info for
        comps. If just_numbers, return only list of strings of numbers."""

        trace_numbers = pd.Series([str(i) for i in range(1, len(comps) + 1)])
        if just_numbers:
            return trace_numbers

        citations = data_set.comparisons.loc[comps, 'Citation']
        treats = data_set.comparisons.loc[comps, 'Treatment'].fillna('')
        return trace_numbers.values + '. ' + citations.values + ', ' + treats.values

    def spawn_boxplot_graph() -> dcc.Graph:

        graph = dcc.Graph(
            id='msgv-gene-boxplots',
            figure=go.Figure(),
            style={'height':'800px', 'width':'1500px'}
        )

        @app.callback(
            Output('msgv-gene-boxplots', 'figure'),
            Output('msgv-order-by', 'options'),
            Output('ordered-comp-store', 'data'),

            Input('msgv-show-boxplots-radio', 'value'),
            Input('msgv-max-boxplots', 'value'),
            Input('msgv-stat-source-selector', 'value'),
            Input('comp-store', 'data'),
            Input('msgv-tabs', 'value'),
            Input('msgv-order-by', 'value'),
            Input('msgv-color-by', 'value'),

            State('msgv-gene-selector', 'value'),
            State('msgv-fdr-threshold', 'value'),
        )
        def update_boxplot_figure(
                show_boxplots_option,
                max_boxplots,
                score_type,
                comps_dict:Dict[str, typing.List[str]],
                selected_tab,
                order_by,
                color_by,

                selected_genes,
                fdr_thresh
        ):
            comps = comps_dict['selected_comps']
            LOG.debug(f"{getfuncstr()}: {comps}")

            if selected_tab != 'msgv-boxplot-tab':
                raise PreventUpdate

            show_boxplots = True
            if show_boxplots_option == 'never':
                show_boxplots = False
            elif show_boxplots_option == 'sometimes':
                if len(comps) > max_boxplots:
                    show_boxplots = False


            data_tabs:Dict[str, pd.DataFrame] = data_set.get_score_fdr(score_type, score_type, )
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
            if len(ordered_comps):
                x_tick_labels = get_numbered_tick_labels(ordered_comps)
            else:
                x_tick_labels = get_numbered_tick_labels(ordered_comps, just_numbers=True)

            # x_tick_labels = pd.Series(x_tick_labels, index=ordered_comps)

            for gn in selected_genes:
                fdrs:pd.Series = data_tabs['fdr'].loc[gn, ordered_comps]
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
            if show_boxplots:
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
            fig.update_layout(xaxis_title='Plot number',
                              yaxis_title=data_set.score_labels[score_type], )
            sort_by_opts = get_lab_val(order_by_categories + selected_genes)
            return fig, sort_by_opts, list(ordered_comps)

        return graph


    def spawn_clustergram() -> Div:
        clustergram_id = 'msgv-clustergram-div'
        clustergram_div = Div(
            id=clustergram_id,
            children=[html.P('')],
            style={'display':'flex'},
        )

        @app.callback(
            Output(clustergram_id, 'children'),

            Input('msgv-stat-source-selector', 'value'),
            Input('comp-store', 'data'),
            Input('msgv-tabs', 'value'),
            Input('msgv-cluster-missing-method', 'value'),
            Input('cluster-height', 'value'),
            Input('cluster-width', 'value'),

            State('msgv-gene-selector', 'value'),
            State('msgv-fdr-threshold', 'value'),
        )
        def update_clustergrams(
                score_type,
                comps_dict: Dict[str, typing.List[str]],
                selected_tab,
                missing_method,
                height, width,

                selected_genes,
                fdr_thresh,

        ):

            data_missing = False
            heatmap_only = False

            if selected_tab != 'msgv-clustergram-tab':
                raise PreventUpdate

            if not selected_genes:
                return [html.P('Select genes above.', className='missing-data-paragraph')]

            LOG.debug(f"{getfuncstr()} updating.")

            comps = comps_dict['selected_comps']
            LOG.debug(f"{getfuncstr()}: {comps}")
            data_tabs = data_set.get_score_fdr(score_type, score_type, )
            filtered_scores:pd.DataFrame = data_tabs['score'].loc[selected_genes, comps]

            comp_label_dict = {c:l for c, l in
                zip(
                    filtered_scores.columns,
                    list(get_numbered_tick_labels(filtered_scores.columns))
                )
            }

            # filter data if there are comps that don't have scores for genes
            if filtered_scores.isna().any().any():
                data_missing = True
                missing_comps = index_of_true(filtered_scores.isna().any())
                missing_genes = index_of_true(filtered_scores.isna().any(1))
                # to be returned in the Div
                misslabs = '\n    '.join([comp_label_dict[c] for c in missing_comps])
                missing_text = f"Sample(s) has missing genes:\n  {misslabs}"
                if missing_method == 'drop_comps':
                    filtered_scores = filtered_scores.drop(missing_comps, axis='columns')
                elif missing_method == 'drop_genes':
                    filtered_scores = filtered_scores.drop(missing_genes, axis='index')
                # Don't filter if heatmap_only
                elif missing_method == 'heatmap_only':
                    heatmap_only = True
                else:
                    raise RuntimeError(f'Unknown missing_method: {missing_method}.')

            # deal with cases where we can't display any graphs.
            if not selected_genes or (len(selected_genes) < 2):
                return [html.P('Select at least 2 genes', className='no-data-paragraph')]
            elif (len(comps) == 0) and (selected_genes > 2):
                if not data_missing:
                    return [html.P(f"None of the selected genes have FDR <= {fdr_thresh}",
                                   className='no-data-paragraph')]
                else:
                    return [html.P(f"None of the selected genes have FDR <= {fdr_thresh}\n"
                                   f"Some samples don't have scores for some genes.",
                                   className='no-data-paragraph')]

            mn, mx = filtered_scores.min().min(), filtered_scores.max().max()
            cmap_slice = bicolour_cmap(mn, mx, )

            if not heatmap_only:
                mainfig = dash_bio.Clustergram(
                    data=filtered_scores.values,
                    row_labels=list(filtered_scores.index),
                    column_labels=[comp_label_dict[c] for c in filtered_scores.columns],
                    # maximum linkage distance to join dendrogram clusters
                    color_threshold={
                        'row': 100,
                        'col': 100
                    },
                    height=height,
                    width=width,
                    color_map=cmap_slice,
                    center_values=False,
                    standardize=False,


                )
            else:
                tab = filtered_scores
                tab.columns = tab.columns.map(comp_label_dict)
                mainfig = px.imshow(
                    tab,
                    color_continuous_scale=cmap_slice,
                    height = height,
                    width = width,
                )



            output_div = [dcc.Graph(figure=mainfig)]
            if data_missing:
                output_div.append(
                    html.P(
                        missing_text,
                        style={'whiteSpace': 'pre-wrap', },
                        className='missing-data-paragraph'
                    )
                )

            return output_div

        return clustergram_div

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
            selected_stats.insert(0, 'Boxplot number',
                                  list(range(1, len(ordered_comps)+1)))
            cols_oi = get_metadata_table_columns(public, PAGE_ID)['comp']

            selected_metadata = data_set.comparisons.loc[selected_stats.index, cols_oi]

            data_table_data = pd.concat([selected_stats, selected_metadata], axis=1)
            LOG.debug(f"{getfuncstr()} datatable data head:\n{str(data_table_data.head())}")

            return create_datatable(data_table_data)

        return table


    # select genes by name, and comparisons FDR in selected samples
    gene_selector = Div(
        children=[
            html.P('Select genes:', style=big_text_style),

            dcc.Dropdown(id='msgv-gene-selector', placeholder='Select genes',
                         value=[],
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

    def spawn_control_cardgroup() -> dbc.CardGroup:
        stat_source_selectr = get_stat_source_selector('msgv', 'Analysis:', 'Analysis:')

        # doing the full boxplots can be slow
        default_max_boxplots = 20
        show_boxplots_ctrl = dbc.Card(
            id='msgv-show-boxplots-card',
            children=[
                dbc.CardHeader('Show boxplots:'),
                dbc.CardBody([
                    dcc.RadioItems(
                        id='msgv-show-boxplots-radio',
                        options=[
                            {'label':'Always (may be very slow)', 'value':'always'},
                            {'label':'Never', 'value':'never'},
                            {'label': 'When no more than', 'value': 'sometimes'},
                        ],
                        value='sometimes',
                    ),
                    dcc.Input(
                        id='msgv-max-boxplots',
                        type='number',
                        value=default_max_boxplots,
                    )
                ])
            ]
        )


        fdr_selectr = dbc.Card([
            dbc.CardHeader('Maximum FDR:'),
            #html.Label('FDR:', htmlFor='msgv-fdr-threshold'),
            dbc.CardBody([
                dcc.Input(id='msgv-fdr-threshold', type='number', min=0, max=2, step=0.01, value=0.1),
            ])
        ], style={'width': '170px'})

        control_order_by = dbc.Card(
            id = 'msgv-order-by-card',
            children=[
                dbc.CardHeader('Order by:'),
                dbc.CardBody([
                    #html.Label('Order by:', htmlFor='msgv-order-by'),
                    dcc.Dropdown(id='msgv-order-by', value=order_by_categories[0],
                             options=get_lab_val(order_by_categories),style={'width':'150px'}),
                ]),
            ],
            style={'width': '250px'})
            #class_name="card col-md-2")

        control_colour_by = dbc.Card(
            id='msgv-color-by-card',
            children=[dbc.CardHeader('Colour by:'),
            dbc.CardBody([
                dcc.Dropdown(
                    id='msgv-color-by', value='Treatment',
                    options=get_lab_val(colourable_categories),
                 style = {'width': '150px'}
                ),
            ]),
        ], style={'width': '250px'})


        missing_opts = [
            ('Drop samples','drop_comps'),
            ('Drop genes','drop_genes',),
            ("Don't cluster",'heatmap_only'),
        ]
        missing_opts = [{'label':l, 'value':v} for l, v in missing_opts]
        clustergram_missing = dbc.Card(
            id='msgv-cluster-missing-method-card',
            children=[dbc.CardHeader('When genes missing:'),
            dbc.CardBody([
                dcc.RadioItems(
                    id='msgv-cluster-missing-method',
                    value='drop_comps',
                    options=missing_opts,
                    style = {'width': '150px'},

                ),
            ]),
        ], style={'width': '170px'})

        clustergram_controls = dbc.Card(
            id='msgv-cluster-controls-card',
            children=[dbc.CardHeader('Figure size:'),
            dbc.CardBody([
                dbc.InputGroup([
                    dbc.InputGroupText('Height:'),
                    dbc.Input(id='cluster-height', type='number', value=600),
                ]),
                dbc.InputGroup([
                    dbc.InputGroupText('Width: '),
                    dbc.Input(id='cluster-width', type='number', value=1200),
                ]),
            ]),
        ], style={'width': '170px'})

        control_panel = dbc.CardGroup(
                [
                    show_boxplots_ctrl,
                    fdr_selectr,
                    stat_source_selectr,
                    control_order_by,
                    control_colour_by,
                    clustergram_missing,
                    clustergram_controls,
                ],
                #class_name='card-group row',
                style={'width':'fit-content'}
            )


        @app.callback(
            Output('msgv-show-boxplots-card', 'style'),
            Output('msgv-order-by-card', 'style'),
            Output('msgv-color-by-card', 'style'),
            Output('msgv-cluster-missing-method-card', 'style'),
            Output('msgv-cluster-controls-card', 'style'),


            Input('msgv-tabs', 'value')
        )
        def hide_irrelevant(selected_tab):
            show = {'width':'250px'}
            hide = {'display':'none'}
            if selected_tab == 'msgv-boxplot-tab':
                return (show, show, show, hide, hide)
            elif selected_tab == 'msgv-clustergram-tab':
                return (hide, hide, hide, show, show)

            LOG.warn(f'MSGV: Unknown tab value {selected_tab}')
            return (show, show, show, show)

        return control_panel


    def get_colour_map():
        # get color map, asssiging colors to the most common values first, so that
        #   common things have different colours.
        def get_colour_map(list_of_things):
            return {thing: colours[i % len(colours)] for i, thing in enumerate(list_of_things)}

        box_colour_maps = {}
        for color_by in colourable_categories:
            ordered_things = data_set.comparisons.loc[:, color_by].value_counts().index
            cm = get_colour_map(ordered_things)
            box_colour_maps[color_by] = cm
        return box_colour_maps

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
                      f"\tgenes={selected_genes}, fdr_thresh={fdr_thresh}\n"
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
    # NOTE: changing any Tab(value="") will require updating callbacks,
    #   check for Input('msgv-tabs', 'value')
    tabs = dcc.Tabs(
        id=f"msgv-tabs",
        value=f"msgv-boxplot-tab",
        className='selector-results-tabs',
        children=[
            dcc.Tab(
                label='Scores',
                value=f'msgv-boxplot-tab',
                className='data-tab', selected_className='data-tab--selected',
                children=[
                    spawn_boxplot_graph()
                ]
            ),
            dcc.Tab(
                label='Clustergram', value='msgv-clustergram-tab',
                className='data-tab', selected_className='data-tab--selected',
                children=[
                    spawn_clustergram()
                ]
            )

        ]
    )

    table = spawn_datatable()
    control_panel = spawn_control_cardgroup()
    box_colour_maps = get_colour_map()

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


    return msgv_layout



if __name__ == '__main__':
    from crispr_screen_viewer.functions_etc import get_cmdline_options
    from crispr_screen_viewer.functions_etc import launch_page
    launch_page(*get_cmdline_options(), name='MSGV', initiate=initiate)


