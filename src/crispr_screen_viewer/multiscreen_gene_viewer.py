#!/usr/bin/env python

import pandas as pd
from dash import dash, dcc, html, Input, Output, State
import dash_bio
from dash.exceptions import PreventUpdate
Div = html.Div
import plotly.graph_objs as go

from loguru import logger

from crispr_screen_viewer.shared_components import (
    create_datatable,
    get_lab_val,
    get_stat_source_selector,
    colours,
    big_text_style,
    logger,
    spawn_filter_dropdowns,
)

from crispr_screen_viewer.functions_etc import (
    get_metadata_table_columns,
    get_selector_table_filter_keys,
    getfuncstr,
    index_of_true,
    bicolour_cmap,

)

from crispr_screen_viewer.dataset import DataSet


border_style = {'border': '4px solid #3DD178',
                'border-radius': '14px'}

import dash_bootstrap_components as dbc
import plotly.express as px

PAGE_ID = 'msgv'


# So the structure here is selected data (comps and genes, tables?) are selected and
#   put into a store, this store is the trigger for updating the tabs (tabs don't also
#   don't update unless selected).


def initiate(app, data_set:DataSet, public=True) -> Div:
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
            style={'height':'1000px', 'width':'1500px'}
        )

        @app.callback(
            Output('msgv-gene-boxplots', 'figure'),
            Output('msgv-order-by', 'options'),
            Output('ordered-comp-store', 'data'),

            Input('msgv-fdr-threshold', 'value'),
            Input('msgv-show-outliers-radio', 'value'),
            Input('msgv-max-boxplots', 'value'),
            Input('msgv-stat-source-selector', 'value'),
            Input('comp-store', 'data'),
            Input('msgv-tabs', 'value'),
            Input('msgv-order-by', 'value'),
            Input('msgv-color-by', 'value'),
            Input('msgv-update-genes-button', 'n_clicks'),
            State('msgv-gene-selector', 'value'),
        )
        def update_boxplot_figure(
                fdr_threshold,
                show_outliers_option,
                max_boxplots,
                score_type,
                filtered_comps:str,
                selected_tab,
                order_by,
                color_by,
                update_genes,
                selected_genes,
        ):
            if (selected_tab != 'msgv-boxplot-tab') or (not selected_genes):
                logger.debug("Not updating")
                raise PreventUpdate



            logger.debug(
                f"{show_outliers_option=}, {max_boxplots=}, {score_type=}, {filtered_comps=}, {selected_tab=}, {order_by=}, {color_by=}, { selected_genes=}, {fdr_threshold=}"
            )


            # *assemble the figure*
            fig = go.Figure()

            hit_table = data_set.get_score_fdr(
                score_type,
                comparisons=filtered_comps,
                genes=selected_genes,
                fdr_max=fdr_threshold,
            )

            hit_table_scores = hit_table['score']

            comps_with_hits = hit_table_scores.columns

            all_genes_score_fdr = data_set.get_score_fdr(
                score_type,
                comparisons=comps_with_hits,
            )

            # todo make no-hits-message always in the middle of the figure
            if (len(comps_with_hits) == 0) and selected_genes:
                fig.add_annotation(x=1.5, y=2.5,
                                   text=f"None of the selected genes have FDR <= {fdr_threshold}",
                                   showarrow=False, )

            show_outliers = True
            if show_outliers_option == 'never':
                show_outliers = False
            elif show_outliers_option == 'sometimes':
                if len(comps_with_hits) > max_boxplots:
                    show_outliers = False



            # data_tabs:Dict[str, pd.DataFrame] = data_set.get_score_fdr(
            #     score_type,
            #     comparisons=comps,
            # )
            # filtered_scores = data_tabs['score'].reindex(index=selected_genes)


            # determine order of plots
            if order_by in selected_genes:
                ordered_comps = hit_table_scores.loc[order_by].sort_values().index.values
            elif order_by == 'Mean score':
                ordered_comps = hit_table_scores.mean().sort_values().index.values
            elif order_by in order_by_categories[1:]:
                # subset the metadata to included comps and then sort by the order_by
                ordered_comps = data_set.comparisons.loc[hit_table_scores.columns, order_by].sort_values().index
            else:  # this shouldn't happen, though maybe I should have a "whatever" order option
                ordered_comps = hit_table_scores.columns.values

            # x tick labels
            if len(ordered_comps):
                x_tick_labels = get_numbered_tick_labels(ordered_comps)
            else:
                x_tick_labels = get_numbered_tick_labels(ordered_comps, just_numbers=True)

            # x_tick_labels = pd.Series(x_tick_labels, index=ordered_comps)
            # note: looping and adding traces was required, using the plotly express functions
            #   resulted in non-overlapping X positions for the gene scatters and boxplots.

            for gn in selected_genes:
                if gn not in all_genes_score_fdr['fdr'].index:
                    logger.debug(f"Skipping gene {gn}")
                    continue
                fdrs:pd.Series = all_genes_score_fdr['fdr'].loc[gn, ordered_comps]
                mrkrs = ['diamond' if (f<fdr_threshold) else 'square' for f in fdrs]


                fig.add_trace(
                    go.Scatter(
                        x=x_tick_labels,
                        y=all_genes_score_fdr['score'].loc[gn, ordered_comps],
                        mode='markers',
                        name=gn,
                        marker_symbol=mrkrs,
                        marker={'size': 15, 'line': {'width': 2, 'color': 'DarkSlateGrey'}},
                        customdata=fdrs.apply(lambda n: f'{float(f"{n:.3g}"):g}' if not pd.isna(n) else ''),
                        hovertemplate=f"{gn}" + "<br>Score: %{y}<br>FDR: %{customdata}<extra></extra>"
                    ),

                )

            # add the boxplot traces
            included = set()
            # Add a boxplot trace for each comparison
            for trace_i, comp in enumerate(ordered_comps):
                # these values define the boxplot
                ys = all_genes_score_fdr['score'][comp]
                fdr = all_genes_score_fdr['fdr'][comp]


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

                if not show_outliers:
                    boxkw['boxpoints'] = False

                # include each treatment/whatever in the legend only once.
                if colorable_value in included:
                    boxkw['showlegend'] = False
                included.add(colorable_value)

                fig.add_trace(
                    go.Box(**boxkw)
                )

            fig.update_xaxes(tickangle=45)

            # labels
            from crispr_screen_viewer.dataset import ANALYSESTYPES
            fig.update_layout(xaxis_title='Plot number',
                              yaxis_title=ANALYSESTYPES[score_type].score_label, )
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

            Input('msgv-fdr-threshold', 'value'),
            Input('msgv-stat-source-selector', 'value'),
            Input('comp-store', 'data'),
            Input('msgv-tabs', 'value'),
            Input('msgv-cluster-missing-method', 'value'),
            Input('cluster-height', 'value'),
            Input('cluster-width', 'value'),
            Input('msgv-update-genes-button', 'n_clicks'),
            State('msgv-gene-selector', 'value'),

        )
        def update_clustergrams(
                fdr_threshold:float,
                score_type:str,
                comps:list[str],
                selected_tab:str,
                missing_method,
                height, width,
                update_genes_button,
                selected_genes,
        ):

            data_missing = False
            heatmap_only = False

            if selected_tab != 'msgv-clustergram-tab':
                logger.debug('not updating.')
                raise PreventUpdate

            if not selected_genes:
                return [html.P('Select genes above.', className='missing-data-paragraph')]

            logger.debug(f"{getfuncstr()} updating.")

            score_fdr = data_set.get_score_fdr(
                score_type,
                genes=selected_genes,
                comparisons=comps,
            )

            # drop genes/comps where there's no significant reading
            sig:pd.DataFrame = score_fdr['fdr'] <= fdr_threshold
            filtered_scores = score_fdr['score'].loc[
                sig.any(axis=1),
                sig.any(axis=0)
            ]

            comp_label_dict = {c:l for c, l in
                zip(
                    filtered_scores.columns,
                    get_numbered_tick_labels(filtered_scores.columns)
                )
            }
            logger.debug(f'Pre-filter for NaN:\n{filtered_scores}')
            # filter data if there are comps that don't have scores for genes
            if filtered_scores.isna().any().any():
                data_missing = True
                missing_comps = index_of_true(filtered_scores.isna().any())
                missing_genes = index_of_true(filtered_scores.isna().any(axis=1))
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

            logger.debug(f'Post-filter for NaN:\n{filtered_scores}')

            # deal with cases where we can't display any graphs.
            if not selected_genes or (len(selected_genes) < 2):
                return [html.P('Select at least 2 genes', className='no-data-paragraph')]
            elif (len(comps) == 0) and (selected_genes > 2):
                if not data_missing:
                    return [html.P(f"None of the selected genes have FDR <= {fdr_threshold}",
                                   className='no-data-paragraph')]
                else:
                    return [html.P(f"None of the selected genes have FDR <= {fdr_threshold}\n"
                                   f"Some samples don't have scores for some genes.",
                                   className='no-data-paragraph')]
            elif (0 in filtered_scores.shape) or (1 in filtered_scores.shape):
                return [html.P(
                    "Not enough observations left after removing missing genes/comparisons",
                    className='no-data-paragraph'
                )]

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

    # data_table_columns = dict(
    #     treatment_label='Treatment',
    #     dose='Dose',
    #     timepoint='Timepoint',
    #     days_grown='Days grown',
    #     cell_line='Cell linen',
    #     ko='KO',
    #     library='Library',
    #     citation='Citation'
    # )

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
            # DataTable structure:
            # First column trace numbers,
            #   then score/FDR for each gene selected,
            #   then the metadata.

            logger.debug(f"{score_type=}, {ordered_comps=}, {selected_genes=}, {selected_tab=}")

            data_tabs = data_set.get_score_fdr(
                score_type,
                score_type,
                genes=selected_genes,
                comparisons=ordered_comps
            )

            filtered_scores = data_tabs['score'].reindex(columns=ordered_comps)
            filtered_fdr = data_tabs['fdr'].reindex(columns=ordered_comps)

            # Create the stat columns
            filtered_scores.index = filtered_scores.index.map(lambda x: x + ' (score)')
            filtered_fdr.index = filtered_fdr.index.map(lambda x: x + ' (FDR)')
            selected_stats = pd.concat([filtered_scores, filtered_fdr], sort=False).T

            selected_stats = selected_stats.map(lambda n: f"{n:.3}")
            selected_stats.insert(0, 'Boxplot number',
                                  list(range(1, len(ordered_comps)+1)))
            cols_oi = get_metadata_table_columns(public, PAGE_ID)['comp']

            selected_metadata = data_set.comparisons.loc[selected_stats.index, cols_oi]

            data_table_data = pd.concat([selected_stats, selected_metadata], axis=1)
            logger.debug(f"{getfuncstr()} datatable data head:\n{str(data_table_data.head())}")

            return create_datatable(data_table_data)

        return table


    # select genes by name, and comparisons FDR in selected samples
    gene_selector = Div(
        children=[
            html.P('Select genes:', style=big_text_style),

            Div([
                dcc.Dropdown(
                    id='msgv-gene-selector', placeholder='Select genes',
                     value=[],
                     options=data_set.dropdown_gene_labels(data_set.genes),
                     multi=True
                ),
                html.Button(
                    'Submit',
                    id='msgv-update-genes-button',
                    n_clicks=0,
                ),
            ])
        ],
        style={'margin-bottom': '15px'}
    )

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
        html.H1("Explore gene results"),
        html.P("Select your gene(s) of interest in the box below. Comparisons that show significant results "
               "(below adjustable FDR) for at least one selected gene are shown below. "
               "Box plots give the overall distribution of scores, and markers show specific genes. "
               "Diamonds indicate significant genes, and squares non-significant genes."),
    ])


    ### CONTROL PANEL for the plot

    def spawn_control_cardgroup() -> dbc.CardGroup:
        stat_source_selectr = get_stat_source_selector('msgv', 'Analysis method:')

        # doing the full boxplots can be slow
        default_max_boxplots = 12
        show_outliers_ctrl = dbc.Card(
            id='msgv-show-boxplots-card',
            children=[
                dbc.CardHeader('Show outliers:'),
                dbc.CardBody([
                    dcc.RadioItems(
                        id='msgv-show-outliers-radio',
                        options=[
                            {'label':'Always (may be very slow)', 'value':'always'},
                            {'label':'Never', 'value':'never'},
                            {'label': 'When boxplots fewer than', 'value': 'sometimes'},
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
                    dcc.Dropdown(
                        id='msgv-order-by',
                        clearable=False,
                        value=order_by_categories[0],
                        options=get_lab_val(order_by_categories),
                        style={'width':'150px'}),
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
                    clearable=False,
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
                    fdr_selectr,
                    stat_source_selectr,
                    control_order_by,
                    control_colour_by,
                    clustergram_missing,
                    clustergram_controls,
                    show_outliers_ctrl,
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

            logger.warn(f'MSGV: Unknown tab value {selected_tab}')
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
            [Input(f'msgv-comp-filter-{cid}', 'value') for cid in filter_cols],
        )
        def filter_store_comps(
                *filters
        ):
            """Filters comparisons by values in the dropdowns, i.e.
            ['Treatment', 'Cell line', 'KO', 'Timepoint',
                'Library', 'Citation',]"""


            # apply filters from the app
            filter_mask = pd.Series(True, index=data_set.comparisons.index)
            for filter_id, values in zip(filter_cols, filters):
                if values:
                    filter_mask = filter_mask & data_set.comparisons[filter_id].isin(values)

            included_comparisons = filter_mask.index[filter_mask]


            return included_comparisons

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


