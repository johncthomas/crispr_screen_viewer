import logging

import pandas as pd
import numpy as np

from dash import dash, dcc, html, Input, Output, State, dash_table, callback_context

from dash.exceptions import PreventUpdate

import plotly.graph_objs as go

import pathlib, os
from dash.dependencies import Input, Output, State
from typing import Collection, Union, Dict, List
from crispr_screen_viewer.functions_etc import DataSet
from crispr_screen_viewer.shared_components import (
    get_lab_val,
    get_treatment_label,
    get_annotation_dicts,
    big_text_style,
    styles,
    LOG,
    get_stat_source_selector
)

from crispr_screen_viewer.functions_etc import (
    DataSet,
    parse_expid
)

Div = html.Div


cell_text_style = {
    'font-family':'Helvetica',
    'font-size':'15px',
    'text-align':'left',
    "whiteSpace": "pre-line",
}


cell_number_style = {
    'font-family':'monospace',
    'font-size':'17px'
}


def initiate(app, data_set:DataSet, public_version=False) -> Div:
    """Source directory should contain the relevant info: metadata.csv,
    screen_analyses and expyaml directories."""

    comparisons = data_set.comparisons
    try:
        comparisons = comparisons.drop('Available analyses', 1)
    except:
        pass

    if public_version:
        source_display = 'none'
    else:
        source_display = 'block'


    # ** helpful functions **
    def column_dict(c, markdown=('DOI',), ):
        """return {'name':k, 'id':k} for all k except "DOI" which
        sets style to markdown"""
        if c in markdown:
            return {"name": c, "id": c, 'type': 'text', 'presentation':'markdown'}
        else:
            return {'name':c, 'id':c}

    def doi_to_link(doi):
        return f"[{doi}](https://doi.org/{doi})"


    def get_gene_selector(idprefix) -> list:
        return [
            html.Label('Select genes:',
                       htmlFor=f'{idprefix}-gene-dropdown',
                       style=big_text_style),
            dcc.Dropdown(
                id=f'{idprefix}-gene-dropdown',
                placeholder='Select genes by name',
                multi=True,
                value=[],
                options=[]),
        ]


    # ################## #
    # ****COMPONENTS**** #
    # ################## #

    # **GRAPHS****
    volcano = dcc.Graph(
        id='volcano0',
        config={
            'editable':True,
            'edits':{'annotationPosition':False},
        },
        style={'height': '800px',
               'padding':'0px'},
        # # note, figure gets replaced in update_volcano callback so add
        # # any changes there.
        # figure={'layout':{'clickmode':'event+select',
        #                   'dragmode':'select'}, }
    )

    volcano_layout = [
        Div(get_stat_source_selector('volcano', 'Analysis:')),
        Div([volcano]),
        Div(get_gene_selector('volc'))
    ]

    # ***TABLES***
    # comptab and exptab are data selectors.
    # Experiments, select experiment, switch to comparisons table with comp.expid == selected_expid
    #   each needs own filter dropdowns
    # The Experiments and Comparisons tables identified using 'exp' and 'comp'

    # Timepoint, renaming to readable
    for val, lab in [('fromstart', 'From experiment start'),
                     ('otherprior', 'From midpoint'),
                     ('endpoints', 'Matched time points')]:
        m = comparisons.Timepoint.str.startswith(val)
        comparisons.loc[m, 'Time point group'] = lab



    comptab_data = comparisons.copy()
    LOG.debug('Comparison tab columns'+str(comptab_data.columns))

    # get the experiment data for these comparisons
    for k in ('DOI', 'Citation'):
        comptab_data.loc[:, k] = comptab_data['Experiment ID'].apply(
            lambda exp: data_set.experiments_metadata.loc[exp, k]
        )
    comptab_data['DOI'] = comptab_data['DOI'].apply(doi_to_link)
    comptab_data.loc[:, 'Comparison ID'] = comptab_data['Comparison ID']

    # # add authors to comparisons, will be used in filtering call back
    # comparisons.loc[:, 'Author'] = comptab_data.Citation.apply(lambda x: x.split(', ')[0])

    # ** table_of_experiments **
    # get the grouped 'Treatment', 'Cell line', 'KO', 'Library' for experiments
    groups_exp = comparisons.groupby('Experiment ID')

    # DF with each cell an ordered list of non-repeating strings
    # used for selecting rows and to construct displayed DataTable
    exptab_data = {expid:{} for expid in comparisons['Experiment ID'].unique()}
    # for filtering, will keep cell values as lists
    exptab_data_searchable = {expid:{} for expid in comparisons['Experiment ID'].unique()}

    # columns determined by what's selected above
    # first group values comparisons
    for col in comptab_data.columns:
        for expid, rows in groups_exp.groups.items():
            vals = sorted(set(comptab_data.loc[rows, col].dropna().astype(str).values))
            exptab_data[expid][col] = ', '.join(vals)
            exptab_data_searchable[expid][col] = vals

    exptab_data = pd.DataFrame(exptab_data).T
    exptab_data_searchable = pd.DataFrame(exptab_data_searchable).T
    exptab_data.index.name = 'Experiment ID'
    exptab_data_searchable.index.name = 'Experiment ID'

    table_dataframes = {'exp':exptab_data, 'comp':comptab_data}

    if public_version:
        tab_columns = {
            'exp':[ 'Citation', 'Treatment', 'Cell line', 'KO',  'Library', 'Experiment ID',],
            'comp':['Comparison ID',  'Treatment', 'Dose', 'Time point group',
                    'Growth inhibition %', 'Days grown', 'Cell line', 'KO',
                    'Library', 'Experiment ID',]
        }
    else:
        tab_columns = {
            'exp':[ 'Experiment ID', 'Treatment', 'Cell line', 'KO',  'Library', 'Citation'],
            'comp':['Comparison ID',  'Treatment', 'Dose', 'Time point group',
                    'Growth inhibition %', 'Days grown', 'Cell line', 'KO',
                    'Library', 'Experiment ID',]
        }

    selctr_tables = {
        tabk:dash_table.DataTable(
            id=f'{tabk}-table',
            # leaving it out of columns doesn't stop it from being in the table data
            columns=[column_dict(c) for c in tab_columns[tabk] if c != 'Comparison ID' ],
            data=table_dataframes[tabk].to_dict('records'),
            sort_action='native',
            sort_mode="multi",
            selected_rows=[],
            row_selectable='single',
            css=[{"selector": "p", "rule": "margin: 0"}],
            style_cell=cell_text_style,
        ) for tabk in ('exp', 'comp')
    }

    # # ** DATA TABLE **
    # # Filled with data when a comparison is selected
    idprfx_res_table = 'gene-results'
    datatable = dash_table.DataTable(
        id=idprfx_res_table+'-table',
        columns=[column_dict(x) for x in ('Effect size', 'FDR', 'Selected')],
        sort_action='native',
        sort_mode='multi',
        filter_action='native',
        # not sure what this does but I put it in the other table for some reason
        css=[{"selector": "p", "rule": "margin: 0"}],
    )

    gene_results_layout = [
        get_stat_source_selector(idprfx_res_table, 'Analysis:'),
        Div([], id=idprfx_res_table+'-title'),
        datatable
    ]


    # update the contents of the gene results tab
    @app.callback(
        Output(idprfx_res_table+'-table', 'columns'),
        Output(idprfx_res_table+'-table', 'data'),
        Output(idprfx_res_table+'-title', 'children'),
        Output('volcano-stat-source-selector', 'value'),

        Input('selected-comp', 'data'),
        Input('volc-gene-dropdown', 'value'),
        Input(idprfx_res_table+'-stat-source-selector', 'value'),
        State('volcano-stat-source-selector', 'value'),
    )
    def update_results_table_data(selected_comp, selected_genes, stat_source,
                                  volcano_stat_source, ):
        LOG.debug(f'CALLBACK: update_results_table_data({selected_comp}, {stat_source})')
        if not selected_comp:
            raise PreventUpdate

        # labels for the analysis type
        ans_lab = data_set.analysis_labels[stat_source]
        score_lab = data_set.score_labels[stat_source]

        # data for the table
        dat = data_set.get_score_fdr(stat_source, stat_source)
        lfc = data_set.get_score_fdr('mag', 'mag')['score'][selected_comp]
        score = dat['score'][selected_comp]
        fdr = dat['fdr'][selected_comp]

        # index is genes
        index = lfc.index.union(score.index).union(fdr.index)
        is_selected = index.map(lambda g: ['❌','Selected'][g in selected_genes])

        # Build the table
        results_tab = pd.DataFrame(
            {
                'Log2(FC)':lfc,
                score_lab:score,
                'FDR':fdr,
                'Gene selected':is_selected,
            },
            index=index
        )
        no_stats = results_tab[[score_lab, 'FDR']].isna().any(1)
        results_tab = results_tab.loc[~no_stats]
        results_tab.index.name = 'Gene'

        # Create the Output objects
        results_data = results_tab.reset_index().to_dict('records')

        columns = [column_dict(x) for x in results_data[0].keys()]
        treatment_label = get_treatment_label(comparisons.loc[selected_comp], ans_lab)

        treatment_para = [html.H3(f"{treatment_label[0]}"),
                          html.H4(f"{treatment_label[1]}")]

        # don't output stat_source to next selector if they all match
        if stat_source == volcano_stat_source:
            stat_source = dash.no_update

        return columns, results_data, treatment_para, stat_source


    # ** filter_dropdowns **
    # Used to filterdown the rows of the experiment and comparisons tables.
    def generate_filter_dropdowns(filter_cols, table_str):
        """Return list of div with dcc.Dropdowns with id=table_str+'-'+col for
        col in filter_cols. Options from comparisons[col]"""
        filter_dropdowns = []
        for col in filter_cols:
            filter_dropdowns.append(
                html.Div(
                    children=[dcc.Dropdown(
                        id=table_str+'-filter-'+col,
                        placeholder='Filter by '+col,
                        multi=True,
                        style={'height':'80px', 'width':'250px'},
                        value=[],
                        options=[{'label':v, 'value':v} for v in sorted(comparisons[col].unique())]
                    )],
                    style={'display':'inline-block'})
            )
        return filter_dropdowns


    filter_keys = {'exp':['Treatment', 'Cell line', 'KO', 'Library'],
                   'comp':['Treatment', 'Cell line', 'KO', 'Library', 'Experiment ID']}
    if not public_version:
        for k, l in filter_keys.items():
            l.append('Source')
    filter_dropdowns = {tabk:generate_filter_dropdowns(comptabk, tabk) for tabk, comptabk in filter_keys.items()}


    # ############## #
    # ****LAYOUT**** #
    # ############## #
    tabs = dcc.Tabs(id='tabs', value='exp-tab', children=[
        # experiment selector
        dcc.Tab(value='exp-tab', label='Select Experiment',
                className='selector-tab', selected_className='selector-tab--selected', children=[
            html.P(['Choose an experiment below to filter the options in "Select Treatment" table. '
                    'Go straight to Select Treatment to see all options.'],
                   style={'margin-top': '15px'}),
            Div(filter_dropdowns['exp'], style={'margin-bottom': '15px', }),
            Div([selctr_tables['exp']])
        ]),
        # comparison selector
        dcc.Tab(
            value='comp-tab', label='Select Treatment',
            className='selector-tab', selected_className='selector-tab--selected', children=[
                html.P(style={'margin-top': '15px'}, children=[
                    'Select a specific treatment using the table below. Click on tabs to '
                    'the right to see results for a selected treatment'
                ]),
            Div(filter_dropdowns['comp'], style={'margin-bottom': '15px', }),
            Div([selctr_tables['comp']])
        ]),
        # volcano
        dcc.Tab(
            label='Volcano plot', value='volcano-tab',
            className='data-tab', selected_className='data-tab--selected',
            children=volcano_layout
        ),
        dcc.Tab(
            label='Results Table', value=idprfx_res_table+'-tab',
            className='data-tab', selected_className='data-tab--selected',
            children=gene_results_layout
        ),

    ])

    se_layout = Div([
        html.H1("Screens explorer"),
        html.P("Select data using the blue tabs, and view results in the green tabs."),
        tabs,
        Div([html.P(id='debug')], ),

        dcc.Store(id='selected-exp', storage_type='session',),
        dcc.Store(id='selected-comp', storage_type='session'),
        dcc.Store(id='volcano-data', storage_type='session',),
    ])


    # ***CALLBACKS***
    for tabk in ('exp', 'comp'):
        table_id = f'{tabk}-table'
        @app.callback(
            Output(table_id, 'data'),
            Output(table_id, 'selected_rows'),
            # Filters from metadata columns, additional filtrs go after
            [Input(f'{tabk}-filter-{filtr_col}', 'value') for filtr_col in filter_keys[tabk]],
            State(table_id, 'selected_rows'),
            State(table_id, 'data'),
        )
        def filter_selector_table(*filters_etc):
            """Filter rows in the exp or comp tables based on value values in the
            filter boxes."""
            LOG.debug(f'CALLBACK: filter_datable()')
            selected_row = filters_etc[-2]
            table_data   = filters_etc[-1]

            if not callback_context.triggered:
                raise PreventUpdate

            selected_table = callback_context.triggered[0]['prop_id'].split('-')[0]
            assert selected_table in ('exp', 'comp')

            # filter the table
            # remove states
            filters = filters_etc[:len(filter_keys[selected_table])]
            filtered_table = table_dataframes[selected_table]
            # go through each of the filters, apply them to the table_of_comparisons
            for col, filtr in zip(filter_keys[selected_table], filters):
                if not filtr:
                    continue
                # select the rows that contain filtered values
                if selected_table == 'comp':
                    filtered_table = filtered_table.loc[
                        filtered_table[col].apply(lambda val: val in filtr)
                    ]
                else: # == 'exp'
                    filtered_table = filtered_table.loc[
                        exptab_data_searchable[col].apply(lambda vals: any([v in filtr for v in vals]))
                    ]

            # if we've switched tables, deselect rows cus the row names don't match anymore
            for trigger in callback_context.triggered:
                if trigger and trigger['prop_id'] == 'table-selector.value':
                    selected_row = []

            # If the selected row still exists in the table, maintain that selection
            #   otherwise deselect the row so we don't show some random graph
            #   (happily, deselecting prevents update so we keep the same graph).
            # Get the pre-filtered compID
            if selected_row:
                index_key = {'exp':'Experiment ID', 'comp':'Comparison ID'}[selected_table]
                correct_compid = table_data[selected_row[0]][index_key]
                # new, if it's still in the table
                if correct_compid in filtered_table.index:
                    # it's a list, cus selected_rows in datatable could be multiple
                    new_selected_row = [filtered_table.index.get_loc(correct_compid)]
                else:
                    new_selected_row = []
            else:
                new_selected_row = []

            # at least one column of comparisons is not for display, so
            filtered_table = filtered_table.loc[:, tab_columns[selected_table]]
            return (filtered_table.to_dict('records'),
                    new_selected_row)


    # Upon selecting experiment, filter the comparisons table (by setting the value of the
    #   filter dropdown) and switch to the comp table tab
    @app.callback(
        Output('tabs', 'value'),
        Output('comp-filter-Experiment ID', 'value'),
        Output('selected-exp', 'data', ),
        Input('exp-table', 'selected_rows'),
        State('exp-table', 'data'),
        State('selected-exp', 'data', ),
    )
    def picked_experiment(selected_row, table_data, previous_exp):
        LOG.debug(f'CALLBACK: picked_experiment (from exp-table): {selected_row}')
        if not selected_row:
            raise PreventUpdate
        print('picked_exp', callback_context.triggered)
        selected_exp = table_data[selected_row[0]]['Experiment ID']
        if selected_exp == previous_exp:
            raise PreventUpdate

        return ('comp-tab', [selected_exp], selected_exp)


    # Store selected comparison. This will trigger update of charts/tables
    # (in future might be multiple ways to select comp, so not attachingn outputs
    # to Input(comp-table, selected-row) is more flexible)
    @app.callback(
        Output('selected-comp', 'data' ),
        Input('comp-table', 'selected_rows'),
        State('comp-table', 'data'),
    )
    def select_comp(selected_row, table_data):

        LOG.debug(f'CALLBACK: select_comp({selected_row}')
        if not selected_row:
            raise PreventUpdate

        # get the comparison ID, select the relevant data from dataset
        compid = table_data[selected_row[0]]['Comparison ID']

        return compid


    # Output data for the volcano chart. X is LFC from mageck
    #  Y is changable source of significance
    #  Not sure this needs to be split from render_volcano
    @app.callback(
        Output('volcano-data', 'data'),
        Output('volc-gene-dropdown', 'options'),
        Output(idprfx_res_table+'-stat-source-selector', 'value'),

        Input('selected-comp', 'data' ),
        Input('volcano-stat-source-selector', 'value'),
        State(idprfx_res_table+'-stat-source-selector', 'value')
    )
    def update_volcano_data(compid, sig_source, table_stat_source):

        args_for_printing = {k:v for k, v in zip(
            'selected_row, sig_source'.split(', '),
            [compid, sig_source]
        )}
        LOG.debug(f'CALLBACK: update_volcano_data with {args_for_printing}')

        if not compid:
            raise PreventUpdate

        # get x, y and genes values
        score_fdr = data_set.get_score_fdr('mag', sig_source)
        score, fdr = [score_fdr[k][compid].dropna() for k in ('score', 'fdr')]

        # some genes may get filtered out
        unified_index = score.index.intersection(fdr.index)
        score,fdr = [xy.reindex(unified_index) for xy in (score,fdr)]

        volcano_data = {'x': score, 'y': fdr, 'genes': score.index}
        gene_options = get_lab_val(score.index)

        LOG.debug(f'End of update_volcano_data with:')
        LOG.debug('     datatable:  '+'\n'.join([f"{k}={volcano_data[k].head()}" for k in ('x', 'y')]))

        # don't output stat_source to next selector if they all match
        if sig_source  == table_stat_source:
            sig_source = dash.no_update

        return (
            volcano_data,
            gene_options,
            sig_source,
        )

    ####
    # render volc should trigger only when a comparison has been selected
    # another call back should deal with selection, look at comparison maker
    # Some more of that should be shared, like the title generator.
    @app.callback(
        Output('volcano0', 'figure'),

        Input('volc-gene-dropdown', 'value'),
        Input('volcano-data', 'data'),

        State('comp-table', 'data'),
        State('comp-table', 'selected_rows'),
    )
    def render_volcano(selected_genes, xy_genes,
                       table_data, selected_row):

        # debug device
        def count_upper():
            n = 0
            while True:
                n+=1
                yield str(n)
        counter = count_upper()
        LOG.debug(f"Rendering volcano {next(counter)}")

        if not selected_row:
            raise PreventUpdate

        # get values and create the figure
        x, fdr, genes = [xy_genes[k] for k in ('x', 'y', 'genes')]
        # dcc.Store converts this to a list...
        x, fdr = [pd.Series(xy, index=genes) for xy in (x,fdr)]
        y = fdr.apply(lambda _x: -np.log10(_x))

        fig = go.Figure(
            data=go.Scattergl(
                x=x.values,
                y=y.values,
                mode='markers',
                customdata=fdr,
                text=genes,
                hovertemplate= ("<b>%{text}</b><br>" +
                                "LFC: %{x:.2f}<br>" +
                                "FDR: %{customdata:.2e}")
            ),
            layout={'clickmode':'event+select',
                              'dragmode':'select'},
        )
        LOG.debug(f"{(len(x), len(y))}")
        LOG.debug(next(counter))
        # add some titles
        # consider using data_set.comparisons.loc[<Input('selected-comp', 'data')>]
        row = table_data[selected_row[0]]
        line1,line2 = get_treatment_label(row)
        title = f"<b>{line1}</b><br>{line2}"
        fig.update_layout(
            title=title,
            xaxis_title='LFC',
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


    # update selected genes by points that are selected on the graph
    # this should only ever add points.
    @app.callback(
        Output('volc-gene-dropdown', 'value'),
        Input('volcano0', 'selectedData'),
        State('volc-gene-dropdown', 'value')
    )
    def add_volcano_selection(selected_data, dropdown_genes):
        LOG.debug(f'add_volcano_selection, {selected_data}')

        if not selected_data:
            raise PreventUpdate

        selected_genes = set()
        for p in selected_data['points']:
            selected_genes.add(p['text'])

        if selected_genes.issubset(dropdown_genes):
            raise PreventUpdate

        return dropdown_genes+list(selected_genes.difference(dropdown_genes))

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
