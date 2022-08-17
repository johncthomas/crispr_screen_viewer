import copy

from dash import Dash

from dash.dependencies import Input, Output, State

from dash.dash_table import DataTable
from dash import dcc, html, callback_context

from crispr_screen_viewer.functions_etc import (
    cell_text_style,
    datatable_column_dict,
    get_selector_table_filter_keys,
    get_metadata_table_columns,
    LOG,
    style_hidden,
    style_comparisons_card,
    style_gene_selector_div
)

from dash.exceptions import PreventUpdate

from shared_components import (
    spawn_filter_dropdowns
)

Div = html.Div
import pandas as pd

from crispr_screen_viewer.functions_etc import (
    DataSet,
    load_dataset
)

import dash_bootstrap_components as dbc

from typing import List, Dict


def get_compid_from_rowi(table_data, rowi):
    """table_data from Input(table_id, 'data').
    rowi from Input(table_id, 'selected_rows')[n]"""
    return table_data[rowi]['Comparison ID']

def prep_exp_data(data_set:DataSet):
    """return DataFrames of (experiment_data, searchable_experiment_data). The
    latter containing lists of values to be searched against."""
    comparisons = data_set.comparisons

    LOG.debug('Comparison tab columns '+str(comparisons.columns))

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
    for col in comparisons.columns:
        for expid, rows in groups_exp.groups.items():
            vals = sorted(set(comparisons.loc[rows, col].dropna().astype(str).values))
            exptab_data[expid][col] = ', '.join(vals)
            exptab_data_searchable[expid][col] = vals

    exptab_data = pd.DataFrame(exptab_data).T
    exptab_data_searchable = pd.DataFrame(exptab_data_searchable).T
    exptab_data.index.name = 'Experiment ID'
    exptab_data_searchable.index.name = 'Experiment ID'

    return exptab_data, exptab_data_searchable


def register_filter_callbacks(app, id_prefix, exp_or_comp, filter_columns,
                              searchable_data, table_dataframes):
    """Registers a callback that filters rows by values put into filter text boxes.

    Input IDs are [f'{id_prefix}-{exp_or_comp}-filter-{col}' for in filter_columns]

    Output ID is  f'{id_prefix}-{exp_or_comp}-table'"""
    table_id = f'{id_prefix}-{exp_or_comp}-table'
    LOG.debug(
        f"Registering data selector table row filters with IDs:\n"
        f"   Output: {table_id}\n"
        f"   Inputs: {id_prefix}-{exp_or_comp}-filter-{{col}} for col in {filter_columns}"
    )

    @app.callback(
        Output(table_id, 'data'),
        Output(table_id, 'selected_rows'),
        # Filters from metadata columns, additional filtrs go after
        [Input(f'{id_prefix}-{exp_or_comp}-filter-{filtr_col}', 'value') for filtr_col in filter_columns],
        State(table_id, 'selected_rows'),
        State(table_id, 'data'),
    )
    def filter_selector_table(*filters_etc):
        """Filter rows in the exp or comp tables based on value values in the
        filter boxes."""
        LOG.debug(f'CALLBACK: filter_selector_table:\n\tTriggered by {callback_context.triggered[0]["prop_id"]}\n\tfilters ={filters_etc[:-2]};')

        selected_row = filters_etc[-2]
        table_data   = filters_etc[-1]
        #print(table_data)
        if not callback_context.triggered:
            raise PreventUpdate

        selected_table = callback_context.triggered[0]['prop_id'].split('-')[1]
        assert selected_table in ('exp', 'comp')

        # filter the table
        # remove states from the inputs
        filters = filters_etc[:len(filter_columns)]

        # go through each of the filters, apply them to the table_of_comparisons
        filtered_table = table_dataframes[selected_table]
        for col, filtr in zip(filter_columns, filters):
            if not filtr:
                continue
            #print(selected_table, col, filtr)
            # select the rows that contain filtered values
            if selected_table == 'comp':

                filtered_table = filtered_table.loc[
                    filtered_table[col].apply(lambda val: val in filtr)
                ]
            else: # == 'exp'
                filtered_table = filtered_table.loc[
                    searchable_data[col].apply(lambda vals: any([v in filtr for v in vals]))
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
            LOG.debug(f"Data of selected row: {table_data[selected_row[0]]}")
            correct_compid = table_data[selected_row[0]][index_key]
            # new, if it's still in the table
            if correct_compid in filtered_table.index:
                # it's a list, cus selected_rows in datatable could be multiple
                new_selected_row = [filtered_table.index.get_loc(correct_compid)]
            else:
                new_selected_row = []
        else:
            new_selected_row = []

        # we reindex here because it's possible that the CSV is missing columns,
        # would rather 'undefined' values than crashing
        #filtered_table = filtered_table.reindex(columns=table_columns)
        return (filtered_table.to_dict('records'),
                new_selected_row)

def register_exptable_filters_comps(app, id_prefix):
    """Register the callback that causes row selections in the exp-table
    to filter the comp table, via the '{id_prefix}-comp-filter-{Experiment ID}'.

    Current implimentation means that any change to selected rows will override
    selections made via the dropdown. Which is probably fine."""
    @app.callback(
        Output(f"{id_prefix}-comp-filter-Experiment ID", 'value'),
        Input(f"{id_prefix}-exp-table", 'selected_rows'),
        State(f"{id_prefix}-exp-table", 'data'),
        #State(f"{id_prefix}-comp-filter-Experiment ID", 'value'),
    )
    def filter_comps_by_exp(selected_rows, table_data, ):
        selected_exps = []
        for rowi in selected_rows:
            exp = table_data[rowi]['Experiment ID']
            selected_exps.append(exp)
        return selected_exps


def get_DataTable(id_prefix, table_type, columns, data) -> DataTable:
    idd = f'{id_prefix}-{table_type}-table'
    LOG.debug(f"get_DataTable(id={idd}, columns={columns})")
    selectable = {'exp':'multi', 'comp':'single'}[table_type]
    return DataTable(
        id=idd,
        # leaving it out of columns doesn't stop it from being in the table data
        columns=[datatable_column_dict(c) for c in columns if c != 'Comparison ID' ],
        data=data.to_dict('records'),
        sort_action='native',
        sort_mode="multi",
        selected_rows=[],
        row_selectable=selectable,
        css=[{"selector": "p", "rule": "margin: 0"}],
        style_cell=cell_text_style,
    )


def spawn_selector_tables(
        app, data_set, filter_columns, public_version, id_prefix,
) -> Dict[str, DataTable]:
    """Generate two DataTables, experiment and comparison (called treatment
    in parlance of the website). Returned in dict with keys 'exp' and 'comp'.
    Registers callbacks with each to handle filtering via values from the
    filter boxes."""


    LOG.debug(f"Spawning selection tables for page {id_prefix}")
    exp_data, exp_data_searchable = prep_exp_data(data_set)
    comp_data = data_set.comparisons

    table_dataframes = {'exp':exp_data, 'comp':comp_data}
    table_columns = get_metadata_table_columns(public_version, id_prefix)
    components = {}
    for tabk in table_dataframes.keys():
        components[tabk] = get_DataTable(
            id_prefix, tabk, table_columns[tabk], table_dataframes[tabk]
        )
        register_filter_callbacks(
            app, id_prefix, tabk,
            filter_columns[tabk],
            exp_data_searchable,
            table_dataframes,
        )

    register_exptable_filters_comps(app, id_prefix)

    return components

def spawn_selector_tabs(
        app, id_prefix:str, filter_dropdowns:Dict, selctr_tables:Dict,
        exp_tab_text, comp_tab_text, comp_choice_panel:list=None
) -> List[dcc.Tab]:
    """Return list of two tabs, comp and exp tabs with tables. Filters and
    tables in dicts keyed with 'exp' and 'comp'. Texts should explain how
    to use the tables, are displayed below the filters.

    Comp_choice_panel is the X/Y setter for Comparison Maker page"""

    @app.callback(
        Output(f'{id_prefix}-results-control-panel', 'style'),
        Output(f'{id_prefix}-gene-selector-div', 'style'),
        Input(f'{id_prefix}-tabs', 'value'),
    )
    def hide_comp_gene_selectors(selected_tab):
        """Set comparison/gene selector dropdown container's
        style={display: none}
        when selection tabs are active, show them when chart or table tabs
        are active."""
        if not selected_tab:
            raise PreventUpdate
        # any of these should be fine, if we add things later
        LOG.debug(f"Selected tab = {selected_tab}")
        show_treat_select = False
        for n in ('chart', 'graph', 'figure', 'table'):
            if n in selected_tab:
                show_treat_select = True
        if show_treat_select:
            return (style_comparisons_card, style_gene_selector_div)
        else:
            return (style_hidden, style_hidden)

    if comp_choice_panel is not None:
        comp_card_body = comp_choice_panel+[Div([selctr_tables['comp']])]
    else:
        comp_card_body = [Div([selctr_tables['comp']])]
    return [
        # selector tabs
        dcc.Tab(
            value=f'{id_prefix}-exp-tab', label='Select Experiment',
            className='selector-tab', selected_className='selector-tab--selected',
            children=[
                html.P(
                    # ['Choose an experiment below to filter the options in "Select Treatment" table. '
                    # 'Go straight to Select Treatment to see all options.'],
                    ' ',
                    style={'margin-top': '15px'}
                ),
                Div(filter_dropdowns['exp'], style={'margin-bottom': '15px', }),
                dbc.Card(
                    [
                        dbc.CardHeader(exp_tab_text),
                        dbc.CardBody(Div([selctr_tables['exp']]))
                    ]
                )
            ]
        ),
        # comparison selector
        dcc.Tab(
            value=f'{id_prefix}-comp-tab', label='Select Treatments',
            className='selector-tab', selected_className='selector-tab--selected',
            children=[
                html.P(
                    style={'margin-top': '15px'},
                    children=[' '
                              # 'Sele'+'ct a specific treatment using the table below. Choose to show scores from a '
                              #        'treatment on the X or Y axis using the buttons below.'
                              ]
                ),
                Div(filter_dropdowns['comp'], style={'margin-bottom': '15px', }),
                dbc.Card(
                    style={'padding-top':'10px'},
                    children=[
                        dbc.CardHeader(comp_tab_text),
                        dbc.CardBody(comp_card_body)
                    ]
                )
            ]
        )
    ]

def spawn_treatment_reselector(id_prefix, is_xy:bool) -> Div:
    """Dropdowns that float next to figures and tables for reselecting
     previously selected comps. Dropdowns with ID '{id_prefix}-{axis}-selector'
     are the component that decides the comp from which data is displayed.

     {axis} is 'x'/'y' for CM page, and 'comp' for SE page.

     Callback that hides this based on selected tab is registered when as
     part of selector-tab spawning. Callbacks using dropdown(s) value as
     Input defined in page scripts as they're pretty specific."""
    def get_xy_selectors(id_prefix, is_xy:bool):

        ilbstyle = {'display':'block'}

        output_selectors = []

        if is_xy:
            dds ='xy'
        else:
            dds = ('comp',)

        for xy in dds:
            label = {'x':'X treatment',
                     'y':'Y treatment',
                     'comp':'Treatment'}[xy]
            id=f'{id_prefix}-{xy}-selector'
            xy_selector = dcc.Dropdown(
                id=id,
                placeholder=f"Select {label}",
                className='comp-selector',
                options=[],
            )
            label = html.Label(f'{label}:',htmlFor=id,style=ilbstyle)
            block = Div([label, xy_selector,])
            output_selectors.append(block)

        return output_selectors


    comparison_selectors_card = Div(
        id=f'{id_prefix}-comparison-selectors-div',
        children=[
            dbc.Card(
                [
                    dbc.CardHeader(
                        "Previously selected treatments are available below."
                    ),
                    dbc.CardBody(
                        get_xy_selectors(id_prefix, is_xy)
                    )
                ]
            )
        ],
    )

    return comparison_selectors_card

def test(testing_page='cm'):
    """testing_page == 'cm' or 'se'"""
    PAGE_ID = testing_page
    p = '/Users/johnc.thomas/Dropbox/crispr/DDRcs/app_data/toy_data'
    dataset = load_dataset(p)
    comparisons = dataset.comparisons
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

    # ** filter_dropdowns **
    # Used to filterdown the rows of the experiment and comparisons tables.
    filter_keys = get_selector_table_filter_keys(public=True, )
    filter_keys = copy.copy(filter_keys)
    filter_keys['comp'].append('Timepoint')
    filter_dropdowns = {tabk: spawn_filter_dropdowns(PAGE_ID, tabk, comptabk, comparisons)
                        for tabk, comptabk in filter_keys.items()}

    selctr_tables = spawn_selector_tables(
        app, dataset, filter_keys, public_version=True, id_prefix=PAGE_ID,
    )

    def get_xy_choice_panel(xy:str) -> List[Div]:
        """Return button and output Divs.
        Button: id="cm-choose-{xy.lower()}"
        Para:   id="cm-chosen-{xy.lower()}"
        """
        xy = xy.lower()
        XY = xy.upper()
        xybutton = html.Button(f'Choose {XY} treatment', id=f'{PAGE_ID}-choose-{xy}', n_clicks=0)
        return [
            Div([xybutton]),
            Div(html.P(
                id=f"{PAGE_ID}-chosen-{xy}",
                children=f'No {XY} treatment selected. Choose from list and press "Choose" button'
            ))
        ]

    tabs = dcc.Tabs(id='tabs', value=f'exp-tab', children=[
        # experiment selector
        dcc.Tab(value=f'{PAGE_ID}-exp-tab', label='Select Experiment',
                className='selector-tab', selected_className='selector-tab--selected', children=[
                html.P(['Choose an experiment below to filter the options in "Select Treatment" table. '
                        'Go straight to Select Treatment to see all options.'],
                       style={'margin-top': '15px'}),
                Div(filter_dropdowns['exp'], style={'margin-bottom': '15px', }),
                Div([selctr_tables['exp']])
            ]),
        # comparison selector
        dcc.Tab(
            value=f'{PAGE_ID}-comp-tab', label='Select Treatments',
            className='selector-tab', selected_className='selector-tab--selected', children=[
                html.P(style={'margin-top': '15px'}, children=[
                    'Select a specific treatment using the table below. Click on tabs to '
                    'the right to see results for a selected treatment'
                ]),
                Div(filter_dropdowns['comp'], style={'margin-bottom': '15px', }),
                Div(get_xy_choice_panel('x')),
                Div(get_xy_choice_panel('y')),
                Div([selctr_tables['comp']])
            ]
        ),
        dcc.Tab(
            value='selected', label='Selected', children=[
                html.P(id='selected-para')
            ]
        )
    ])

    app.layout = tabs

    # for updating what the chosen treatments are
    @app.callback(
        Output(f'{PAGE_ID}-chosen-x', 'children'),
        Output(f'{PAGE_ID}-chosen-y', 'children'),

        Input(f'{PAGE_ID}-choose-x', 'n_clicks'),
        Input(f'{PAGE_ID}-choose-y', 'n_clicks'),
        State(f'{PAGE_ID}-comp-table', 'selected_rows'),
        State(f'{PAGE_ID}-chosen-x', 'children'),
        State(f'{PAGE_ID}-chosen-y', 'children'),
        State(f'{PAGE_ID}-comp-table', 'data'),
    )
    def select_treat_for_cm(
            xbutton, ybutton, selected_row, state_x, state_y, table_data
    ):
        try:
            triggered = callback_context.triggered_id
        except AttributeError:
            # v<2.4
            triggered = callback_context.triggered[0]['prop_id'].split('.')[0]
        if not triggered:
            raise PreventUpdate

        xy = triggered.split('-')[-1]

        compid = get_compid_from_rowi(table_data, selected_row[0])
        LOG.debug(f'Choosing {compid} for axis {xy}')
        if xy == 'x':
            return [compid, state_y]
        else:
            return [state_x, compid]

    register_exptable_filters_comps(app, PAGE_ID)

    app.run_server(debug=True)

if __name__ == '__main__':
    LOG.setLevel('DEBUG')
    test('cm')