# select any two samples that share some genes and analysis method
#   Warnings if they are inappropriate comparisons, sample deets should be visable on the screen
#   Maybe make the table a seperate tab
# quantify distance by euc, mageck fdr, JACKS p difference as appropriate
# distance from x=y or regression
#   regression uses y-y' as the distance

# dropdown to select distance type
# checkbox to do regression (performed at the time if not available)

# to check for the existence of MAGecK derived p values we need to access
# the origin of these, sample identity must include screen information

# todo sample options don't seemt ot properly update when changing ctrl groups
# todo?? make it possible to save the state of the whole thing

from copy import copy
from functools import partial
import dash
import dash_table
from dash.exceptions import PreventUpdate
from dash_table.Format import Format, Scheme
import dash_core_components as dcc
import dash_html_components as html
Div = html.Div
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
from typing import Iterable, List, Union, Tuple, Dict
import numpy as np
import statsmodels.api as sm
from scipy import stats

#from . import util


from crispr_screen_viewer.util import tabulate_score, tabulate_mageck, load_mageck_tables

import pandas as pd
from pandas import MultiIndex


def valid_comparisons(ctrldict:dict, lax=True) -> dict:
    """pass {<ctrlgrp>: [{<ctrlsamp:[samps]}, ...],
             ...}

    returns {"<ctrlgrp>-<ctrlsamp>:[<samps>, ...], ...}
    or (if lax) {"<ctrlgrp>":[<samps>, ...], ...}

    i.e. the 'controls' dict from expdict/yaml

    These should be used to construct the sample selector in the scatter
    function.

    if `lax` is True then all samples in the same control group will be
    assigned to the same key
    if `lax` is False then samples will only be assigned to a specific shared
    sample.

    """

    outdict = {}

    if not lax:
        for ctrlgrp, sampd in ctrldict.items():
            for ctrlsamp, samps in sampd.items():
                # Sample specific
                outdict[f"{ctrlgrp}-{ctrlsamp}"] = samps
    else:
        for ctrlgrp, sampd in ctrldict.items():
            comps = []
            for ctrlsamp, samps in sampd.items():
                comps.extend(samps)
            # control group specific.
            outdict[f"{ctrlgrp}"] = comps


    return outdict



def get_jacks_stats(df, sampx, sampy):
    """returns ps, fdr and size of difference between X and Y distributions.
    Negative values for those below the line. Is log10.

    a dataframe with columns:
        'p-value', 'FDR', 'Difference'
        """
    # if you add more, update STAT_KEYS
    X = df[sampx]
    Y = df[sampy]
    # We can assume the subcol names i guess
    diff = abs(X.jacks_score - Y.jacks_score)
    sumstd = X.stdev + Y.stdev
    # get bool array, index will be
    negative = (Y.jacks_score < X.jacks_score).values
    # cdf value at zero
    ps = stats.norm(diff, sumstd).cdf(0)
    fdr = sm.stats.multipletests(ps, 0.1, 'fdr_bh')[1]
    ps, fdr = [-np.log10(s) for s in (ps, fdr)]
    ps[negative] = 0-ps[negative]
    fdr[negative] = 0-fdr[negative]

    DISTANCES = {
        'p-value': ps,
        'FDR': fdr,
        'Difference': diff
    }

    return pd.DataFrame(DISTANCES)




def get_mageck_stats(df:pd.DataFrame, sampx:str, sampy:str, extra:pd.DataFrame) -> pd.DataFrame:
    """get DF with p & FDR from the control to x & t, and from the EXTRA
    analyses if present"""
    # select columns that contain the two samples of interest, if present
    stats_df = pd.DataFrame(index=df.index)
    stat_ks = 'fdr_log10', 'p_log10'
    for stat in stat_ks:
        stats_df.loc[:, f"ctrl->X {stat}"] = df[(sampx, stat)]
        stats_df.loc[:, f"ctrl->Y {stat}"] = df[(sampy, stat)]
    comp_present = extra.columns.levels[0][extra.columns.levels[0].isin([f"{sampx}-{sampy}", f"{sampy}-{sampx}"])]
    if comp_present.shape[0] > 0:
        if f"{sampx}-" in comp_present[0]:
            order = "X->Y"
        elif f"{sampy}-" in comp_present[0]:
            order = "Y->X"
        else:
            raise ValueError("Something wrong with the EXTRA column headers " + str(comp_present))

        for stat in stat_ks:
            stats_df.loc[:, f"{order} {stat}"] = extra[(comp_present[0], stat)]

    return stats_df


# put selections in df with single leveled column index
def df_with_xy_scores(df, xy_samples:Iterable[str]):
    """A DF with `x/y {stat}` for each column in the df, except
    Gene that will only appear once."""
    nu_df = {'Gene':df.loc[:, ('Gene', 'Gene')]}
    for stat_hdr in df.columns.levels[1]:
        if stat_hdr == 'Gene':
            continue
        for xy, samp_hdr in zip('xy', xy_samples):
            # dash is alt-`–`, en-dash
            #print(samp_hdr, stat_hdr)
            nuhdr = xy+' '+stat_hdr
            nu_df[nuhdr] = df.loc[:, (samp_hdr, stat_hdr)]
    return pd.DataFrame(nu_df)

def get_annotation_dicts(df, xkey, ykey, annote_kw=None):
    """A list of dict defining annotations from specified columns
    of a dataframe."""
    annotations = []
    if annote_kw is None:
        annote_kw = {}
    # NOTE x and y swapped, I don't know why
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
            'name': nice_headers[hdr],
            # formatting isn't working' here so doing it with strings
            # 'type': 'numeric',
            # 'format': Format(precision=2, scheme=Scheme.decimal)
        }



def get_opt_dict(values):
    return [{'label':v, 'value':v} for v in values]



##########
## so... I basically wrote this as a single script, then wrapped it into a function
## to avoid issues with globals. Then I've moved some subfunction
## definitions out of this spawn function, but not others.
##   Basically, it's a mess
def spawn_scatter(tables:Dict[str, pd.DataFrame], analysis_type:str, expd:dict, lax=True,
                  distance_filter=0):
    """"""


    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    if analysis_type == 'jacks':
        SCOREKEY = 'jacks_score'
        STAT_KEYS = ['p-value', 'FDR', 'Difference']
        get_stats = get_jacks_stats
    elif analysis_type == 'mageck':
        print(tables.keys())
        get_stats = partial(get_mageck_stats, extra=tables['EXTRA'])
        SCOREKEY = 'lfc'
    else:
        raise ValueError("Unknown analysis type "+str(analysis_type))

    # the table component ignores the index labels.
    # **this `tab` instance is used in later code.**
    for tab in tables.values():
        tab.loc[:, ('Gene', 'Gene')] = tab.index

    # AVAILABLE_SAMPLES = DF.columns.levels[0]
    # print(AVAILABLE_SAMPLES)
    #selected_samples = ('nt22a', '10gy')

    valid_comps = valid_comparisons(expd['controls'], lax)

    class COMPONENTS():
        pass
    ######################
    ## **********
    ## COMPONENTS
    ## **********
    graph_config = {'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d',
                                               'resetScale2d'],
                    'editable':True}
    graph = dcc.Graph(
        id='scatter',
        config=graph_config,
        style={'width':'900px', 'height': '800px', 'padding':'0px', },
        figure={'layout':{'clickmode':'event+select'}},
    )


    def sample_dropdown(n):
        return dcc.Dropdown(
            id='sample-dropdown-'+str(n),
            #options=[{'label': v, 'value': v} for v in AVAILABLE_SAMPLES],
            options=[],
            value=None,  # available_samples[0],
            placeholder=f'Select a {n} sample',
        )

    ctrl_dropdown = dcc.Dropdown(
        id='ctrl-dropdown',
        options=get_opt_dict(
            valid_comps.keys()
        ),
        value=None,
        placeholder='Select control group'
    )

    def selection_panel(attr,val=0):
        """Tool panel with distance type dropdown, num boxes/ranger for minimums.
        Used to control the appearence of color gradient or labels.

        `attr` defines the type, labels based on type..."""
        attr = str(attr)
        # select panel label
        labels = dict(gradient = ['Colour gradient', 'Min value', 'Max value'],
                      labels   = ['Labeling threshold', 'Deplet. threshold', 'Enrich threshold'])
        attr_label, lowlab, hilab = labels[attr]


        distdrop = dcc.Dropdown(
            id='dist-dropdown-' + attr,
            options=[{'label':'Disabled', 'value':'Disabled'}],
            value='Disabled',
            clearable=False,
            style={'height':'25px', 'width':'95%', 'display':'inline-block'},
            placeholder='Distance metric'
        )

        low_slider = dcc.Slider(
            id='lowslider-' + attr,
            min=0,
            max=100,
            step=0.01,
            value=val
        )

        hi_slider = dcc.Slider(
            id='hislider-' + attr,
            min=0,
            max=100,
            step=0.01,
            value=100
        )

        low_box = dcc.Input(
            id='lowbox-' + attr,
            type='number', #+ attr,
            size='5', # in char i think

            #inputMode='numeric',
            value=0,
        )
        hi_box = dcc.Input(
            id='hibox-' + attr,
            type='number', # + attr,
            size='5',

            #inputMode='numeric',
            value=0,
        )

        panel = Div([
            Div(
                [
                    html.P([attr_label],), #style=dict(position='absolute', top=0, left=0)),
                    Div([distdrop],), #style=dict(position='absolute', top=20, left=0, width="100px", height="20px")),
                    Div([html.P([lowlab],),low_box, low_slider], ),#style=dict(position='relative')),#, top=60, left=0,) ),
                    Div([html.P([hilab],),hi_box, hi_slider ], ),#style=dict(position='absolute', top=90,)),
                    #Div(ranger, style=dict(position='absolute', top=50, left=0, width="100%"))
                ],
                style=dict(position='relative', padding="10px")
            )],
            style=dict(width="250px",  position="relative", float='left', top="30px", padding="20px", display='inline-block'),
            id="panel-" + attr
        )

        return panel

    ########
    ## Define the DataTable
    #
    # Define the column headers for the table in the Dash format (list of dicts),
    # "x/y <stat>" for stats in the table
    column_headers = []
    # using the tab instance that was created in a loop way above
    for stat_hdr in tab.columns.levels[1]:
        if stat_hdr == 'Gene':
            continue
        column_headers.append('x ' + stat_hdr)
        column_headers.append('y ' + stat_hdr)

    column_headers = [
        {'id':h, 'name':h, 'type':'numeric','format':Format(precision=2)} for h in column_headers
    ]
    column_headers = [{'id':h, 'name':h} for h in ('Selected', 'Gene')]+column_headers

    empty_datatable = dash_table.DataTable(
        id='table0',
        columns=column_headers,
        data=[],
        sort_action="native",
        sort_mode="multi",
        #pagination_settings={"page_size": 20,},

    )

    def gene_dropdown(selected_genes):
        dropdown = dcc.Dropdown(
            id='gene-dropdown',
            options=[{'label': v, 'value': v} for v in tab.index],
            placeholder='Select genes by name',
            multi=True,
            style=dict(height='200px', width="100%"),
            value=selected_genes
        )
        return dropdown

    ##################
    ## ******
    ## LAYOUT
    ## ******
    class LAYOUT():
        pass

    samp_style=dict(width="45%", position="relative", float="left", padding="2%",
                   display="inline-block")

    main_chunk = Div([
        Div(graph, style=dict(display='inline-block', float='left')),
        selection_panel('labels', val=100),
        selection_panel('gradient'),
        Div(
            gene_dropdown([]),
            style=dict(width="350px", position="relative", float='left',
                       top="30px", padding="20px",display='inline-block'),
            id="genes-box"
        ),
        html.Div([empty_datatable],
                 style={'width':'900px', 'display':'inline-block'}),
        Div([], id='selected-genes',   style={'display':'none'}),
        Div([], id='selected-samples', style={'display':'none'}),
    ])

    app.layout = Div(
        [
            Div(ctrl_dropdown),
            Div(sample_dropdown('x'), style=samp_style),
            Div(sample_dropdown('y'), style=samp_style),
            main_chunk
        ]
    )


    class CALLBACKS():
        pass
    ##############
    ## *********
    ## CALLBACKS
    ## *********
    # Panel distance type and slider callback:
    #   When a new distance type is selected, the values in the top/bottom number boxes
    #   is set to reflect the maximum values for that distance.
    #
    #   For the gradient the default min is zero, and max is max(max(pos), max(neg))
    #   For the labels it is the maximum value in each direction
    #
    #   I guess they get one callback each since they output different things
    @app.callback(
        [Output('lowbox-gradient', 'value'), Output('hibox-gradient', 'value'),
         Output('lowslider-gradient', 'max'), Output('hislider-gradient', 'max'),],

        [Input('dist-dropdown-gradient', 'value'),
         Input('lowslider-gradient', 'value'),
         Input('hislider-gradient', 'value'),],

        [State('lowbox-gradient', 'value'), State('hibox-gradient', 'value'),
         State('ctrl-dropdown', 'value'), State('selected-samples', 'children')]
    )
    def set_gradient_range(dist_type, _lownotused, _hinotused, lowbox, hibox,
                           ctrl_samp, selected_samples):

        if dist_type == 'Disabled' or dist_type is None:
            return [0,0,1,1]

        # return current box values if they aren't updated
        boxout = {'hi': hibox, 'low': lowbox}
        triggers = dash.callback_context.triggered

        # if distance type has been changed, set hi box to maximum value
        # else, figure out which slider has moved and update boxes
        ctrl_grp = ctrl_samp.split('-')[0]
        samp_x, samp_y = selected_samples
        df = tables[ctrl_grp]
        d = get_stats(df, samp_x, samp_y)[dist_type]
        mx = max(np.abs(d))
        if 'dist-dropdown-gradient.value' in [trig['prop_id'] for trig in triggers]:
            boxout['hi'] = mx
        else:
            boxout = process_triggers(triggers,
                                  '{}slider-gradient.value',
                                  boxout)

        return [boxout['low'], boxout['hi'], mx, mx]


    def process_triggers(triggers, keystr, boxout):
        """Return first value from triggered props that match keystr.
        Used by set_<thing>_range() callbacks."""

        # format of triggers is list of dict with 'prop_id' and 'value' keys
        # Need to deal with possibility that len(triggers) > 1 (not sure it ever would be, but...)

        #print('triggers', triggers)
        trig_props = [trig['prop_id'] for trig in triggers]
        for hilow in 'hi', 'low':
            trig_k = keystr.format(hilow)
            if trig_k not in trig_props:
                continue
            ti = trig_props.index(trig_k)
            nubox = triggers[ti]['value']
            #print('skfuhsidfuh', nubox)
            boxout[hilow] = nubox

        return boxout

    @app.callback(
        [Output('lowbox-labels', 'value'), Output('hibox-labels', 'value'),
         Output('lowslider-labels', 'max'), Output('hislider-labels', 'max'), ],

        [Input('dist-dropdown-labels', 'value'),
         Input('lowslider-labels', 'value'),
         Input('hislider-labels', 'value'), ],

        [State('lowbox-labels', 'value'), State('hibox-labels', 'value'),
         State('ctrl-dropdown', 'value'), State('selected-samples', 'children')]
    )
    def set_label_range(dist_type,  _lownotused, _hinotused, lowbox, hibox,
                        ctrl_samp, selected_samples):
        """low label threshold is a postive number that will be interpreted by
        the render function to match negative distances"""
        if dist_type == 'Disabled' or dist_type is None:
            return [100,100,100,100]

        # return current box values if they aren't updated
        boxout = {'hi': hibox, 'low': lowbox}
        triggers = dash.callback_context.triggered

        # if distance type has been changed, set hi box to maximum value
        # else, figure out which slider has moved and update boxes
        ctrl_grp = ctrl_samp.split('-')[0]
        samp_x, samp_y = selected_samples
        df = tables[ctrl_grp]
        d = get_stats(df, samp_x, samp_y)[dist_type]
        #print(d)
        mx_hi = max(d)
        mx_lo = abs(min(d))
        if 'dist-dropdown-labels.value' in [trig['prop_id'] for trig in triggers]:
            boxout['hi'] = mx_hi
            boxout['low'] = mx_lo
        else:
            boxout = process_triggers(triggers,
                                  '{}slider-labels.value',
                                  boxout)
        return [boxout['low'], boxout['hi'], mx_lo, mx_hi]


    @app.callback(
        [Output('sample-dropdown-x', 'options'),
         Output('sample-dropdown-y', 'options')],
        [Input( 'ctrl-dropdown', 'value')]
    )
    def set_sample_options(ctrl_samp):

        if ctrl_samp is None:
            raise PreventUpdate

        opts  = valid_comps[ctrl_samp]

        # if we have labels, use them, otherwise just use internal samp names
        try:
            fancy_opts = []
            for o in opts:
                # todo deal with no labels (indexError here)
                lab = expd['labels'][o]
                fancy_opts.append(
                    {'label': "{} ({})".format(lab, o), 'value': o}
                )
            return fancy_opts, fancy_opts
        except KeyError:
            if expd['labels']:
                print('probably a fucked up label,', o)
            return [get_opt_dict(o) for o in opts], [get_opt_dict(o) for o in opts]



    @app.callback(
        [Output('selected-samples', 'children'),
         Output('dist-dropdown-gradient', 'options'),
         Output('dist-dropdown-labels', 'options')],

        [Input('sample-dropdown-x', 'value'),
         Input('sample-dropdown-y', 'value')],

        [State('ctrl-dropdown', 'value')]
    )
    def select_samples(samp_x, samp_y, ctrl_samp):
        """Distance stat options are set here as with mageck not every stat is
        available for every sample pair."""

        #print('** select_samples callback')

        if not samp_y or not samp_x:
            raise PreventUpdate

        if samp_x == samp_y:
            raise PreventUpdate

        #todo: we also call get_mageck_stats in render_scatter, could probably not repeat that
        if analysis_type == 'mageck':

            keys = get_mageck_stats(tables[ctrl_samp], samp_x, samp_y, tables['EXTRA']).columns
            #print('**** ', samp_x, samp_y, keys)
        else:
            keys = STAT_KEYS
        dist_opts = get_opt_dict(keys)

        return [[samp_x, samp_y], dist_opts, dist_opts]


    # Main chart callback
    panel_tups = []
    for attr in 'gradient', 'labels':
        for s in "hibox- lowbox-".split(' '):
            panel_tups.append((s+attr, 'value'))
    @app.callback(
        [Output('scatter', 'figure'),
         Output('selected-genes', 'children'),
         Output('genes-box', 'children'),
         Output('table0', 'data'),],

        [Input('selected-samples', 'children'),
         Input('scatter', 'selectedData')]+[Input(*t) for t in panel_tups],

        [State('dist-dropdown-gradient', 'value'),
         State('dist-dropdown-labels', 'value'),
         State('ctrl-dropdown', 'value'),
         State('selected-genes', 'children')]
    )
    def render_scatter(selected_samples, selectedData,
                       hi_gradient, low_gradient,
                       # inputs
                       hi_lab_lim, low_lab_lim,
                       # states
                       dist_type_gradient,
                       dist_type_labels,
                       ctrl_samp,
                       selected_genes
                       ):

        #print('** render_scatter callback')

        if not all(selected_samples) or not selected_samples or not ctrl_samp:
            raise PreventUpdate
        samp_x, samp_y = selected_samples
        ctrl_grp = ctrl_samp.split('-')[0]
        df = tables[ctrl_grp]
        # get [2]series of x/y scores o
        scores = [df[k][SCOREKEY] for k in (samp_x, samp_y)]
        scores_tab = pd.DataFrame({k:df[k][SCOREKEY] for k in  (samp_x, samp_y)})


        distances = get_stats(df, samp_x, samp_y)

        if distance_filter:

            # get mask for filtering, based on distance from median, maintaining previously
            # selected genes.
            d = scores[0] - scores[1]
            md = np.median(d)
            dist_mask = (d > md + distance_filter) | (d < md - distance_filter) | d.index.isin(selected_genes)

            # filter. This doesn't effect the table
            scores = [s.loc[dist_mask] for s in scores]





        ctx = dash.callback_context
        # print(ctx.triggered)

        # ctx triggered events formated like
        #[{'prop_id': 'dist-dropdown-gradient.value', 'value': 'Difference'}]

        # deal with colour gradient
        marker_dict = dict(size=9,
                            line=dict(width=1),
                            color='#4da6ff',
                            #colorscale='Viridis',
                            showscale=True)

        if dist_type_gradient is not None \
                and dist_type_gradient != 'Disabled' \
                and dist_type_gradient in distances.columns:

            #print(hi_gradient, low_gradient)
            hi_gradient, low_gradient = float(hi_gradient), float(low_gradient)
            # Gradient uses the absolute values, distance from zero
            gradient = np.abs(distances[dist_type_gradient])
            gradient[gradient > hi_gradient] = hi_gradient
            above_min = gradient > low_gradient

            # coloured with gradient
            grad_marker_dict = copy(marker_dict)
            grad_marker_dict['color'] = gradient[above_min]
            grad_marker_dict['colorscale'] = 'Viridis'

            # not coloured, below min value set by user
            dull_marker_dict = copy(marker_dict)
            dull_marker_dict['color'] = '#bfbfbf'

            # plot uncoloured and coloured
            g = [
                go.Scatter(x=scores[0][~above_min],
                           y=scores[1][~above_min],
                           mode='markers',
                           marker=dull_marker_dict,
                           text=df.index),
                go.Scatter(x=scores[0][above_min],
                           y=scores[1][above_min],
                           mode='markers',
                           marker=grad_marker_dict,
                           text=df.index)
            ]
        else:
            #gradient = np.array([])
            g = go.Scatter(
                x=scores[0],
                y=scores[1],
                mode='markers',
                marker=marker_dict,
                text=df.index
            )
            g = [g]


        # decipher what is selected, add annotations to the layout
        # from the low/hi boxes and from box selection
        # These alter the annotation via layout, and updates selected-genes.children
        #   which is used to update selectedData when the gene-dropdown callback happens
        #   (gene-dropdown callback always happens because the component is recreated by
        #   this callback).
        genes_to_label = []
        if selectedData and selectedData is not None:
            box_selected_genes = [x['text'] for x in selectedData['points']]
            genes_to_label.extend(box_selected_genes)

        if dist_type_labels is not None and dist_type_labels != 'Disabled':
            low_lab_lim, hi_lab_lim = float(low_lab_lim), float(hi_lab_lim)
            lab_dist = distances[dist_type_labels]
            # things below the -low limit, and above the hi limit should be labeled
            low_labels = df.index[lab_dist < -low_lab_lim]
            hi_labels  = df.index[lab_dist > hi_lab_lim]
            for labels in (low_labels, hi_labels):
                genes_to_label.extend(labels)

        genes_to_label = list(set(genes_to_label))

        #score_df = df.loc[genes_to_label, (slice(None, None), SCOREKEY)]

        layout = go.Layout(
            annotations=get_annotation_dicts(scores_tab.loc[genes_to_label], *scores_tab.columns),
        )

        # Finally, the table

        table_df = df_with_xy_scores(df, selected_samples)

        table_df.loc[:, 'Selected'] = 'No'
        table_df.loc[genes_to_label, 'Selected'] = '!Yes'

        cols = list(table_df.columns)
        #print(cols)
        # put is selected first
        table_df = table_df.reindex(columns=[cols[-1]] + cols[:-1])

        #print('********')
        return (
            dict(data=g, layout=layout),
            sorted(genes_to_label),
            gene_dropdown(genes_to_label),
            table_df.to_dict("rows")
        )


    @app.callback(
        Output('scatter', 'selectedData'),
        [Input('gene-dropdown', "value")],
        [State('scatter', 'selectedData'), State('selected-samples', 'children'),
         State('ctrl-dropdown', 'value'), State('selected-genes', 'children'),
         State('scatter', 'figure'), ]
    )
    def gene_drop_callback(dropdown_values, selectedData, selected_samples,
                          ctrl_samp, selected_genes,
                          chart_data,):
        #print('dropdown callback')

        if selectedData is None:
            selectedData = {'points':[], 'range':{'x':[0,0], 'y':[0,0]}}

        # don't update if nothing has changed (prevent recurssion)
        if sorted(dropdown_values) == sorted(selected_genes):
            raise PreventUpdate

        ctrl_grp = ctrl_samp.split('-')[0]
        df = tables[ctrl_grp]

        new_selected_points = []
        for gn in dropdown_values:
            dati = chart_data['data'][0]['text'].index(gn)
            xy = df.loc[gn, (selected_samples, SCOREKEY)]
            point = dict(
                curveNumber=0,
                pointNumber=dati,
                pointIndex=dati,
                x=xy[selected_samples[0]],
                y=xy[selected_samples[1]],
                text=gn,
            )
            new_selected_points.append(point)
        selectedData['points'] = new_selected_points
        #print('****end dropdown callback')
        return selectedData

    return app
    #app.run_server(debug=True, port=8051)

# spawn_scatter(
#     tabulate_score('/Users/johnc.thomas/thecluster/jct61/database/screens/ramsay_528-620/take1_all/jacks_median/files/ram_528-620.pretreat.'),
#     'jacks'
# )



if __name__ == '__main__':

    import yaml
    #_expd = yaml.safe_load(open('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/david_756-7^ddrV2/dav756-7.yaml'))
    #f = '/Users/johnc.thomas/Dropbox/crispr/screens_analysis/david_756-7^ddrV2/dav_756-7/take1-firstruns/jacks_median/files/dav_756-7.'
    _expd = yaml.safe_load(open('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/ramsay_759+/ram_759-80.3.repmap.yaml'))

    _analysis_type = 'mageck'


    if _analysis_type == 'jacks':
        f = '/Users/johnc.thomas/Dropbox/crispr/pkg/dash_charts/ram_759-80/take3/jacks_median/files/ram_759-80.'
        _tables = {}
        for _ctrl_grp in _expd['controls'].keys():
            _tables[_ctrl_grp] = tabulate_score(f+_ctrl_grp+'.')
    elif _analysis_type == 'mageck':
        f = "/Users/johnc.thomas/Dropbox/crispr/pkg/dash_charts/ram_759-80/take3/mageck/files/ram_759-80."
        _tables = load_mageck_tables(f, list(_expd['controls'].keys())+['EXTRA'])
    else:
        exit('pick an analysis type')


    app = spawn_scatter(_tables, _analysis_type, _expd)
    app.run_server(debug=True, port=8051)

