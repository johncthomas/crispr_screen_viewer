import copy
import inspect
import typing
import unicodedata

import loguru
#from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
import os
from pathlib import Path
from scipy.stats import norm
from scipy import odr
from scipy.stats import linregress
from typing import Union, List, Dict, Iterable, Collection, Sequence
from loguru import logger
from dash import html
parse_expid = lambda comp: comp.split('.')[0]
import sys


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

timepoint_labels = dict([
    ('fromstart', 'From experiment start'),
    ('otherprior', 'From midpoint'),
    ('endpoints', 'Matched time points')
])


# Styles for things that get hidden - these should not be transferred to
# .css file because they are dynamically changed.
style_comparisons_card = {'padding-top':'98px',
                          'display':'inline-block',
                          #'border':'2px red dotted',
                          'width':'350px'}

style_hidden = {'display':'none'}
style_gene_selector_div = {}

def set_loguru_level(logger, level='INFO'):
    logger.remove(0)
    logger.add(sys.stderr, level=level)

def get_resource_path(relative_path):
    fn = os.path.join(os.path.dirname(__file__), relative_path)
    return fn

def load_stats_csv(
        fn,
        drop_controls_with=('control', 'Control', 'CONTROL', 'Non-target', 'Cutting'),
) -> pd.DataFrame:
    """Load the double-headered stats CSV. Drop genes with any substrings
    defined in `drop_controls_with`."""
    df = pd.read_csv(fn, index_col=0, header=[0,1])
    if df.index.isna().any():
        df = df.loc[~df.index.isna()]
        logger.warning(f"NaN genes found in (and dropped from) {fn}")

    if drop_controls_with:
        # get Falses
        drop_mask = df.index != df.index
        # Make True where index contains a control indicating substring
        for ctrl_indicator in drop_controls_with:
            drop_mask = drop_mask | df.index.str.contains(ctrl_indicator)
        df = df.loc[~drop_mask]

    return df



def normalise_text(s:str):
    """Removes formatting bytes, literally `unicodedata.normalize('NFKD' s)`"""
    return unicodedata.normalize('NFKD', s)

def datatable_column_dict(c,):
    """return {'name':k, 'id':k} for all k except "DOI" which
    sets style to markdown. List of returned dicts to be passed to
    DataTable columns arg.

    The idea is that more columns can get specific options that are shared
    across all the tables in all the tools."""
    markdown=('DOI',)
    if c in markdown:
        return {"name": c, "id": c, 'type': 'text', 'presentation':'markdown',}
    else:
        return {'name':c, 'id':c}

def doi_to_link(doi):
    """Return string formated as a Markdown link to doi.org/{doi}"""
    # should be formated to not be a hyperlink, but sometimes it is
    if pd.isna(doi) or (not doi):
        return ''
    doi = doi.replace('https', '').replace('http', '').replace('://', '').replace('doi.org/', '')
    return f"[{doi}](https://doi.org/{doi})"

def orthoregress(x, y):
    """Orthogonal Distance Regression.

    Returns: [slope, offset]"""

    def f(p, x):
        return (p[0] * x) + p[1]

    model = odr.Model(f)
    data = odr.Data(x, y)
    od = odr.ODR(data, model, beta0=linregress(x, y)[0:2])
    res = od.run()

    return list(res.beta)


def index_of_true(bool_mask):
    return bool_mask[bool_mask].index


def get_selector_table_filter_keys(public=False) -> Dict[str, List[str]]:
    """Get dictionary used to filter the exp and comp tables. """
    filter_keys = {
        'exp':['Treatment', 'Cell line', 'KO', 'Library'],
        'comp':['Treatment', 'Cell line', 'KO', 'Timepoint',
                'Library', 'Citation',]
    }
    if not public:
        for k, l in filter_keys.items():
            l.append('Source')
    return filter_keys


def get_metadata_table_columns(public, page_id) -> Dict[str, List[str]]:
    # this is set up so that different pages can recieve different columns but
    # at the time of writing they all use the same...
    tab_columns = {
        'exp':['Citation', 'Treatment', 'Cell line', 'KO',  'Library', 'DOI', ],
        'comp':['Treatment', 'Dose', 'Timepoint',
                'Growth inhibition %', 'Days grown', 'Cell line', 'KO',
                'Library', 'Citation', 'DOI', ]
    }
    # tab_columns_private = {
    #     'exp':['Treatment', 'Cell line', 'KO',  'Library', 'Citation', 'DOI'],
    #     'comp':['Comparison ID',  'Treatment', 'Dose', 'Timepoint',
    #             'Growth inhibition %', 'Days grown', 'Cell line', 'KO',
    #             'Library', 'Citation', 'DOI']
    # }

    #tab_columns = tab_columns_public if public else tab_columns_private

    if page_id == 'msgv':
        return tab_columns
    elif page_id in ('cm', 'se'):
        return tab_columns
    else:
        raise RuntimeError(f"page_id={page_id} not recognised, only 'cm', 'se' or 'msgv'")



def get_cmdline_options() -> typing.Tuple[str, str, bool]:
    """[script.py] SOURCE PORT [DEBUG]"""
    print('[script.py] SOURCE PORT [DEBUG]')
    import sys
    args = sys.argv
    if (len(args) == 1) or (args[1] in ('-h', '--help')):
        print('usage: comparison_maker.py source_dir port [debug]\n    Any value in the debug position means True.')
        sys.exit(0)
    source = sys.argv[1]
    port = sys.argv[2]
    if len(sys.argv) > 3:
        debug = True
    else:
        debug = False

    return source, port, debug

import dash, pathlib
import dash_bootstrap_components as dbc
def launch_page(source:Union[pathlib.Path, str],
                port:int,
                debug:bool,
                name:str,
                initiate:typing.Callable):
    """Create the app, set debug levels, call `initiate` """
    from dataset import DataSet # not top level to avoid circular import

    app = dash.Dash(name, external_stylesheets=[dbc.themes.BOOTSTRAP])

    if debug:
        logger.level('DEBUG')
    logger.debug(source)
    source_directory = pathlib.Path(source)
    data_set = DataSet(source_directory) #todo update this
    app.layout = initiate(app, data_set, public=True)
    app.run_server(debug=debug, host='0.0.0.0', port=int(port), )


def html_small_span(s):
    """Wrap s in html tags specifying small text using <span>

    Lit: f'<span style="font-size: small;">{s}</span>'
    """
    return f'<span style="font-size: small;">{s}</span>'

def getfuncstr():
    try:
        return inspect.currentframe().f_back.f_code.co_name
    except Exception as e:
        return f'{__name__}.getfuncstr() failed with error:\n\t{e}'


def p_adjust_bh(p:Collection) -> np.ndarray:
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""

    p = np.asfarray(p)
    if any(np.isnan(p)):
        import warnings
        warnings.warn(f'p_adjust_bh(): NaNs in p-values! '
                      f'Called by {inspect.currentframe().f_back.f_code.co_name}')
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def bicolour_cmap(mn, mx, minmax_value=(-3, 3), clr_intensity=170):
    """Return a plotly formatted list of colours that translates to a
    colour map with 0=white, min(mn, -3)=red, max(mx, 3)=blue, so if
    mn > 0 and mx > 0 returned colours will be light-blue -> darker-blue."""
    colours = []
    thresholds = list(minmax_value)
    if mn < -3:
        thresholds[0] = mn
    if mx > 3:
        thresholds[1] = mx

    # ==clr_intensity at threshold, 255 in centre.
    get_primary_int = lambda p: 255 - (p * (255 - clr_intensity))

    for x in (mn, mx):
        if x < 0:
            prop = x / thresholds[0]
            bg_int = int((1 - prop) * 255)
            clr = f"rgb({get_primary_int(prop)}, {bg_int}, {bg_int})"
        else:
            prop = x / thresholds[1]
            print(prop)
            bg_int = int((1 - prop) * 255)
            clr = f"rgb({bg_int}, {bg_int}, {get_primary_int(prop)})"
        colours.append(clr)
    if (mn < 0) and (mx > 0):
        span = mx - mn
        central = (0 - mn) / span
        colours = [(0, colours[0]),
                   (central, 'rgb(255, 255, 255)'),
                   (1, colours[1])]
    else:
        colours = [(i, c) for i, c in enumerate(colours)]
    return colours

def get_treatment_label(row:dict, analysis_label='', inline_style=True) -> typing.Tuple[str, str]:
    """Pass comparison row (either from data_set.comparisons.loc[compid] or
    from dashtable data), return a pair of strings.

    First string comparison specific, second line library, experiment ID."""
    if '-KO' not in row['Treatment']:
        if row['KO'] == 'WT':
            ko = ''
        else:
            ko = f" {row['KO']}-KO"
    else:
        ko = ''

    if analysis_label:
        analysis_label = f"{analysis_label}, "

    if inline_style:
        idstr = f'<span style="font-size: small;">(ID: {row["Comparison ID"]})</span>'
    else:
        idstr = f'(ID: {row["Comparison ID"]})'

    title = (f"Effect of {row['Treatment']} in {row['Cell line']}{ko} cells ({analysis_label}{row['Timepoint']})",
             f"{row['Library']} library {idstr}")

    return title

def get_table_title_text(comp_row, analysis_lab):
    treatment_label = get_treatment_label(
        comp_row,
        analysis_lab,
        inline_style=False,
    )

    return [html.H3(f"{treatment_label[0]}"),
            html.P(f"{treatment_label[1]}")]


def reciprocate_dict(d:dict) -> None:
    """For every key:value, create value:key mapping.

    Checks for conflicts. Mutates in place."""
    # values that are also keys would result in overwriting a key
    #  but {'x':'x'} and {'x':'y', 'y':'x'} are fine.
    v_in_d = [
        v for v in d.values()
            if (v in d.keys())
            and (d[v] != v)
            and (d[v] != d[d[v]])
    ]
    if v_in_d:
        raise RuntimeError(
            f"Can't reciprocate dict as values {v_in_d} "
            f"present in both keys and values of dict."
        )

    for k, v in list(d.items()):
        d[v] = k


def df_rename_columns(df:pd.DataFrame, newcols=dict, inplace=False, axis='columns') -> pd.Index:
    """Return index with renamed columns. Columns missing from newcols will be
    unchanged -- this is the main difference to using df.columns.map(newcols).

    Args:
        df: DataFrame with columns we want to change
        newcols: Mapping of current column names to new ones. Only those we
            want to change need be included.
        inplace: df.columns updated in place. Still returns the new Index.
        axis: "columns"|0 or "index"|1

    Column labels not found as keys in newcols will be retained.

    """

    if axis in ('columns', 0):
        axis = 'columns'
    elif axis in ('index', 1):
        axis = 'index'
    else:
        raise ValueError('axis needs to be "columns"|0 or "index"|1')

    mapper = {k: k for k in getattr(df, axis)}
    mapper.update(newcols)

    nucols = getattr(df, axis).map(mapper)

    if not inplace:
        return nucols

    setattr(df, axis, nucols)
    return nucols

def get_ith_from_all(arr:Sequence[Sequence], index=0):
    """
    Literally: [a[index] for a in arr]
    """
    return [a[index] for a in arr]