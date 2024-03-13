import copy
import inspect
import typing
import unicodedata
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



def load_mageck_tables(prefix:str, controls:Iterable[str]):
    """Get a dict of DF of mageck results keyed to control groups.
    Returned sample index only gives the treatment sample, not ctrl-treat,
    except for the EXTRA group, which needs both"""
    tables = {}
    for ctrl_group in controls:
        tab = tables[ctrl_group] = tabulate_mageck(prefix+ctrl_group+'.')
        if ctrl_group != 'EXTRA':
            tab.columns.set_levels([c.split('-')[1] for c in tab.columns.levels[0]], 0, inplace=True)
    return tables

def get_annotation_dicts(df, xkey, ykey, annote_kw=None):
    """A list of dict defining annotations from specified columns
    of a dataframe."""
    annotations = []
    if annote_kw is None:
        annote_kw = {}
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


def tabulate_mageck(prefix):
    """
    :param prefix: Input file prefix, including path
    :return: pd.DataFrame
    """
    prefix = Path(prefix)
    tables = {}
    tab = None
    for fn in os.listdir(prefix.parent):
        if not fn.endswith('.gene_summary.txt') or \
                prefix.parts[-1] not in fn:
            continue
        #mtab from the mageck output, reformatted into tab
        mtab = pd.read_csv(prefix.parent / fn, '\t', index_col=0)
        tab = pd.DataFrame(index=mtab.index)
        tab.loc[:, 'lfc'] = mtab.loc[:, 'neg|lfc']
        # turn sep pos|neg columns into one giving only the appropriate LFC/FDR
        pos = mtab['pos|lfc'] > 0
        for stat in 'fdr', 'p-value':
            tabk = stat.replace('-value', '')
            tab.loc[pos, tabk] = mtab.loc[pos, f'pos|{stat}']
            tab.loc[~pos, tabk] = mtab.loc[~pos, f'neg|{stat}']
            tab.loc[:, f'{tabk}_log10'] = tab[tabk].apply(lambda x: -np.log10(x))

        sampnm = fn.split(prefix.stem)[1].split('.gene_s')[0]
        tables[sampnm] = tab
    if tab is None:
        raise FileNotFoundError('Failed to find any .gene_summary.txt files with prefix '+str(prefix))

    # create multiindex using the sample names and stats
    tbcolumns = pd.MultiIndex.from_product(
        [sorted(tables.keys()), ['lfc', 'fdr', 'fdr_log10', 'p', 'p_log10']],
        1
    )
    table = pd.DataFrame(index=tab.index, columns=tbcolumns)
    for exp, tab in tables.items():
        table[exp] = tab
    return table

def tabulate_jacks(prefix):
    """Return 3 tables giving the:
        1. jacks_score|fdr_pos|fdr_neg|std,
        2. guide efficacy data,
    and 3. fold changes for each gene.

    Tables are multiindexed by sample name and then results columns for those
    samples.

    fdr in the

    Prefix is the used to identify the results files. So prefix
    should contain the path to the files if they aren't in os.getcwd()"""

    kwtab = dict(sep='\t', index_col=0)

    sig_df = tabulate_score(prefix)
    samples = sig_df.columns.levels[0]
    # get guide data, foldchange and efficacies
    guide_cols = pd.MultiIndex.from_product((samples, ['foldchange', 'fold_std', 'eff', 'eff_std']),
                                            names=['exp', 'stat'])
    fchange_df = pd.DataFrame(columns=guide_cols)
    foldchange = pd.read_csv(prefix + '_logfoldchange_means.txt', **kwtab)
    foldstd = pd.read_csv(prefix + '_logfoldchange_std.txt', **kwtab)
    eff_tab = pd.read_csv(prefix + '_grna_JACKS_results.txt', **kwtab)

    for exp in samples:
        fchange_df.loc[:, (exp, 'lfc')] = foldchange[exp]
        fchange_df.loc[:, (exp, 'fold_std')] = foldstd[exp]
    fchange_df.loc[:, 'gene'] = foldchange['gene']

    efficacies = pd.DataFrame(columns=('eff', 'eff_std'))
    efficacies.loc[:, 'eff'] = eff_tab['X1']
    efficacies.loc[:, 'eff_std'] = (eff_tab['X2'] - eff_tab['X1'] ** 2) ** 0.5
    efficacies.loc[:, 'gene'] = fchange_df['gene']

    return sig_df, efficacies, fchange_df


def tabulate_score(prefix, return_ps=False):
    """Return a multiindexed DF of JACKS results.

    Table columns are sample names as given in the repmap at level 0,
    and then 'jacks_score, fdr_log10, fdr_neg, fdr_pos, stdev' at level 1"""
    # othertab = pd.DataFrame(columns=("IC10","IC90","D14"), index=essen['D14'].index)
    # Tables produced by jacks have columns that are the groups
    genes = pd.read_csv(prefix + '_gene_JACKS_results.txt', sep='\t', index_col=0)
    genes_index = sorted(genes.index)
    genes = genes.reindex(genes_index)
    genesstd = pd.read_csv(prefix + '_gene_std_JACKS_results.txt', sep='\t', index_col=0)
    genesstd = genesstd.reindex(genes_index)
    ps = genes / genesstd
    ps = ps.apply(norm.cdf)

    # multiindex DF for each experiment giving results
    sig_cols = pd.MultiIndex.from_product((ps.columns, ['jacks_score', 'fdr_pos', 'fdr_neg', 'fdr_log10', 'stdev']),
                                          names=('exp', 'stat'))
    sig_df = pd.DataFrame(index=genes_index, columns=sig_cols)

    for exp in ps.columns:
        sig_df.loc[:, (exp, 'p_neg')] = ps[exp]
        sig_df.loc[:, (exp, 'p_pos')] = 1-ps[exp]
        # sig_df.loc[:, (exp, 'fdr_neg')] = multipletests(ps[exp], method='fdr_bh')[1]
        # sig_df.loc[:, (exp, 'fdr_pos')] = multipletests(1 - ps[exp], method='fdr_bh')[1]
        sig_df.loc[:, (exp, 'fdr_neg')] = p_adjust_bh(ps[exp])
        sig_df.loc[:, (exp, 'fdr_pos')] = p_adjust_bh(1-ps[exp])
        sig_df.loc[:, (exp, 'jacks_score')] = genes[exp]
        sig_df.loc[:, (exp, 'stdev')] = genesstd[exp]
        score_table = sig_df[exp].copy()
        # get one FDR
        pos = score_table['jacks_score'] > 0
        neg = ~pos
        # get lowest not zero and set zeroes to 1/10 that value
        min_pos = min(score_table.loc[score_table['fdr_pos'] > 0, 'fdr_pos'])
        min_neg = min(score_table.loc[score_table['fdr_neg'] > 0, 'fdr_neg'])

        for fdri, fdr in enumerate(['fdr_pos', 'fdr_neg']):
            score_table.loc[score_table[fdr] == 0, fdr] = (min_pos, min_neg)[fdri] / 10

        for mask, posneg in (pos,'pos'), (neg, 'neg'):
            score_table.loc[mask, 'fdr'] = score_table.loc[mask, 'fdr_'+posneg]
            score_table.loc[mask, 'p'] = score_table.loc[mask, 'p_' + posneg]

        sig_df.loc[:, (exp, 'fdr_log10')] = score_table['fdr'].apply(lambda x: -np.log10(x))
        sig_df.loc[:, (exp, 'p_log10')] = score_table['p'].apply(lambda x: -np.log10(x))
        # put the columns in a good order
        if not return_ps:
            sig_df = sig_df.reindex(['jacks_score', 'fdr_log10', 'stdev'], axis=1, level=1, )

    return sig_df

def iter_files_by_prefix(prefix:Union[str, Path], req_suffix=None):
    prefix = Path(prefix)
    check_suffix = lambda s: s.endswith(req_suffix) if req_suffix is not None else True

    for fn in os.listdir(prefix.parent):
        # filter incorrect files, the '._' are mac files that are not ignored automatically on unix filesystem
        if not check_suffix(fn) or \
                prefix.parts[-1] not in fn or \
                fn.startswith('._'):
            continue
        yield fn

def tabulate_drugz_files(file_names, prefix, compjoiner='-'):
    prefix=Path(prefix)
    tables = {}
    for fn in file_names:
        fn_prefix = prefix.name
        # remove fn parent (directories), remove file prefix, remove file suffix
        #   what's left is the comp
        comp = Path(fn).name.replace(fn_prefix, '')[:-4]
        comp = comp.replace('-', compjoiner)

        tab = pd.read_csv(fn, sep='\t', index_col=0)

        # sort out the column names to be consistent with other results tables
        tab.index.name = 'gene'
        tab = tab.loc[:, ['normZ', 'pval_synth', 'fdr_synth', 'pval_supp', 'fdr_supp']]
        stats_cols = 'normZ neg_p neg_fdr pos_p pos_fdr'.split()

        tab.columns = stats_cols

        # get the minimum value significance stats
        for stat in 'p', 'fdr':
            min_stat = tab.loc[:, [f'neg_{stat}', f'pos_{stat}']].min(1)
            tab.loc[:, f'{stat}'] = min_stat
            tab.loc[:, f'{stat}_log10'] = min_stat.apply(lambda p: -np.log10(p))
        tables[comp] = tab

    tbcolumns = pd.MultiIndex.from_product(
        [sorted(tables.keys()), tab.columns],
        1
    )
    table = pd.DataFrame(index=tab.index, columns=tbcolumns)
    for exp, tab in tables.items():
        table[exp] = tab
    return table


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