from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
import os
from pathlib import Path
from scipy.stats import norm
from typing import Iterable

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
        sig_df.loc[:, (exp, 'fdr_neg')] = multipletests(ps[exp], method='fdr_bh')[1]
        sig_df.loc[:, (exp, 'fdr_pos')] = multipletests(1 - ps[exp], method='fdr_bh')[1]
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