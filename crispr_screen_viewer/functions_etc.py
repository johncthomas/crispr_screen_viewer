import typing

from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
import os
from pathlib import Path
from scipy.stats import norm
from scipy import odr
from scipy.stats import linregress
from typing import Union, List, Dict, Iterable, Collection
import logging, pickle
parse_expid = lambda comp: comp.split('.')[0]

logging.basicConfig()
LOG = logging.getLogger('screen_viewers')

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
    """Return string formated as a markdown link to doi.org/{doi}"""
    # should be formated to not be a hyperlink, but sometimes it is
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

class DataSet:
    """Class for holding, retrieving screen data and metadata.

    Storage and retreval of amalgamated results tables. Experiment tables as
    single analysis type/statistic with filename {ans_type}_{stat}.csv
    Currently supports MAGeCK ('mag') and DrugZ ('drz').

    Attributes:
        exp_data: data tables keyed first by analysis type then 'score'|'fdr'
        comparisons: descriptions of exp_data samples. Indexed by comparison keys
        score/analysis_labels: text labels for supported analysis types
        genes: Index used in tables in exp_data

    Methods:
        get_results: get dict of score and fdr for specified analysis types and
            datasets."""
    def __init__(self, source_directory, print_validations=True):
        source_directory = Path(source_directory)

        # shorthand internal name: label name
        avail_analyses = []
        for ans in ('drz', 'mag'):
             if os.path.isfile(source_directory/f"{ans}_fdr.csv"):
                 avail_analyses.append(ans)

        self.available_analyses = avail_analyses
        self.analysis_labels = {'drz':'DrugZ', 'mag':'MAGeCK'}
        self.score_labels = {'mag':'Log2(FC)', 'drz':'NormZ'}

        # put the data tables in {analysis_type:{score/fdr:pd.DataFrame}} format dictionary
        exp_data = {ans:{stt:pd.read_csv(source_directory/f"{ans}_{stt}.csv", index_col=0)
                         for stt in ('score', 'fdr')}
                    for ans in self.available_analyses}

        # unify the indexes
        genes = pd.Index([])
        # use an index from a table from each analysis
        for analysis in self.available_analyses:
            genes = genes.union(exp_data[analysis]['fdr'].index)
        self.genes = genes

        # reindex with the union of genes
        self.exp_data = {ans:{stt:exp_data[ans][stt].reindex(genes)
                            for stt in ('score', 'fdr')}
                       for ans in self.available_analyses}

        comparisons = pd.read_csv(f'{source_directory}/comparisons_metadata.csv', )
        # this is sometimes put in wrong...
        m = comparisons['Timepoint'] == 'endpoint'
        comparisons.loc[m, 'Timepoint'] = 'endpoints'
        comparisons.loc[m, 'Control group'] = comparisons.loc[m, 'Control group'] \
            .apply(lambda x: x.replace('endpoint', 'endpoints'))

        comparisons = comparisons.set_index('Comparison ID', drop=False)
        #comparisons.loc[:, 'Available analyses'] = comparisons['Available analyses'].str.split('|')
        try:
            comparisons = comparisons.drop('Available analyses', axis=1)
        except KeyError:
            pass
        comparisons.loc[comparisons.Treatment.isna(), 'Treatment'] = 'No treatment'
        # these cols could be blank and aren't essential to have values
        for col in  ['Cell line', 'Library', 'Source']:
            if col not in comparisons.columns:
                comparisons.loc[:, col] = 'Unspecified'
            comparisons.loc[comparisons[col].isna(), col] = 'Unspecified'

        # replace some values with ones that read better
        for old, new in timepoint_labels.items():
            comparisons.loc[comparisons.Timepoint == old, 'Timepoint'] = new

        # list of all datasources for filtering
        self.data_sources = comparisons.Source.fillna('Unspecified').unique()
        # main metadata tables
        self.comparisons = comparisons
        self.experiments_metadata = pd.read_csv(f'{source_directory}/experiments_metadata.csv', )
        # rename "Experiment name" to "Experiment ID" for consistency
        colmap = {k:k for k in self.experiments_metadata}
        colmap['Experiment name'] = 'Experiment ID'
        self.experiments_metadata.columns = self.experiments_metadata.columns.map(colmap)
        self.experiments_metadata.set_index('Experiment ID', drop=False, inplace=True)


        # add formated DOI to the comparisons metadata
        dois = self.experiments_metadata.loc[
            self.comparisons['Experiment ID'],
            'DOI'
        ].apply(doi_to_link).values

        self.comparisons.insert(2, 'DOI', dois)

        if print_validations:
            # Print information that might be helpful in spotting data validity issues
            # Check for comparisons present in the metadata/actual-data but missing
            #   in the other
            all_good = True
            for ans in self.available_analyses:
                score_comps = self.exp_data[ans]['score'].columns
                meta_comps = self.comparisons.index

                meta_in_score = meta_comps.isin(score_comps)
                missing_in_data = meta_comps[~meta_in_score]
                # todo log.warning
                # todo check experiments metadata
                if missing_in_data.shape[0] > 0:
                    all_good = False
                    print(
                        f"Comparisons in comparisons metadata but not in {ans}_score.csv:"
                        f"\n    {', '.join(missing_in_data)}\n"
                    )
                score_in_meta = score_comps.isin(meta_comps)
                missing_in_score = score_comps[~score_in_meta]
                if missing_in_data.shape[0] > 0:
                    all_good = False
                    print(
                        f"Comparisons in {ans}_score.csv, but not in comparisons metadata:"
                        f"\n    {', '.join(missing_in_score)}\n"
                    )
            comps = self.comparisons.index
            if comps.duplicated().any():
                all_good = False
                print('Duplicate comparisons found - this will probably stop the server from working:')
                print('   ' , ', '.join(sorted(comps.index[comps.index.duplicated(keep=False)])))
            if all_good:
                print(f'All comparisons data in {source_directory} are consistent')

            # check all comparison ExpID appear in experiments metadata
            # the reverse isn't fatal
            expids = self.comparisons['Experiment ID'].unique()
            found = [xi in self.experiments_metadata.index for xi in expids]
            if not all(found):
                not_found = [x for (x, b) in zip(expids, found) if not b]
                print('Experiment IDs used in comparisons_metadata not found in experiments_metadata:\n'
                      f'   {", ".join(not_found)}')



    def get_score_fdr(self, score_anls:str, fdr_anls:str=None,
                      data_sources:Collection= 'all') -> Dict[str, pd.DataFrame]:
        """Get score and FDR tables for the analysis types & data sets.
        Tables give the per gene values for included comparisons.

        Arguments:
            score_anls: The analysis type from which to get the score values
                per gene
            fdr_anls: Optional. As score_anslys
            data_sources: Data sources (i.e. SPJ, or other peoples papers) to
                include in the returned DFs. Any comparison that comes from a
                dataset that does not have both fdr/score analysis types
                available will not be present in the table.

        Returns {'score':pd.DataFrame, 'fdr':pd.DataFrame}"""

        # if only one type supplied, copy it across
        if fdr_anls is None:
            fdr_anls = score_anls

        score_fdr = {stt:self.exp_data[ans][stt] for ans, stt in ((score_anls, 'score'), (fdr_anls, 'fdr'))}

        if data_sources == 'all':
            return score_fdr

        # Filter returned comparisons (columns) by inclusion in data sources and having
        #   results for both analysis types
        comps_mask = self.comparisons.Source.isin(data_sources)
        # for analysis_type in (score_anls, fdr_anls):
        #     m = self.comparisons['Available analyses'].apply(lambda available: analysis_type in available)
        #     comps_mask = comps_mask & m
        comparisons = index_of_true(comps_mask)
        score_fdr = {k:tab.reindex(columns=comparisons) for k, tab in score_fdr.items()}

        return score_fdr


def load_dataset(paff):
    """If paff is a dir, the dataset is constructed from the files
    within, otherwise it is assumed to be a pickle."""
    if os.path.isfile(paff):
        LOG.info('args.data_path is a file, assuming pickle and loading.')
        with open(paff, 'rb') as f:
            data_set = pickle.load(f)
    else:
        data_set = DataSet(Path(paff))
    return data_set

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
                'Library', 'Experiment ID',]
    }
    if not public:
        for k, l in filter_keys.items():
            l.append('Source')
    return filter_keys


def get_metadata_table_columns(public, page_id) -> Dict[str, List[str]]:
    # this is set up so that different pages can recieve different columns but
    # at the time of writing they all use the same...
    tab_columns_public = {
        'exp':[ 'Citation', 'Treatment', 'Cell line', 'KO',  'Library', 'DOI', 'Experiment ID', ],
        'comp':['Comparison ID',  'Treatment', 'Dose', 'Timepoint',
                'Growth inhibition %', 'Days grown', 'Cell line', 'KO',
                'Library', 'Experiment ID', 'DOI']
    }
    tab_columns_private = {
        'exp':['Treatment', 'Cell line', 'KO',  'Library', 'Citation', 'Experiment ID', 'DOI'],
        'comp':['Comparison ID',  'Treatment', 'Dose', 'Timepoint',
                'Growth inhibition %', 'Days grown', 'Cell line', 'KO',
                'Library', 'Experiment ID', 'DOI']
    }

    tab_columns = tab_columns_public if public else tab_columns_private

    if page_id == 'msgv':
        return tab_columns
    elif page_id in ('cm', 'se'):
        return tab_columns
    else:
        raise RuntimeError(f"page_id={page_id} not recognised, only 'cm', 'se' or 'msgv")



def get_cmdline_options() -> typing.Tuple[str, str, bool]:
    """[script.py] SOURCE PORT [DEBUG]"""
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