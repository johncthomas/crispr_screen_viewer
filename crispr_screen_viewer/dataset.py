from pathlib import Path
import pandas as pd
import os, pickle
from typing import Dict, Collection
from crispr_screen_viewer.functions_etc import (
    timepoint_labels,
    doi_to_link,
    index_of_true,
    LOG,
)

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
        # for future ref
        self.source_directory = source_directory

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

        comparisons = pd.read_csv(source_directory/'comparisons_metadata.csv', )
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

        # fill in blank treatments
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

        # add citation to comparisons table
        try:
            cites = self.experiments_metadata.loc[
                self.comparisons['Experiment ID'],
                'Citation'
            ].values
            self.comparisons.loc[:, 'Citation'] =  cites
        except:
            LOG.warning('Citations column missing from exeriments_metadata')
            self.comparisons.loc[:, 'Citation'] = ''


        # DF of previous symbols and IDs for currently used.
        try:
            self.previous_and_id = pd.read_csv(
                os.path.join(source_directory, 'previous_and_id.csv'), index_col=0
            )
            self.previous_and_id.fillna('', inplace=True)

        except FileNotFoundError:
            LOG.warning("file 'previous_and_id.csv' is missing.")
            # when .loc fails to find a name in the table it just uses the current name.
            self.previous_and_id = pd.DataFrame()

        self.previous_and_id.fillna('', inplace=True)

        if print_validations:
            self.validate_comparisons()
            self.validate_previous_and_id()

    def validate_comparisons(self):
        """Print information that might be helpful in spotting data validity issues
        Check for comparisons present in the metadata/actual-data but missing
          in the other"""
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
            print(f'All comparisons data in {self.source_directory} are consistent')

        # check all comparison ExpID appear in experiments metadata
        # the reverse isn't fatal
        expids = self.comparisons['Experiment ID'].unique()
        found = [xi in self.experiments_metadata.index for xi in expids]
        if not all(found):
            not_found = [x for (x, b) in zip(expids, found) if not b]
            print('Experiment IDs used in comparisons_metadata not found in experiments_metadata:\n'
                  f'   {", ".join(not_found)}')

    def validate_previous_and_id(self):
        """Check all stats index in datasets are in previous_and_id"""
        for ans in self.available_analyses:
            score_index = self.exp_data[ans]['score'].index
            m = score_index.isin(self.previous_and_id.index)
            print(m.sum(), 'of', len(m), f'gene symbols have record in previous_and_id.csv, in file {ans}_score.csv')


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


        # todo surely we don't need to do this every time?
        #   write tables for each score/stat, store in a dict

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

    def dropdown_gene_label(self, gn):
        try:
            row = self.previous_and_id.loc[gn]
        except:
            return gn

        if not row.HGNC_ID:
            return gn

        s = f"{row.Symbol}  ({row.HGNC_ID}"
        s_end = ')'
        if row.Previous_symbol:
            s_end = f"; {row.Previous_symbol})"

        return s + s_end

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