from pathlib import Path
import pandas as pd
import os, pickle, typing
from typing import Dict, Collection, TypedDict
from types import MappingProxyType # a kind of frozen dict, though it's slow.

import sqlalchemy
from sqlalchemy import (
    ForeignKey, ForeignKeyConstraint, select, insert, Engine
)
from sqlalchemy import orm
from sqlalchemy.orm import Session, mapped_column, Mapped

from crispr_screen_viewer.functions_etc import (
    timepoint_labels,
    doi_to_link,
    index_of_true,
    LOG,
    get_ith_from_all,
)

from crispr_screen_viewer.database import *

class ScoreFDR(TypedDict):
    score:pd.DataFrame
    fdr:pd.DataFrame

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
    def __init__(self, source_directory, db_engine:Engine=None, print_validations=True):
        LOG.debug(source_directory)
        source_directory = Path(source_directory)
        # for future ref
        self.source_directory = source_directory

        self._use_db = False
        self.engine = None
        if db_engine is not None:
            self._use_db = True
            self.engine = db_engine

        # shorthand internal name: label name
        avail_analyses:list[str] = []
        for ans in ANALYSESTYPES:
             if os.path.isfile(source_directory/f"{ans.shortname}_fdr.csv"):
                 avail_analyses.append(ans.name)
        self.available_analyses = avail_analyses

        # put the data tables in {analysis_type:{score/fdr:pd.DataFrame}} format dictionary
        exp_data = {ans:{stt:pd.read_csv(source_directory/f"{ANALYSESTYPES[ans].shortname}_{stt}.csv", index_col=0)
                         for stt in ('score', 'fdr')}
                    for ans in self.available_analyses}

        # unify the indexes
        genes = pd.Index([])
        # use an index from a table from each analysis
        for analysis in self.available_analyses:
            genes = genes.union(exp_data[analysis]['fdr'].index)
        self.genes = genes

        # reindex with the union of genes
        if not self._use_db:
            self.exp_data = {ans:{stt:exp_data[ans][stt].reindex(genes)
                                for stt in ('score', 'fdr')}
                           for ans in self.available_analyses}
        else:
            self.exp_data = None

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
            pidf = pd.read_csv(
                os.path.join(source_directory, 'previous_and_id.csv'), index_col=0
            )
            self.previous_and_id = pidf.fillna('')

        except FileNotFoundError:
            LOG.warning("file 'previous_and_id.csv' is missing.")
            # when .loc fails to find a name in the table it just uses the current name.
            self.previous_and_id = pd.DataFrame()

        self.previous_and_id.fillna('', inplace=True)

        if print_validations and (not self._use_db):
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


    def comparisons_with_hit(
            self, fdr_max: float,
            genes:Collection[str],
            comparisons:Collection[str],
            analysis_type: str
    ) -> list[str]:
        """List of all comparisons with FDR<fdr_max in given genes."""
        ans_id = ANALYSESTYPES.str_to_id(analysis_type)
        with (Session(self.engine) as session):
            comp_ids = session.query(StatTable.comparison_id).where(
                StatTable.fdr < fdr_max,
                StatTable.gene_id.in_(genes),
                StatTable.comparison_id.in_(comparisons),
                StatTable.analysis_type_id == ans_id,
            ).distinct().all()

        return get_ith_from_all(comp_ids, 0)


    def get_score_fdr(
            self,
            score_anls:str,
            fdr_anls:str=None,
            comparisons:Collection[str]='ALL',
            genes:Collection[str]='ALL',
            #data_sources:Collection= 'NOT IMPLEMENTED'
    ) -> Dict[str, pd.DataFrame]:
        """Get score and FDR tables for the analysis types & data sets.
        Tables give the per gene values for included comparisons.

        Arguments:
            score_anls: The analysis type from which to get the score values per gene
            fdr_anls: Optional. Default = score_anls. If set to False, get_score_fdr['fdr'] is None.
            comparisons: list of comparisons to include
            genes: list of genes to include

        Returns {'score':pd.DataFrame, 'fdr':pd.DataFrame}"""
        LOG.debug(f"getting score fdr, n genes={len(genes)}, n comps={len(comparisons)}, fdr={fdr_anls}, score={score_anls}")
        print(genes, comparisons)

        def run_query(analysis_type: str) -> list:
            ans_id = ANALYSESTYPES.str_to_id(analysis_type)
            with (Session(self.engine) as session):
                constraints = [StatTable.analysis_type_id == ans_id]
                if comparisons != 'ALL':
                    constraints.append(
                        StatTable.comparison_id.in_(comparisons)
                    )

                if genes != 'ALL':
                    constraints.append(
                        StatTable.gene_id.in_(genes)
                    )

                # if fdr_max is not None:
                #     constraints.append(
                #         StatTable.fdr < fdr_max
                #     )

                query = session.query(
                    StatTable.gene_id, StatTable.comparison_id,
                    StatTable.score, StatTable.fdr
                ).where(
                    *constraints
                )

                results = query.all()

            return results

        def pivot_results(results, stat: str):

            data = [(r.gene_id, r.comparison_id, getattr(r, stat)) for r in results]

            # Convert to DataFrame
            df = pd.DataFrame(data, columns=['gene_id', 'comparison_id', 'score'])

            # Pivot the DataFrame to get genes as index and comparisons as columns
            pivot_df = df.pivot(index='gene_id', columns='comparison_id', values='score')

            return pivot_df

        query_res = run_query(score_anls)
        scores = pivot_results(query_res, 'score')

        if fdr_anls is False:
            return ScoreFDR(score=scores, fdr=None)

        if fdr_anls is None:
            fdrs = pivot_results(query_res, 'fdr')
        else:
            query_res = run_query(fdr_anls)
            fdrs = pivot_results(query_res, 'fdr')

        return ScoreFDR(score=scores, fdr=fdrs)


    def dropdown_gene_label(self, gn:str) -> str:
        """Gene symbol plus HGNC ID and previous symbols, if available."""
        try:
            row = self.previous_and_id.loc[gn]
        except:
            return gn

        if not row.hgnc_id:
            return gn

        s = f"{gn}  ({row.hgnc_id}"
        s_end = ')'
        if row.prev_symbol:
            s_end = f"; {row.prev_symbol})"

        return s + s_end

def load_dataset(paff, db_engine:Engine=None):
    """If paff is a dir, the dataset is constructed from the files
    within, otherwise it is assumed to be a pickle."""
    if os.path.isfile(paff):
        LOG.info('args.data_path is a file, assuming pickle and loading.')
        with open(paff, 'rb') as f:
            data_set = pickle.load(f)
    else:
        data_set = DataSet(Path(paff), db_engine=db_engine)

    return data_set

import dataclasses
@dataclasses.dataclass(frozen=True)
class AnalysisType:
    id: int
    name: str
    shortname: str
    label: str
    score_label: str


class _AnalysesTypes:
    """Container class with available analyses. Access analyses with
    name or ID (starts at 1), or iterate through.

    Use the singleton ANALYSESTYPES"""
    def __init__(self, analyses_types:Collection[AnalysisType], default_type='drugz'):

        by_str = {v.name: v for v in analyses_types}
        by_str.update({v.shortname:v for v in analyses_types})
        self.by_str = MappingProxyType(by_str)
        self.by_id = MappingProxyType({v.id:v for v in analyses_types})
        self.list = tuple(analyses_types)
        self.default = self[default_type]

        rev_enumerator = sorted(enumerate(analyses_types), key=lambda x: x[0], reverse=True)

        self.binary_values = {s.name: 2 ** i for i, s in rev_enumerator}

    def encode_bitmask(self, analyses_names: list[str]):
        return sum([self.binary_values[an] for an in analyses_names])

    def decode_bitmask(self, val:int) -> set[str]:
        present = set()
        for s, i in self.binary_values.items():
            if val >= i:
                present.add(s)
                val -= i
        return present

    def str_to_id(self, name):
        return self.by_str[name].id

    def __getitem__(self, item:typing.Union[str, int]) -> AnalysisType:
        if isinstance(item, str):
            return self.by_str[item]
        elif isinstance(item, int):
            return self.by_id[item]
        else:
            raise ValueError(f"Access analyses by str name or int id not {type(item)}")

    def __iter__(self):
        return iter(self.list)

ANALYSESTYPES = _AnalysesTypes([
    AnalysisType(
        id=1, name='mageck', shortname='mag', label='MAGeCK', score_label='Log2(FC)'
    ),
    AnalysisType(
        id=2, name='drugz', shortname='drz', label='DrugZ',  score_label='NormZ'
    )
])

def comps_with_analysis_type(
        ans_name_id: typing.Union[str, int],
        engine: Engine
) -> list[int]:
    """List comparison numeric IDs that have results for given analysis type.

    Args:
        ans_name_id: string or numeric ID of analysis type
        engine: SQLAlchemy engine."""
    name = ANALYSESTYPES[ans_name_id].name

    b = ANALYSESTYPES.binary_values[name]

    with Session(engine) as S:
        comps = S.execute(
            select(ComparisonTable.stringid)
            .where(
                ComparisonTable.analyses_bitmask.bitwise_and(b) == b
            )
        ).all()

    return [c[0] for c in comps]


def sample_dataset(inpath,
                    outpath,
                    n_exps=10,
                    max_comps=20,
                    n_genes=1000):
    """Create a smaller dataset sampled from a larger one"""
    short_data = {}

    p = Path(inpath)

    data: dict[str, pd.DataFrame] = {fn.replace('.csv', ''): pd.read_csv(p / fn, index_col=0) for fn in os.listdir(p)}

    stat_tbl_keys = ['mag_pos_p',
                     'mag_neg_p',
                     'mag_fdr',
                     'mag_score',
                     'drz_pos_p',
                     'drz_neg_p',
                     'drz_fdr',
                     'drz_score']

    # shorten the data
    seed = 810613

    exp_ids = data['experiments_metadata'].sample(n_exps, replace=False, random_state=seed).index
    comps = data['comparisons_metadata'].loc[data['comparisons_metadata']['Experiment ID'].isin(exp_ids)]
    comps = comps.sample(max_comps, replace=False, random_state=seed)
    exps = data['experiments_metadata'].loc[comps['Experiment ID'].unique()]

    short_data['experiments_metadata'] = exps
    short_data['comparisons_metadata'] = comps

    for k in stat_tbl_keys:
        tbl = data[k]
        tbl = tbl.loc[:, comps.index]
        tbl = tbl.loc[~tbl.isna().all(1)]
        tbl = tbl.head(n_genes)
        short_data[k] = tbl

    print(tbl.shape)
    os.makedirs(outpath, exist_ok=True)
    for p, d in short_data.items():
        d.to_csv(outpath+'/'+p)

    return short_data