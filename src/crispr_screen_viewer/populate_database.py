#!/usr/bin/env python
import dataclasses
import logging
import os.path, typing
import pathlib
from pathlib import Path

import pandas as pd
from sqlalchemy import Engine, create_engine, select
from sqlalchemy.orm import Session

from database import *

from crispr_screen_viewer import functions_etc
from crispr_screen_viewer.functions_etc import df_rename_columns
from crispr_screen_viewer.dataset import ANALYSESTYPES, AnalysisType

from crispr_tools.data_classes import AnalysisWorkbook, CrisprCounts
from crispr_tools.crispr_pipeline import analysis_tabulate
analysis_tabulate:dict[str, typing.Callable]

import numpy as np

def create_engine_with_schema(destination="sqlite://", echo=False) -> Engine:
    engine = create_engine(destination, echo=echo, connect_args={'timeout':4000})

    # create tables in the actual SQL database
    Base.metadata.create_all(engine)

    return engine


def prep_deets_columns(tbl:pd.DataFrame, cmap:dict):
    """drop not defined columns, rename columns.

    args:
        tbl: a comparison or experiment details dataframe
        cmap: mapper for renaming columns
    """
    tbl = tbl[cmap.keys()]
    cols = functions_etc.df_rename_columns(tbl, cmap)
    tbl.columns = cols
    return tbl


def insert_df_rows(tablecls, df:pd.DataFrame, session:Session):
    """Take row values from df, and add them as rows to the table.

    NaN are skipped over.

    Assumes everything is formated correctly with df.columns matching field names.
    Fails if they've already been added, dunno how you update them

    Does not commit"""
    # https://docs.sqlalchemy.org/en/20/tutorial/data_update.html

    # filter the records so key:value pairs of null value are removed
    session.add_all(
        [tablecls(**{k: v for k, v in kw.items() if not pd.isna(v)})
         for kw in df.to_dict(orient='records')]
    )



# def format_experiments_df(
#         experiments_metadata:pd.DataFrame,
# ) -> pd.DataFrame:
#     """Reformat experiments DF for insertion into database"""
#     # experiment
#     experiment_column_mapping = {
#         'Experiment ID': 'stringid',
#         'Date published': 'date',
#         'Library': 'library',
#         'Multiplicity of infection': 'moi',
#         'Representation average': 'representation',
#         'Experiment description (a few sentences)': 'description',
#         'Notes': 'notes',
#         'DOI': 'doi',
#         'Reference': 'reference',
#         'Source': 'source',
#         'Citation': 'citation',
#     }
#
#     exp = experiments_metadata.reset_index()
#     exp = prep_deets_columns(exp, experiment_column_mapping)
#     d = exp.date.str.replace('/', '-')
#     exp.date = pd.to_datetime(d, format='mixed', dayfirst=True)
#
#     return exp


def get_exp_ids(engine:Engine) -> dict[str, int]:
    with Session(engine) as S:
        rows = S.execute(select(ExperimentTable.stringid, ExperimentTable.id)).all()

    return dict(rows) # noqa

def get_comp_ids(engine:Engine) -> dict[str, int]:
    with Session(engine) as S:
        rows = S.execute(select(ComparisonTable.stringid, ComparisonTable.id)).all()

    return dict(rows) # noqa


def format_comparisons_df(
        comparisons_metadata:pd.DataFrame,
        experiment_ids:dict[str, int]=None,
) -> pd.DataFrame:

    comparison_column_mapping = {
        'Comparison ID': 'stringid',
        'Treatment': 'treatment_label',
        'Timepoint': 'timepoint',
        'Cell line': 'cell',
        'Ctrl samp': 'ctrl',
        'Treat samp': 'treat',
        'KO': 'ko',
        'Dose': 'dose',
        'Growth inhibition %': 'gi',
        'Days grown': 'days_grown',
        'Library': 'library',
        'analyses_bitmask': 'analyses_bitmask',
        'Experiment ID': 'experiment'
    }

    comps = comparisons_metadata.reset_index()
    #comps.index.name = 'id'
    comps.head()

    comps.loc[:, 'analyses_bitmask'] = comps['Available analyses'].apply(
        lambda s: ANALYSESTYPES.encode_bitmask([ss for ss in s.split('|')])
    )

    comps = prep_deets_columns(comps, comparison_column_mapping)
    # comps = comps[comparison_column_mapping.keys()]
    # cols = functions_etc.df_rename_columns(comps, comparison_column_mapping, )
    # comps.columns = cols

    if experiment_ids is not None:
        comps.experiment = comps.experiment.map(experiment_ids)

    return comps


def format_statistics_dfs(
        statistics:dict[str, pd.DataFrame],
        comparison_ids:dict[str, int]=None,
        analyses=('mageck', 'drugz'),
        stats_names=('pos_p', 'neg_p', 'fdr', 'score')
) -> list[pd.DataFrame]:
    """Convert wide tables in dict keyed by things like "mag_pos_p" or "drz_score"
    to list of long tables for insertion into DB.
    """

    stats_dfs = []

    for method_name in analyses:
        method = ANALYSESTYPES[method_name]

        stt_values = {}
        for stt in stats_names:
            d = statistics[f"{method.shortname}_{stt}"]

            # ...NaN in the index can cause crashes. Why can there be a NaN in the index?
            #   found in a few mageck results tables.
            d = d.loc[~d.index.isna()]
            stt_values[stt] = d.unstack()

        df = pd.DataFrame(stt_values).reset_index()

        functions_etc.df_rename_columns(
            df,
            {'Comparison ID': 'comparison_id', 'id':'gene_id', 'gene':'gene_id', 'level_0':'comparison_id'},
            inplace=True
        )

        df.loc[:, 'fdr10'] = df.fdr.apply(lambda p: -np.log10(p))

        if comparison_ids is not None:
            df.comparison_id = df.comparison_id.map(comparison_ids)

        df.loc[:, 'analysis_type_id'] = method.id

        stats_dfs.append(df)

    return stats_dfs


# def create_database(
#         experiments_metadata:pd.DataFrame,
#         comparisons_metadata:pd.DataFrame,
#         stats_tables:dict[str, pd.DataFrame],
#         db_url:str,
#         use_autogen_ids=False,
#         skip_existing_experiments=True,
# ):
#     """Create a SQL database experimental results from dataframes.
#
#     `stats` organised keyed by experiment method short name ('mag', 'drz') then
#     by stat type ('pos_p', 'fdr', etc)."""
#     if use_autogen_ids:
#         raise NotImplementedError("Currently using string IDs from experiment/sample names, "
#                                   "using autogenerated IDs will require the schema to be changed.")
#     exp = format_experiments_df(experiments_metadata)
#
#     engine = create_engine_with_schema(
#         db_url
#     )
#
#     # needs to be inserted first to generate IDs
#     insert_df_rows(ExperimentTable, exp, engine)
#
#
#     exp_ids = get_exp_ids(engine) if use_autogen_ids else None
#     comps = format_comparisons_df(comparisons_metadata, exp_ids)
#     insert_df_rows(ComparisonTable, comps, engine)
#
#     comp_ids = get_comp_ids(engine) if use_autogen_ids else None
#     stats = format_statistics_dfs(stats_tables, comp_ids)
#     for tbl in stats:
#         insert_df_rows(StatTable, tbl, engine)


def split_data(data:dict[str, pd.DataFrame]):
    out_data = {}
    for k in ('experiments_metadata', 'comparisons_metadata'):
        out_data[k] = data[k]
        del data[k]

    out_data['stats_tables'] = data

    return out_data

def __create_test_db():
    os.chdir('/mnt/m/stuff/ddrcs/')

    dbpaff = '/home/jthomas/tmp/test_db1.medium.db'
    dspaff = 'dataset/2023-02-23_external/'
    smlpaff = '/home/jthomas/tmp/meddataset'
    if os.path.exists(dbpaff):
        os.remove(dbpaff)
    from dataset import sample_dataset
    short_data = sample_dataset(dspaff, smlpaff, 20, 30, 5000)

    create_database(**split_data(short_data), db_url=f"sqlite:///{dbpaff}", use_autogen_ids=False)

def load_csv(paff:str):
    """Load all .csv in paff, should be experiments_metadata.csv, comparisons_metadata.csv and
    individual tables for the statistics. The stats tables get put under "stats_tables" key. """
    data: dict[str, pd.DataFrame] = {
        fn.replace('.csv', ''): pd.read_csv(pathlib.Path(paff) / fn, index_col=0)
        for fn in os.listdir(paff) if fn.endswith('.csv')}
    return split_data(data)

def run_cli():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Create SQL database from compiled CSVs.")

    parser.add_argument(
        'data-directory',
        help='Directory containing CSVs of data and metadata.'
    )
    parser.add_argument(
        'database-url',
        help="e.g. sqlite:////some/dir/database.db.\nSee https://docs.sqlalchemy.org/en/20/core/engines.html#database-urls"
    )
    parser.add_argument(
        '-r', '--remove-exant',
        help='Set flag to delete exant database file before creating a new version.',
        dest='remove', action='store_true', default=False
    )
    args = vars(parser.parse_args())
    print(args)
    db_path = args['database-url'].split('//', maxsplit=1)[-1]
    if os.path.exists(db_path) and os.path.isfile(db_path):
        if args['remove']:
            os.remove(db_path)
        else:
            print(f'Output database file {db_path} exists. Run with -r or use a different file name.')
            exit(1)
    _data = load_csv(args['data-directory'])
    create_database(**_data, db_url=args['database-url'])

def tabulate_experiments_metadata(experiment_details:list[pd.DataFrame]) \
        -> pd.DataFrame:
    """Takes information from the "Experiment details" sheet, returns DF formated
    to be added to the database."""

    column_renamer = {
        "Experiment name":"stringid",
        'Analysis name': 'stringid', # old label for experiment name
        'Library':'library',
        'Multiplicity of infection':'moi',
        'Representation average':'representation',
        'Experiment description (a few sentances)': 'description',
        'Experiment description (a few sentences)':'description',
        'Experiment description': 'description',
        'Notes':'notes',
        'DOI':'doi',
        # Citation in the workbook is mislabeled. Citation in the app generated from this
        "Citation":"reference",
        "Reference":'reference',
        'Date screen completed (yyyy-mm-dd)':'date',
        'Date screen completed':'date',
        'Date published':'date',
    }

    for expmet in experiment_details:
        df_rename_columns(expmet, column_renamer, inplace=True, axis='index')

    experiment_details_table = pd.DataFrame(experiment_details)

    # drop columns not in the mapper
    cols = list(set(column_renamer.values()))
    experiment_details_table = experiment_details_table.loc[:, cols]


    #experiment_details_table.set_index('stringid', inplace=True, )

    # Short citation is generated from the full reference, automatic detection of clashing cites
    def find_date(reference):
        """Pull the date from a reference in the MLA style."""
        pos_dates = set()
        for d in reference.split(' (')[1:]:
            d = d.split(')')[0]
            try:
                d = int(d)
                if 1900 < d < 2100:
                    pos_dates.add(d)
            except:
                pass

        pos_dates = list(pos_dates)
        if len(pos_dates) == 1:
            return pos_dates[0]
        logging.warning(f"Multiple possible dates for reference '{reference}'. No citation created.")
        return '????'

    def short_cite_str(cite):
        year = find_date(cite)
        auth = cite.split(',')[0]
        return f"{auth} ({year})"

    # if there's no external data, and the details excel were put together with older template
    #   it's possible that we won't have a reference column in the whole table
    if "reference" not in experiment_details_table.columns:
        experiment_details_table.loc[:, 'reference'] = np.nan

    # should only be internal
    isinternal = lambda s: 'internal' in str(s).lower()
    noref = experiment_details_table.reference.isna() | experiment_details_table.reference.apply(isinternal)
    # fill with exp ID which is the index
    experiment_details_table.loc[noref, 'reference'] = experiment_details_table.loc[noref].index
    experiment_details_table.loc[noref, 'citation'] = experiment_details_table.loc[noref].index

    experiment_details_table.loc[~noref, 'citation'] = experiment_details_table.loc[~noref, 'reference'].apply(short_cite_str)

    d = experiment_details_table['date'].str.replace('/', '-')
    experiment_details_table['date'] = pd.to_datetime(d, format='mixed', dayfirst=True)

    experiment_details_table = experiment_details_table.infer_objects()

    return experiment_details_table


@dataclasses.dataclass
class DataPaths:
    experiment_id: str
    details: Path|str
    #counts: Path|str
    results: dict[AnalysisType, Path | str]



def get_paths_simons_structure_v1(
        source_dirs:list[str|Path],
        count_dir='./cts',
        details_dir='./det',
        results_dir='./res',
        analysis_filename_prefix='result.'
) -> list[DataPaths]:
    """Return list of DataPaths.

    Uses the experiment specific directory name as the expierment ID.

    Assumes a specific directory structure which was true in Feb 2024. Assumes .gz for some
    files.

    File naming scheme is inconsistent, so it's assumed that a single file in the details
    directory is the right file, if there's zero or >1 files it'll raise an exception.

    Pulls only the reannotated results if they exist, harmonised otherwise. Warns if neither
    are found."""

    logging.info(str(source_dirs))

    def get_fn_or_crash(p:Path):
        fns = os.listdir(p)
        if len(fns) != 1:
            raise ValueError(f"Wrong number of files in {p}")
        return p/fns[0]

    data = []
    for source in source_dirs:
        source = pathlib.Path(source)
        expid = source.stem
        basepath = source/'re'

        tablepath = basepath / results_dir/'tables'
        if not os.path.isdir(tablepath):
            tablepath = basepath / 'harm'/results_dir/'tables'
            if not os.path.isdir(tablepath):
                logging.warning(f'No "harm" or "re" dir found in {str(basepath)}, skipping')
                continue
            logging.warning(f'No "re" results found in {str(basepath)}, using "harm" instead')

        resultspaths = {}
        for ans in ANALYSESTYPES:
            # since this is Simon's structure, assume .gz
            tabfn = tablepath / (analysis_filename_prefix+f'{ans.name}_table.csv.gz')
            if os.path.isfile(tabfn):
                resultspaths[ans] = tabfn
            else:
                logging.warning(f"No results table found at {str(tabfn)}")
        d = DataPaths(
            experiment_id=expid,
            details=get_fn_or_crash(basepath/details_dir),
            #CrisprCounts(get_fn_or_crash(basepath/'cts'))
            results=resultspaths,
        )
        data.append(d)

    return data

def add_statistics(table_paths:dict[AnalysisType, Path|str], exeriment_id:str, session:Session):
    """Take a table, as output by crispr_pipeline, add rows to StatTable."""
    for anstype, fn in table_paths.items():

        multi_table = pd.read_csv(fn, header=[0, 1], index_col=0)
        for cmp in multi_table.columns.levels[0]:
            table = multi_table[cmp].reset_index()
            df_rename_columns(table, {'id': 'gene_id', 'gene': 'gene_id', 'lfc': 'score', 'normZ': 'score',
                                      'fdr_log10': 'fdr10'}, inplace=True)
            table.loc[:, 'comparison_id'] = f"{exeriment_id}.{cmp}"
            table.loc[:, 'analysis_type_id'] = anstype.id
            table.loc[:, 'experiment_id'] = exeriment_id

            insert_df_rows(
                StatTable,
                table.loc[:, ['gene_id', 'score', 'fdr', 'fdr10', 'pos_p', 'neg_p',
                              'comparison_id', 'analysis_type_id','experiment_id']],
                session
            )

def add_data(datapaths:list[DataPaths], engine:Engine, overwrite=False):
    #todo check for exant expid and skip if skipping
    analyses_metadata = [AnalysisWorkbook(d.details) for d in datapaths]
    experiment_details_table = tabulate_experiments_metadata(
        [wb.experiment_details for wb in analyses_metadata]
    )

    with Session(engine) as session:

        insert_df_rows(ExperimentTable, experiment_details_table, session)

        for xp in datapaths:
            add_statistics(xp.results, xp.experiment_id, session)

        session.commit()



    # things should be added in a single transaction, per experiment at the very least



if __name__ == '__main__':
    # Input is list of directory paths that lead to dir with "cts", 'det' and 'res' folders.
    # function takes one folder path and
    #   validates metadata
    #   creates all rows for tables
    #   1.
    # Symbol mapping: I'm going to use what Exorcise says and just pull prev and id from HGNC table.
    #   future might make sense for exorcise to supply all that information, or have script to update
    #   symbols in database from HGNC_table.
    # If experiment ID exists, by default skip and do not alter DB, alternatively drop all rows associated with that experiment
    # Add rows

    import os, pathlib
    import pandas as pd

    rute = pathlib.Path('/Users/thomas03/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofCambridge/Simon Lam - ddrcs/runs')
    srces = [rute/p for p in os.listdir(rute,) if not p.startswith('.')]
    srces = srces[:2]
    print(srces)

    ngin = create_engine_with_schema(echo=False)

    dta = get_paths_simons_structure_v1(srces)
    add_data(dta, ngin)

    #
    with Session(ngin) as S:
        results_path = S.query(ExperimentTable).limit(5).all()
        for result in results_path:
            print(result)






