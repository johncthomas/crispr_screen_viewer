#!/usr/bin/env python

import os.path, typing
import pathlib

import pandas as pd
from sqlalchemy import Engine, create_engine, select
from sqlalchemy.orm import Session

from database import *

from crispr_screen_viewer import functions_etc
from crispr_screen_viewer.dataset import ANALYSESTYPES

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



def insert_df_rows(tablecls, df, engine):
    """Take row values from df, and add them as rows to the table.

    NaN are skipped over.

    Assumes everything is formated correctly with df.columns matching field names.
    Fails if they've already been added, dunno how you update them"""
    # https://docs.sqlalchemy.org/en/20/tutorial/data_update.html
    with Session(engine) as S:
        # filter the records so key:value pairs of null value are removed
        S.add_all(
            [tablecls(**{k: v for k, v in kw.items() if not pd.isna(v)})
             for kw in df.to_dict(orient='records')]
        )
        S.commit()


def format_experiments_df(
        experiments_metadata:pd.DataFrame,
) -> pd.DataFrame:
    """Reformat experiments DF for insertion into database"""
    # experiment
    experiment_column_mapping = {
        'Experiment ID': 'stringid',
        'Date published': 'date',
        'Library': 'library',
        'Multiplicity of infection': 'moi',
        'Representation average': 'representation',
        'Experiment description (a few sentences)': 'description',
        'Notes': 'notes',
        'DOI': 'doi',
        'Reference': 'reference',
        'Source': 'source',
        'Citation': 'citation',
    }

    exp = experiments_metadata.reset_index()
    exp = prep_deets_columns(exp, experiment_column_mapping)
    d = exp.date.str.replace('/', '-')
    exp.date = pd.to_datetime(d, format='mixed', dayfirst=True)

    return exp


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


def create_database(
        experiments_metadata:pd.DataFrame,
        comparisons_metadata:pd.DataFrame,
        stats_tables:dict[str, pd.DataFrame],
        db_url:str,
        use_autogen_ids=False,
):
    """Create a SQL database experimental results from dataframes.

    `stats` organised keyed by experiment method short name ('mag', 'drz') then
    by stat type ('pos_p', 'fdr', etc)."""
    if use_autogen_ids:
        raise NotImplementedError("Currently using string IDs from experiment/sample names, "
                                  "using autogenerated IDs will require the schema to be changed.")
    exp = format_experiments_df(experiments_metadata)

    engine = create_engine_with_schema(
        db_url,
    )

    # needs to be inserted first to generate IDs
    insert_df_rows(ExperimentTable, exp, engine)


    exp_ids = get_exp_ids(engine) if use_autogen_ids else None
    comps = format_comparisons_df(comparisons_metadata, exp_ids)
    insert_df_rows(ComparisonTable, comps, engine)

    comp_ids = get_comp_ids(engine) if use_autogen_ids else None
    stats = format_statistics_dfs(stats_tables, comp_ids)
    for tbl in stats:
        insert_df_rows(StatTable, tbl, engine)


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

if __name__ == '__main__':
    #__create_test_db()
    run_cli()





