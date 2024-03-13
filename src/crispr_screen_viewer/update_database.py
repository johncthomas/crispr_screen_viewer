#!/usr/bin/env python
import dataclasses
from glob import glob

import requests, json
import os.path, typing
import pathlib
from pathlib import Path
import unicodedata

from typing import Literal

import pandas as pd
from sqlalchemy import Engine, create_engine, select
from sqlalchemy.orm import Session

from crispr_screen_viewer.database import *

from crispr_screen_viewer import functions_etc
from crispr_screen_viewer.functions_etc import (
    df_rename_columns,
    normalise_text,
    load_stats_csv,
    get_resource_path,
    set_loguru_level
)
from crispr_screen_viewer.dataset import ANALYSESTYPES, AnalysisType

from crispr_tools.data_classes import AnalysisWorkbook, CrisprCounts
from crispr_tools.crispr_pipeline import analysis_tabulate
analysis_tabulate:dict[str, typing.Callable]

import numpy as np

from loguru import logger


def is_nt(s):
    nts = {'None', 'DMSO', '', 'WT', 'NT'}
    return pd.isna(s) or (s in nts)

def get_treatment_str(samp_deets:pd.DataFrame, ctrl:str, treat:str):
    """Return a string describing the treatment performed between the
    control and test samples.

    Args:
        samp_deets: dataframe, rows indexed by samp IDs, giving sample
            metadata.
        ctrl, treat: Sample IDs of a comparison.


    Divides possible treatments into chemical or KO;
    outputs a string for a chemical treatment, a gene-KO
    treatment, a chemical treatment in a KO background,
    or a KO in a chemical background. Depending o """

    fatarrow = '➤'

    t_deets, c_deets = samp_deets.loc[treat], samp_deets.loc[ctrl]

    # get chem treat str
    if not is_nt(t_deets.Treatment):
        if is_nt(c_deets.Treatment) or (c_deets.Treatment == t_deets.Treatment):
            chem_str = t_deets.Treatment
        # if both ctrl/treat samps are different treatments
        else:
            chem_str = f"{c_deets.Treatment}{fatarrow}{t_deets.Treatment}"
    else:
        chem_str = ''

    # Get the KO str.
    #  We assume that we'll never have a KO control and WT treat samp
    if is_nt(t_deets.KO):
        ko_str = ''
    else:
        if is_nt(c_deets.KO) or (c_deets.KO == t_deets.KO):
            ko_str = t_deets.KO+'-KO'
            # if it's just in a TP53 background we don't care that much
            if ko_str == 'TP53-KO':
                ko_str = ''
        else:
            ko_str = f"{c_deets.KO}-KO{fatarrow}{t_deets.KO}-KO"

    if chem_str and not ko_str:
        treatment_str = chem_str
    elif ko_str and not chem_str:
        treatment_str = ko_str
    elif not ko_str and not chem_str:
        treatment_str = 'No treatment'
    else:
        if t_deets.Treatment == c_deets.Treatment:
            # it's a KO in a chemical background
            treatment_str = f"{ko_str} (with {chem_str})"
        else:
            treatment_str = f"{chem_str} (in {ko_str} cells)"

    return treatment_str

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


def insert_records(tablecls, records:list[dict], session:Session):
    """Insert rows into given table. Record dict.keys must be present in tablecls
    columns. Use DataFrame.to_dict(orient='records') to pass DF info.

    NaN are skipped over.

    Fails if they've already been added.

    Does not commit"""
    # https://docs.sqlalchemy.org/en/20/tutorial/data_update.html

    # filter the records so key:value pairs of null value are removed
    session.add_all(
        [tablecls(**{k: v for k, v in kw.items() if not pd.isna(v)})
         for kw in records]
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


# def get_exp_ids(engine:Engine) -> dict[str, int]:
#     with Session(engine) as S:
#         rows = S.execute(select(ExperimentTable.stringid, ExperimentTable.id)).all()
#
#     return dict(rows) # noqa
#
# def get_comp_ids(engine:Engine) -> dict[str, int]:
#     with Session(engine) as S:
#         rows = S.execute(select(ComparisonTable.stringid, ComparisonTable.id)).all()
#
#     return dict(rows) # noqa


# def format_comparisons_df(
#         comparisons_metadata:pd.DataFrame,
#         experiment_ids:dict[str, int]=None,
# ) -> pd.DataFrame:
#
#     comparison_column_mapping = {
#         'Comparison ID': 'stringid',
#         'Treatment': 'treatment_label',
#         'Timepoint': 'timepoint',
#         'Cell line': 'cell',
#         'Ctrl samp': 'ctrl',
#         'Treat samp': 'treat',
#         'KO': 'ko',
#         'Dose': 'dose',
#         'Growth inhibition %': 'gi',
#         'Days grown': 'days_grown',
#         'Library': 'library',
#         'analyses_bitmask': 'analyses_bitmask',
#         'Experiment ID': 'experiment'
#     }
#
#     comps = comparisons_metadata.reset_index()
#     #comps.index.name = 'id'
#     comps.head()
#
#     comps.loc[:, 'analyses_bitmask'] = comps['Available analyses'].apply(
#         lambda s: ANALYSESTYPES.encode_bitmask([ss for ss in s.split('|')])
#     )
#
#     comps = prep_deets_columns(comps, comparison_column_mapping)
#     # comps = comps[comparison_column_mapping.keys()]
#     # cols = functions_etc.df_rename_columns(comps, comparison_column_mapping, )
#     # comps.columns = cols
#
#     if experiment_ids is not None:
#         comps.experiment = comps.experiment.map(experiment_ids)
#
#     return comps


# def format_statistics_dfs(
#         statistics:dict[str, pd.DataFrame],
#         comparison_ids:dict[str, int]=None,
#         analyses=('mageck', 'drugz'),
#         stats_names=('pos_p', 'neg_p', 'fdr', 'score')
# ) -> list[pd.DataFrame]:
#     """Convert wide tables in dict keyed by things like "mag_pos_p" or "drz_score"
#     to list of long tables for insertion into DB.
#     """
#
#     stats_dfs = []
#
#     for method_name in analyses:
#         method = ANALYSESTYPES[method_name]
#
#         stt_values = {}
#         for stt in stats_names:
#             d = statistics[f"{method.shortname}_{stt}"]
#
#             # ...NaN in the index can cause crashes. Why can there be a NaN in the index?
#             #   found in a few mageck results tables.
#             d = d.loc[~d.index.isna()]
#             stt_values[stt] = d.unstack()
#
#         df = pd.DataFrame(stt_values).reset_index()
#
#         functions_etc.df_rename_columns(
#             df,
#             {'Comparison ID': 'comparison_id', 'id':'gene_id', 'gene':'gene_id', 'level_0':'comparison_id'},
#             inplace=True
#         )
#
#         df.loc[:, 'fdr10'] = df.fdr.apply(lambda p: -np.log10(p))
#
#         if comparison_ids is not None:
#             df.comparison_id = df.comparison_id.map(comparison_ids)
#
#         df.loc[:, 'analysis_type_id'] = method.id
#
#         stats_dfs.append(df)
#
#     return stats_dfs


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


# def split_data(data:dict[str, pd.DataFrame]):
#     out_data = {}
#     for k in ('experiments_metadata', 'comparisons_metadata'):
#         out_data[k] = data[k]
#         del data[k]
#
#     out_data['stats_tables'] = data
#
#     return out_data

# def load_csv(paff:str):
#     """Load all .csv in paff, should be experiments_metadata.csv, comparisons_metadata.csv and
#     individual tables for the statistics. The stats tables get put under "stats_tables" key. """
#     data: dict[str, pd.DataFrame] = {
#         fn.replace('.csv', ''): pd.read_csv(pathlib.Path(paff) / fn, index_col=0)
#         for fn in os.listdir(paff) if fn.endswith('.csv')}
#     return split_data(data)

# def run_cli():
#     from argparse import ArgumentParser
#     parser = ArgumentParser(description="Create SQL database from compiled CSVs.")
#
#     parser.add_argument(
#         'data-directory',
#         help='Directory containing CSVs of data and metadata.'
#     )
#     parser.add_argument(
#         'database-url',
#         help="e.g. sqlite:////some/dir/database.db.\nSee https://docs.sqlalchemy.org/en/20/core/engines.html#database-urls"
#     )
#     parser.add_argument(
#         '-r', '--remove-exant',
#         help='Set flag to delete exant database file before creating a new version.',
#         dest='remove', action='store_true', default=False
#     )
#     args = vars(parser.parse_args())
#     print(args)
#     db_path = args['database-url'].split('//', maxsplit=1)[-1]
#     if os.path.exists(db_path) and os.path.isfile(db_path):
#         if args['remove']:
#             os.remove(db_path)
#         else:
#             print(f'Output database file {db_path} exists. Run with -r or use a different file name.')
#             exit(1)
#     _data = load_csv(args['data-directory'])
#     create_database(**_data, db_url=args['database-url'])


@dataclasses.dataclass
class AnalysisInfo:
    experiment_id:str
    analysis_workbook: AnalysisWorkbook
    #counts: Path|str
    results_path: dict[AnalysisType, Path | str]



def get_paths_exorcise_structure_v1(
        source_dirs:list[str|Path],
        details_dir='./det',
        results_dir='./res',
        #count_dir='./cts',
        analysis_filename_prefix='result.'
) -> list[AnalysisInfo]:
    """Return list of DataPaths.

    Splits out the experiment ID with `expid.split('/')[1] if the ID in AnalysisWorkbook
    contains '/'.

    Assumes a specific directory structure which was true in Feb 2024. Assumes .gz for some
    files.

    It's assumed that a single file in the details directory is the right file, if there's zero
    or >1 files it'll raise an exception.

    Pulls only the reannotated resultsm from dir 're' if it exists, harmonised from 'harm'
    otherwise. Warns if neither are found."""

    logger.debug(str(source_dirs))

    def get_fn_or_crash(p:Path):
        fns = os.listdir(p)
        if len(fns) != 1:
            raise ValueError(f"Wrong number of files in {p}")
        return p/fns[0]

    data = []

    good_dirs = set()
    for d in source_dirs:
        d = Path(d)

        # catch bad files, add other checks as elif
        if d.stem.startswith('.'):
            pass
        elif os.path.isfile(d):
            logger.info(f'Skipping non-directory: {d}')
        # if it's not a file or directory how did it get here?
        elif not os.path.isdir(d):
            logger.warning(f"Path doesn't exist: {d}")
        # else everything's good
        else:
            good_dirs.add(d)


    for source in list(good_dirs):
        try:
            assert os.path.isdir(source)
        except AssertionError as e:
            print(source)
            raise e


        basepath = source / 're'
        if not os.path.isdir(basepath):
            basepath = source / 'harm'
            if not os.path.isdir(basepath):
                logger.warning(f'No "harm" or "re" dir found in {str(basepath)}, skipping')
                continue
            logger.info(f'No "re" results found in {str(basepath)}, using "harm" instead')

        tablepath = basepath/results_dir/'tables'

        resultspaths = {}
        for ans in ANALYSESTYPES:
            # since this is Simon's structure, assume .gz
            tabfn = tablepath / (analysis_filename_prefix+f'{ans.name}_table.csv.gz')
            if os.path.isfile(tabfn):
                resultspaths[ans] = tabfn
            else:
                tabfn = tablepath / (analysis_filename_prefix + f'{ans.name}_table.csv')
                if os.path.isfile(tabfn):
                    resultspaths[ans] = tabfn
                else:
                    logger.warning(f"No results table found at {str(tabfn)}")
        analysiswb = AnalysisWorkbook(get_fn_or_crash(basepath/details_dir))

        xpid = analysiswb.experiment_details['Experiment name']
        if '/' in xpid:
            xpid = xpid.split('/')[1]
            logger.debug(f'New experiment ID is {xpid}')
            analysiswb.experiment_details['Experiment name'] = xpid
            analysiswb.expd['experiment_id'] = xpid

        d = AnalysisInfo(
            experiment_id=xpid,
            analysis_workbook=analysiswb,
            results_path=resultspaths,
        )
        data.append(d)

    return data

def tabulate_experiments_metadata(experiment_details:list[pd.DataFrame]) \
        -> pd.DataFrame:
    """Takes information from the "Experiment details" sheet, returns DF formated
    to be added to the database."""

    column_renamer = {
        "Experiment name":"Experiment ID",
        'Analysis name': 'Experiment ID', # old label for experiment name
        #'Library':'Library',
        #'Multiplicity of infection':'MOI',
        #'Representation average':'Representation',
        'Experiment description (a few sentances)': 'Experiment description',
        'Experiment description (a few sentences)':'Experiment description',
        #'Experiment description': 'Experiment description',
        # 'Notes':'notes',
        # 'DOI':'doi',
        # Citation in the workbook is mislabeled. Citation in the app generated from this
        "Citation":"Reference",
        # "Reference":'reference',
        'Date screen completed (yyyy-mm-dd)':'Date',
        'Date screen completed':'Date',
        'Date published':'Date',
    }

    experiment_details = [xpmet.copy() for xpmet in experiment_details]
    logger.debug(experiment_details)
    for expmet in experiment_details:
        df_rename_columns(expmet, column_renamer, inplace=True, axis='index')

    experiment_details_table = pd.DataFrame(experiment_details)
    experiment_details_table = experiment_details_table.reset_index(drop=True)

    # drop columns not in the mapper
    #cols = list(set(column_renamer.values()))
    #experiment_details_table = experiment_details_table.loc[:, cols]

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
        logger.warning(f"Multiple or zero possible dates for reference '{reference}', possibilities '{pos_dates}'. No citation created.")
        return '????'

    def short_cite_str(reference):
        """Get "{first_author} ({year})" from reference, if possible."""
        # remove weird formating characters.
        logger.debug(f"Reference {reference}")
        reference = normalise_text(reference)
        year = find_date(reference)
        auth = reference.split(',')[0]
        return f"{auth} et al ({year})"

    # if there's no external data, and the details excel were put together with older template
    #   it's possible that we won't have a reference column in the whole table
    if "Reference" not in experiment_details_table.columns:
        experiment_details_table.loc[:, 'Reference'] = ''

    # should only be internal
    isinternal = lambda s: 'internal' in str(s).lower()
    noref = experiment_details_table.Reference.isna() | experiment_details_table.Reference.apply(isinternal)
    # fill with exp ID which is the index
    experiment_details_table.loc[noref, 'Reference'] = experiment_details_table.loc[noref].index
    experiment_details_table.loc[:, 'Citation'] = ''
    experiment_details_table.loc[noref, 'Citation'] = experiment_details_table.loc[noref].index

    experiment_details_table.loc[~noref, 'Citation'] = experiment_details_table.loc[~noref, 'Reference'].apply(short_cite_str)

    bad_cite = experiment_details_table.loc[experiment_details_table.Citation.str.contains('?', regex=False), 'Experiment ID'].values

    if len(bad_cite) > 0:
        logger.warning(f"The following experiments have bad citation info:\n{'\n\t'.join(bad_cite)}")

    # deduplicate citations
    cites = experiment_details_table['Citation']
    if cites.duplicated().any():
        for ref, indicies in cites.groupby(cites).groups.items():
            if len(indicies) > 1:
                for n, i in enumerate(indicies):
                    # starting from 'a'...
                    letter = chr(97 + n)
                    cites.loc[i] = ref.replace(')', letter + ')')

    # d = experiment_details_table['date'].str.replace('/', '-')
    # experiment_details_table['date'] = pd.to_datetime(d, format='mixed', )

    experiment_details_table = experiment_details_table.infer_objects()

    return experiment_details_table

def add_experiments(data_paths:list[AnalysisInfo], session:Session):
    experiment_details = [d.analysis_workbook.experiment_details for d in data_paths]
    table = tabulate_experiments_metadata(experiment_details)
    logger.info(f"Adding {table.shape[0]} rows to ExperimentTable")
    insert_records(ExperimentTable,  table.to_dict(orient='records'), session)

def tabulate_comparisons(analysis_wb:AnalysisWorkbook):
    # comparison_column_mapping = {
    #     'Treatment': 'treatment_label',
    #     'Timepoint': 'timepoint',
    #     'Cell line': 'cell',
    #     'Ctrl samp': 'ctrl',
    #     'Treat samp': 'treat',
    #     'KO': 'ko',
    #     'Dose': 'dose',
    #     'Growth inhibition %': 'gi',
    #     'Days grown': 'days_grown',
    #     'Library': 'library',
    #     'Experiment ID': 'experiment',
    #     'Notes':'notes'
    # }
    already_warned_of_missing_column = set()
    # this will become the returned table
    comparisons_metadata = []

    # iterate through comparisons, pull the sample data for controls/treatments and build the table
    for group_name, (ctrl, treat) in analysis_wb.iter_comps():
        comp_row = {}
        # print(ctrl, treat)
        comp_row['control_sample'] = treat
        comp_row['test_sample'] = ctrl
        comp_row['Treatment'] = get_treatment_str(
            analysis_wb.wb['Sample details'],
            ctrl, treat
        )
        ctrl_row = analysis_wb.samples.loc[ctrl]
        comp_row['control_treatment'] = ctrl_row['Treatment']
        comp_row['control_ko'] = ctrl_row['KO']

        for k in ('Dose', 'Growth inhibition %', 'Days grown', 'Cell line', 'KO', 'Notes'):
            # Deal with columns being dropped from the input data.
            if k in analysis_wb.samples.columns:
                comp_row[k] = analysis_wb.samples.loc[treat, k]
            else:
                if k not in already_warned_of_missing_column:
                    already_warned_of_missing_column.add(k)
                    logger.warning(f"Missing columns '{k}' in workbook for experiment {analysis_wb.expd['experiment_id']}")
                comp_row[k] = np.nan

        exp_id = analysis_wb.expd['experiment_id']
        comp_row['Experiment ID'] = exp_id
        comp_row['Library'] = analysis_wb.experiment_details['Library']
        comp_row['stringid'] = f"{exp_id}.{ctrl}-{treat}"

        # format strings
        if not pd.isna(comp_row['Dose']):
            comp_row['Dose'] = str(comp_row['Dose']).replace('uM', 'μM')
        if pd.isna(comp_row['KO']) or (comp_row['KO'] == ''):
            comp_row['KO'] = 'WT'

        comp_row['Timepoint'] = group_name.split('_')[0]

        comparisons_metadata.append(comp_row)

    return pd.DataFrame(comparisons_metadata)


def add_comparisons(analyses_info:list[AnalysisInfo], session:Session):
    for info in analyses_info:
        comparisons_metadata = tabulate_comparisons(info.analysis_workbook)
        logger.info(f"Adding {comparisons_metadata.shape[0]} rows to ComparisonTable from {info.experiment_id}")
        insert_records(ComparisonTable, comparisons_metadata.to_dict(orient='records'), session)


def get_gene_symbols_db(session) -> set[str]:
    gns = set([g[0] for g in session.query(GeneTable.symbol).distinct().all()])
    return gns

def gene_info_from_refseq_by_symbols(
        symbols:list[str],
        organism:Literal['Human']|Literal['Mouse']|int,
) -> list[dict[str,str]]:
    """Get list of information from RefSeq, formatted to pass to GeneTable."""
    if type(organism) is int:
        orgid = str(organism)
    else:
        orgid = {
            'Human':"9606",
            'Mouse':"10090",
        }[organism]

    #might be a set or something
    symbols = list(symbols)

    # seems to have ~1000 gene query limit
    gene_info: list[dict[str, str]] = []
    done = False
    batch_size = 1000
    chunk = 0
    n_symbols = len(symbols)
    logger.info(f"Querying refseq with {n_symbols} genes in batches of {batch_size}.")
    while not done:
        start, stop = batch_size*chunk, batch_size*(chunk+1)
        chunk += 1
        symb_chunk = symbols[start:stop]
        if stop > n_symbols:
            done = True

        payload = {
            "symbols_for_taxon": {
                "symbols": symb_chunk,
                "taxon": orgid}
        }

        logger.debug(payload)

        refseqres = requests.request(
            'POST',
            'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene',
            data=json.dumps(payload),
        )

        refseqres.raise_for_status()



        for r in refseqres.json()['reports']:
            gn_res = r['gene']
            symbol:str = gn_res['symbol']
            query = r['query'][0]

            # Skip synonyms.
            # We are working with the assumption that the symbol comes from refseq
            if not symbol == query:
                continue

            if 'synonyms' in gn_res:
                symonyms = gn_res['synonyms']
            else:
                symonyms = []

            try:
                official_id = gn_res['nomenclature_authority']['identifier']
                symonyms = gn_res['synonyms']
            except KeyError:
                logger.debug(f"Key error getting information from refseq results, for gene symbol {gn_res['symbol']}.")
                continue

            gninfo = {
                'id':symbol,
                'symbol':symbol,
                'ncbi':'NCBI: '+str(gn_res['gene_id']),
                'official_id':official_id,
                'synonyms_str':str(symonyms),
                'organism':organism
            }

            gene_info.append(gninfo)

    return gene_info

def add_genes_from_symbols(
        symbols: list[str],
        organism: Literal['Human'] | Literal['Mouse'] | int,
        session:Session,
        query_refseq_by_symbol=False,
):
    """Look up the symbol in refseq, get some IDs"""

    symbols_to_add = set(symbols).difference(get_gene_symbols_db(session))
    if query_refseq_by_symbol:
        #logger.debug('Adding: '+str(symbols_to_add))
        records = gene_info_from_refseq_by_symbols(list(symbols_to_add), organism)
        insert_records(GeneTable, records, session)

    # above queries will not add rows for genes that aren't in refseq
    new_current_symbols = get_gene_symbols_db(session)
    no_refseq_genes = symbols_to_add.difference(new_current_symbols)
    logger.debug('Num ID-less genes being added: '+str(len(no_refseq_genes)))
    empty_records = [dict(id=s, symbol=s, organism=organism) for s in no_refseq_genes]
    #logger.debug(empty_records)
    insert_records(GeneTable, empty_records, session)

def tabulate_statistics(info:AnalysisInfo) -> pd.DataFrame:
    experiment_id = info.experiment_id
    tables = []
    for analysis_type, fn in info.results_path.items():
        logger.debug(f"tabulating from {fn}")
        stats_table = load_stats_csv(fn)
        stats_table.index.name = 'gene_id'

        for cmp in stats_table.columns.levels[0]:
            table = stats_table[cmp].reset_index()
            df_rename_columns(table, {'lfc': 'score', 'normZ': 'score',
                                      'fdr_log10': 'fdr10'}, inplace=True)
            table.loc[:, 'comparison_id'] = f"{experiment_id}.{cmp}"
            table.loc[:, 'analysis_type_id'] = analysis_type.id
            table.loc[:, 'experiment_id'] = experiment_id
            table = table.loc[:, ['gene_id', 'score', 'fdr', 'fdr10', 'pos_p', 'neg_p',
                                  'comparison_id', 'analysis_type_id', 'experiment_id']]
            tables.append(table)
    return pd.concat(tables)

def add_statistics(analysesinfo:list[AnalysisInfo], session:Session, query_refseq_by_symbol=False):
    """Take a table, as output by crispr_pipeline, add rows to StatTable."""
    for info in analysesinfo:
        table = tabulate_statistics(info)

        add_genes_from_symbols(
            table.gene_id,
            info.analysis_workbook.experiment_details.Organism,
            session,
            query_refseq_by_symbol=query_refseq_by_symbol
        )

        logger.info(f"Adding {table.shape[0]} rows to StatTable from {info.experiment_id}")
        insert_records(
            StatTable,
            table.to_dict(orient='records'),
            session
        )

def write_metadata_tables(analysesinfo:list[AnalysisInfo], outdir):
    # This is a bit messy currently - the tabulate functions are written with database names
    #    for the columns, all lower case, no spaces, but currently the code expects dataframes
    #    with the columns named as displayed on the website, so here we build the tables and
    #    then rename the columns. Later all functions should be rewritten to use the SQL DB
    outdir = Path(outdir)
    comparisons_metadata = []
    for info in analysesinfo:
        comparisons_metadata.append(tabulate_comparisons(info.analysis_workbook))
    comparisons_metadata = pd.concat(comparisons_metadata)
    db_to_df_header = {
        'treatment_label': 'Treatment',
         'timepoint': 'Timepoint',
         'cell': 'Cell line',
         'ctrl': 'Ctrl samp',
         'treat': 'Treat samp',
         'ko': 'KO',
         'dose': 'Dose',
         'gi': 'Growth inhibition %',
         'days_grown': 'Days grown',
         'library': 'Library',
         'experiment': 'Experiment ID',
         'notes': 'Notes',
        'stringid':'Comparison ID'
    }

    df_rename_columns(comparisons_metadata, db_to_df_header, inplace=True)

    experiments_metadata = tabulate_experiments_metadata(
            [d.analysis_workbook.experiment_details for d in analysesinfo]
    )

    df_rename_columns(
        experiments_metadata,
        {'stringid': 'Experiment ID',
         'library': 'Library',
         'moi': 'Multiplicity of infection',
         'representation': 'Representation average',
         'description': 'Experiment description',
         'notes': 'Notes',
         'doi': 'DOI',
         'reference': 'Reference',
         'date': 'Date published'},
        inplace=True,
    )

    comparisons_metadata.set_index('Comparison ID', drop=False, inplace=True)
    experiments_metadata.set_index('Experiment ID', drop=False, inplace=True)


    comparisons_metadata.to_csv(outdir/'comparisons_metadata.csv.gz')
    experiments_metadata.to_csv(outdir/'experiments_metadata.csv.gz')

def add_data_to_database(
        analysesinfo:list[AnalysisInfo],
        engine:Engine,
        outdir:str|Path,
        overwrite=False,
        query_refseq_by_symbol=True):
    #todo check for exant expid and skip if skipping
    #todo convert comparisons and experiment details to use SQL

    # if overwrite is True:
    #     raise NotImplementedError("Not currently able to update the database in place")

    outdir = Path(outdir)

    write_metadata_tables(analysesinfo, outdir)

    with Session(engine) as session:
        # add_experiments(analysesinfo, session) # note, only actual thing this is currently used for is doi and reference
        # add_comparisons(analysesinfo, session)
        add_statistics(analysesinfo, session, query_refseq_by_symbol=query_refseq_by_symbol)
        logger.info("Commiting changes")
        session.commit()

def create_database(
        outdir,
        analysis_infos:list[AnalysisInfo],
        refseq=True,
        ask_before_deleting=True,
        run_server=False,
        port=8050,
        max_analyses_for_testing:int=None,
) -> None:
    outdir = Path(outdir)
    import glob

    outdir.mkdir(exist_ok=True)

    files = [fn for fn in glob.glob(str(outdir) + '/*') if os.path.isfile(fn)]
    filestr = '\n'.join(files)
    if files:
        if ask_before_deleting:
            print(f"Output dir has files, these will ALL be deleted: \n {filestr}")
            input("Press enter to continue. Ctrl+C to cancel")
        for f in files:
            os.remove(f)

    engine_url = f'sqlite:///{str(outdir)}/database.db'

    sql_engin = create_engine_with_schema(
        engine_url
    )

    analysis_infos=analysis_infos[:max_analyses_for_testing]

    add_data_to_database(analysis_infos, sql_engin, outdir, query_refseq_by_symbol=refseq)

    write_metadata_tables(analysis_infos, outdir)

    if run_server:
        from crispr_screen_viewer.launch import init_app
        app = init_app(
            str(outdir),
            engine_url,
            debug_messages=True,

        )

        app.run_server(debug=False, host='0.0.0.0', port=port, )


def create_test_database(**kwargs):
    outdir = get_resource_path('data/test_db')
    g = glob(f'{get_resource_path("tests/test_data/exorcise_style")}/*')
    logger.debug(g)
    infos = get_paths_exorcise_structure_v1(g)
    logger.debug(infos)

    actual_kwargs = dict(refseq=False, ask_before_deleting=False) | kwargs
    create_database(outdir, infos, **actual_kwargs)

def __create_database_20240228():
    stemd = '2024-02-28'
    outd = Path(f'/Users/thomas03/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/ddrcs/app_data/{stemd}')

    xclude = ['ParrishBerger2021', 'SchleicherMoldovan2020_Ca', ]

    src = Path(
        '/Users/thomas03/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofCambridge/Simon Lam - ddrcs/runs'
    )

    drs = [src/d for d in os.listdir(src) if d not in xclude]

    # print(drs)
    analysis_infos = get_paths_exorcise_structure_v1(drs)

    create_database(outd, analysis_infos[:], refseq=True, run_server=False)


# def check_for_null_genes(analysis_infos:list[AnalysisInfo]):
#     for a in analysis_infos:
#         tab = tabulate_statistics(a)
#         if tab.gene_id.isna().any():
#             print(a.experiment_id)



if __name__ == '__main__':
    set_loguru_level(logger, 'INFO')
    __create_database_20240228()
    #__create_database_20240228()
    # exptab = tabulate_experiments_metadata(
    #     [xp.analysis_workbook.wb for xp in get_paths_simons_structure_v1(['./tests/test_data/exorcise_style/test1'])]
    # )
    # print(exptab)
    pass


