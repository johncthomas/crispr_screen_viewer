#!/usr/bin/env python
from dataclasses import dataclass
from glob import glob

import requests, json
import os.path, typing
import pathlib
from pathlib import Path
import unicodedata

from typing import Literal, Tuple

from argparse import ArgumentParser

import pandas as pd
import sqlalchemy
from sqlalchemy import Engine, create_engine, select
from sqlalchemy.orm import Session

from crispr_screen_viewer.database import *

from crispr_screen_viewer import functions_etc
from crispr_screen_viewer.functions_etc import (
    df_rename_columns,
    normalise_text,
    load_stats_csv,
    get_resource_path,
    set_loguru_level,
    get_ith_from_all
)
from crispr_screen_viewer.dataset import (
    ANALYSESTYPES,
    AnalysisType,
    MetadataTables,
    DB_FILES,
    get_db_url
)

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
    or a KO in a chemical background. """

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


def load_test_db_data(d='data/test_db') -> Tuple[Engine, MetadataTables]:
    test_db_dir = get_resource_path(d)
    url = get_db_url(test_db_dir)
    engine = create_engine(url)
    metadata = MetadataTables.from_files(test_db_dir)

    return engine, metadata

def create_test_database(**kwargs):
    """Create a database from ./tests/test_data, output to ./data/test_db
    by default. **kwargs passed to create_database()"""
    outdir = get_resource_path('data/test_db')
    g = glob(f'{get_resource_path("tests/test_data/exorcise_style")}/*')
    logger.debug(g)
    infos = get_paths_exorcise_structure_v1(g)
    logger.debug(infos)

    default_kwargs = dict(
        outdir=outdir,
        analysis_infos=infos,
        #refseq=False,
        ask_before_deleting=False)
    create_database(**default_kwargs | kwargs)

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


@dataclass
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
    #for group_name, (ctrl, treat) in analysis_wb.iter_comps():
    for _, comprow in analysis_wb.control_groups.iterrows():

        comparison_info = {}
        # print(ctrl, treat)
        comparison_info['ControlSample'] = ctrl = comprow['Control sample']
        comparison_info['TestSample'] = treat =  comprow['Test sample']

        # get treatment string
        treatment_label:str = None
        if 'Contrast' in comprow.index:
            contrast = comprow['Contrast']
            if (not pd.isna(contrast)) and contrast:
                treatment_label = contrast
                logger.debug(f"Using given Contrast name '{treatment_label}'")

        if treatment_label is None:
            treatment_label = get_treatment_str(
                analysis_wb.wb['Sample details'],
                ctrl, treat
            )

        comparison_info['Treatment'] = treatment_label

        ctrl_row = analysis_wb.samples.loc[ctrl]
        comparison_info['ControlTreatment'] = ctrl_row['Treatment']
        comparison_info['ControlKO'] = ctrl_row['KO']

        for k in ('Dose', 'Growth inhibition %', 'Days grown', 'Cell line', 'KO', 'Notes'):
            # Deal with columns being dropped from the input data.
            if k in analysis_wb.samples.columns:
                comparison_info[k] = analysis_wb.samples.loc[treat, k]
            else:
                if k not in already_warned_of_missing_column:
                    already_warned_of_missing_column.add(k)
                    logger.warning(f"Missing columns '{k}' in workbook for experiment {analysis_wb.expd['experiment_id']}")
                comparison_info[k] = np.nan

        exp_id = analysis_wb.expd['experiment_id']
        comparison_info['Experiment ID'] = exp_id
        comparison_info['Library'] = analysis_wb.experiment_details['Library']
        comparison_info['Comparison ID'] = f"{exp_id}.{ctrl}-{treat}"

        # format strings
        if not pd.isna(comparison_info['Dose']):
            comparison_info['Dose'] = str(comparison_info['Dose']).replace('uM', 'μM')
        if pd.isna(comparison_info['KO']) or (comparison_info['KO'] == ''):
            comparison_info['KO'] = 'WT'

        comparison_info['Timepoint'] = comprow['Group'].split('_')[0]

        comparisons_metadata.append(comparison_info)

    return pd.DataFrame(comparisons_metadata)


# def add_comparisons(analyses_info:list[AnalysisInfo], session:Session):
#     for info in analyses_info:
#         comparisons_metadata = tabulate_comparisons(info.analysis_workbook)
#         logger.info(f"Adding {comparisons_metadata.shape[0]} rows to ComparisonTable from {info.experiment_id}")
#         insert_records(ComparisonTable, comparisons_metadata.to_dict(orient='records'), session)


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

def create_metadata_tables(analysesinfo:list[AnalysisInfo]) \
        -> MetadataTables:


    comparisons_metadata = []
    for info in analysesinfo:
        comparisons_metadata.append(tabulate_comparisons(info.analysis_workbook))
    comparisons_metadata = pd.concat(comparisons_metadata)
    # db_to_df_header = {
    #     'treatment_label': 'Treatment',
    #      'timepoint': 'Timepoint',
    #      'cell': 'Cell line',
    #      'ctrl': 'Ctrl samp',
    #      'treat': 'Treat samp',
    #      'ko': 'KO',
    #      'dose': 'Dose',
    #      'gi': 'Growth inhibition %',
    #      'days_grown': 'Days grown',
    #      'library': 'Library',
    #      'experiment': 'Experiment ID',
    #      'notes': 'Notes',
    #     'stringid':'Comparison ID'
    # }
    #
    # df_rename_columns(comparisons_metadata, db_to_df_header, inplace=True)

    experiments_metadata = tabulate_experiments_metadata(
            [d.analysis_workbook.experiment_details for d in analysesinfo]
    )

    # df_rename_columns(
    #     experiments_metadata,
    #     {'stringid': 'Experiment ID',
    #      'library': 'Library',
    #      'moi': 'Multiplicity of infection',
    #      'representation': 'Representation average',
    #      'description': 'Experiment description',
    #      'notes': 'Notes',
    #      'doi': 'DOI',
    #      'reference': 'Reference',
    #      'date': 'Date published'},
    #     inplace=True,
    # )

    comparisons_metadata.set_index('Comparison ID', drop=False, inplace=True)
    experiments_metadata.set_index('Experiment ID', drop=False, inplace=True)

    #
    # comparisons_metadata.to_csv(outdir/COMP_CSV)
    # experiments_metadata.to_csv(outdir/EXP_CSV)
    return MetadataTables(comparisons=comparisons_metadata, experiments=experiments_metadata)

def write_db_files(
        outdir:str|Path,
        analysis_infos:list[AnalysisInfo],
        metadata:MetadataTables,
        session:Session
) -> None:
    add_statistics(analysis_infos, session, )
    logger.info("Commiting changes")
    session.commit()
    logger.info("Writing metadata tables")
    metadata.to_files(outdir)

def remove_experiments(
        metadata_tables:MetadataTables,
        exp_ids: typing.Collection[str],
        session:Session,
) -> MetadataTables:

    delete_stmt = sqlalchemy.delete(
        StatTable
    ).where(
        StatTable.experiment_id.in_(exp_ids)
    )
    session.execute(delete_stmt)

    comparisons = metadata_tables.comparisons.loc[
        ~metadata_tables.comparisons['Experiment ID'].isin(exp_ids)
    ]
    experiments = metadata_tables.experiments.loc[
        ~metadata_tables.experiments['Experiment ID'].isin(exp_ids)
    ]

    return MetadataTables(comparisons=comparisons, experiments=experiments)

def remove_experiments_from_db(
        db_dir:str|Path,
        experiments_to_remove:list[str]
):
    db_dir = Path(db_dir)
    metadata = MetadataTables.from_files(db_dir)
    engine = create_engine(get_db_url(db_dir))
    with Session(engine) as session:
        mod_metadata = remove_experiments(metadata, experiments_to_remove, session)
        session.commit()

    mod_metadata.to_files(db_dir)

def update_database(
        db_dir,
        analysis_infos: list[AnalysisInfo],
        refseq:bool=False,
        update_experiments=False,
) -> None:
    db_dir = Path(db_dir)
    metadata = MetadataTables.from_files(db_dir)
    engine = create_engine(get_db_url(db_dir))

    # get exp already in the db,
    new_expid = set([ans.experiment_id for ans in analysis_infos])
    exant_expid = set(metadata.experiments['Experiment ID'].unique())
    overlapping_expid = new_expid.intersection(exant_expid)
    olxp_str = ', '.join(sorted(overlapping_expid))

    with Session(engine) as session:
        #  if update_experiments==False, remove those experiments before adding,
        #    otherwise, remove those exps (do not save metatadatas till the end
        if not update_experiments:
            logger.info(f"Experiments already in the database will be skipped: {olxp_str}")
            analysis_infos = [ans for ans in analysis_infos if ans.experiment_id not in overlapping_expid]
        else:
            logger.info(f"Experiments already in the database will be removed and re-added: {olxp_str}")
            modified_metadata = remove_experiments(
                session=session,
                metadata_tables=metadata,
                exp_ids=list(overlapping_expid)
            )
            recreated_metadata = create_metadata_tables(analysis_infos)
            metadata = modified_metadata.join(recreated_metadata)

        write_db_files(db_dir, analysis_infos, metadata, session)




def create_database(
        outdir,
        analysis_infos:list[AnalysisInfo],
        ask_before_deleting=True,
        run_server=False,
        port=8050,
        max_analyses_for_testing:int=None,
) -> None:
    outdir = Path(outdir)
    import glob

    outdir.mkdir(exist_ok=True, parents=True)

    old_files = [fn for fn in glob.glob(str(outdir) + '/*') if (Path(fn).name in DB_FILES)]

    if old_files:
        if ask_before_deleting:

            print(f"Output dir old database files, these will be deleted before continuing.")
            input("Press enter to continue. Ctrl+C to cancel")
        for f in old_files:
            os.remove(f)

    engine_url = get_db_url(outdir)
    engine = create_engine_with_schema(
        engine_url
    )

    analysis_infos=analysis_infos[:max_analyses_for_testing]

    metadata = create_metadata_tables(analysis_infos)

    with Session(engine) as session:
        write_db_files(outdir,  analysis_infos, metadata, session)

    if run_server:
        from crispr_screen_viewer.launch import init_app
        app = init_app(
            str(outdir),
            engine_url,
            debug_messages=True,
        )

        app.run_server(debug=False, host='0.0.0.0', port=port, )


def run_test_server(port=8050):
    datadir = get_resource_path('data/test_db')
    from crispr_screen_viewer.launch import init_app
    app = init_app(
        datadir,
        debug_messages=False
    )
    app.run_server(debug=False, host='0.0.0.0', port=port, )

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

    create_database(outd, analysis_infos[:],  run_server=False)


def parse_cli_args(args):
    from argparse import ArgumentParser

    parser = ArgumentParser(
        'CrSV database updater',
        description='Write database files for launching the CRISPR Screen Explorer.'
    )
    parser.add_argument(
        'indirs', metavar='DIRECTORIES',
        nargs='+',
        help='Input directories that contain data in Exorcise structure'
    )
    parser.add_argument(
        '--dbdir', '-d', metavar='DIRECTORY',
        help='Location to which output files will be written',
        required=True
    )
    parser.add_argument(
        '--new-db', '-n',
        action='store_true',
        help='Create a new database in the output directory. If database files are present '
             'you will be asked before they are replaced, unless -f is also set.'
    )
    parser.add_argument(
        '--update-existing', '-u',
        action='store_true',
        help='Experiment IDs already present in the database will be updated with new data. By default they are skipped.'
    )
    parser.add_argument(
        '--force-overwrite', '-f',
        action='store_true',
        help='When creating new database, overwrite existing database files without asking.'
    )
    parser.add_argument(
        '--verbosity', '-v', metavar='N',
        type=int, default=1,
        help='Set verbosity level: 0=warning, 1=info (default), 2=debug'
    )

    @dataclass
    class CLIArgs:
        indirs: list[str]
        dbdir: str
        update_existing: bool
        new_db: bool
        force_overwrite: bool
        verbosity: int

        def verbosity_str(self) -> str:
            return ['WARNING', 'INFO', 'DEBUG'][self.verbosity]

    args = CLIArgs(**vars(parser.parse_args(args)))

    for d in args.indirs:
        if not os.path.isdir(d):
            raise ValueError(f"{d} is not a directory")
    analysis_infos = get_paths_exorcise_structure_v1(args.indirs)

    set_loguru_level(logger, args.verbosity_str())

    if args.new_db:
        create_database(
            args.dbdir,
            analysis_infos,
            ask_before_deleting=not args.force_overwrite,
        )
        return 0

    # else we're adding to existing DB
    update_database(
        args.dbdir,
        analysis_infos,
        refseq=False,
        update_experiments=args.update_existing,
    )

def test_cli_args():
    data_path = Path(get_resource_path('tests/exorcise_style/'))
    paths = [data_path / d for d in os.listdir() if os.path.isdir(d)]
    logger.debug(paths)
    parse_cli_args(
        *paths,


    )

if __name__ == '__main__':
    import sys
    parse_cli_args(sys.argv[1:])



