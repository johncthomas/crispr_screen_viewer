#!/usr/bin/env python
from dataclasses import dataclass

import requests, json
import os.path, typing
from pathlib import Path

from typing import Literal, Type

import pandas as pd
import sqlalchemy
from sqlalchemy import Engine, create_engine
from sqlalchemy.orm import Session

from crispr_screen_viewer.database import *

from crispr_screen_viewer import functions_etc
from crispr_screen_viewer.functions_etc import (
    df_rename_columns,
    normalise_text,
    load_stats_csv,
    set_loguru_level,
    is_temp_file,
    is_nt,
    maybe_its_gz,
    timepoint_labels
)
from crispr_screen_viewer.dataset import (
    ANALYSESTYPES,
    AnalysisType,
    MetadataTables,
    DB_FILES,
    get_db_url
)

from crispr_tools.data_classes import AnalysisWorkbook

analysis_tabulate:dict[str, typing.Callable]

import numpy as np

from loguru import logger


def doi_to_link(doi:str):
    """Return string formated as a Markdown link to doi.org/{doi}.

    Args:
        doi: can be formatted as full URl or like 1.2.3/whatever"""
    if pd.isna(doi) or (not doi):
        return ''
    # if the string is a proper URL, strip it down to the DOI
    doi = doi.replace('https', '').replace('http', '').replace('://', '').replace('doi.org/', '')
    return f"[{doi}](https://doi.org/{doi})"

def get_treatment_str(samp_deets:pd.DataFrame, ctrl:str, treat:str):
    """Return a string describing the treatment performed between the
    control and test samples.

    Args:
        samp_deets: dataframe, rows indexed by samp IDs, giving sample
            metadata.
        ctrl, treat: Sample IDs of a comparison.
    """

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
    """For creating a new database"""
    engine = create_engine(destination, echo=echo, connect_args={'timeout':4000})

    # create tables in the actual SQL database
    TableBase.metadata.create_all(engine)

    return engine


def insert_records(table:Type[TableBase], records:list[dict], session:Session):
    """Insert rows into given database table.
    Record dict.keys must be present in table columns.
    DF to records: DataFrame.to_dict(orient='records')

    Column values that are pd.nan are skipped over.

    Fails if they've already been added (i.e. primary key value exists)."""

    # filter the records so key:value pairs of null value are removed
    session.add_all(
        [table(**{k: v for k, v in kw.items() if not pd.isna(v)})
         for kw in records]
    )

def upsert_records(
        records:list[dict],
        session:Session,
        table:Type[TableBase],
        primary_key='id',
) -> None:
    """Update existing records, or add if primary_id not found in the table."""
    logger.debug("First 5 records:\n" + f"{records[:5]}")
    for record in records:
        try:
            # Try to fetch the existing record by 'id'
            existing_record = session.query(
                table
            ).filter(
                getattr(table, primary_key) == record[primary_key]
            ).one()

            # If found, update the existing record
            for key, value in record.items():
                setattr(existing_record, key, value)

        except sqlalchemy.exc.NoResultFound:
            # If not found, create a new record
            new_record = table(**record)
            session.add(new_record)


@dataclass
class AnalysisInfo:
    """Collection of information used to create tables,
    analysis_workbook and paths."""
    experiment_id:str
    analysis_workbook: AnalysisWorkbook
    counts_path: str
    results_paths: dict[AnalysisType, Path | str]


def get_paths_branched_structure_v1(
        source_dirs:list[str|Path],
        details_dir='./det',
        results_dir='./res',
        count_dir='./cts',
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
        fns = list(filter(
            lambda x: not (x.startswith('.') or x.startswith('~$')),
            fns
        ))
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

        count_filepath = get_fn_or_crash( basepath/count_dir)

        d = AnalysisInfo(
            experiment_id=xpid,
            analysis_workbook=analysiswb,
            results_paths=resultspaths,
            counts_path=count_filepath
        )
        data.append(d)

    return data

def get_paths(
        details_xlsx:list[str],
        results_dir:str|Path,
        count_dir:str|Path,
        analysis_filename_prefix:str='result'
) -> list[AnalysisInfo]:
    """Get paths to files, assuming flat structure where files of same type of
    the same type are placed together in directories.
    """
    infos = []
    for fn in details_xlsx:
        if is_temp_file(fn):
            continue

        analysiswb = AnalysisWorkbook(fn)

        xpid = analysiswb.experiment_details['Experiment name']

        cfns = analysiswb.analylses['Counts file'].dropna().unique()
        if len(cfns) != 1:
            raise RuntimeError(f"Only a single counts file is actually supported. Error in workbook {fn}, "
                               f"multiple (or zero) counts files found: {cfns}")

        # get path for counts file
        counts_path = maybe_its_gz(
            os.path.join(count_dir, cfns[0])
        )
        assert os.path.isfile(counts_path), f"Counts file for experiment {xpid} not found, {counts_path}"

        # get paths for results tables
        resultspaths = {}

        for ans in ANALYSESTYPES:
            tabfn = maybe_its_gz(os.path.join(
                results_dir, xpid, 'tables', f"{analysis_filename_prefix}.{ans.name}_table.csv"
            ))
            resultspaths[ans] = tabfn
            
        infos.append(
            AnalysisInfo(
                experiment_id=xpid,
                analysis_workbook=analysiswb,
                counts_path=counts_path,
                results_paths=resultspaths
            )
        )

    return infos


def tabulate_experiments_metadata(experiment_details:list[pd.DataFrame]) \
        -> pd.DataFrame:
    """Takes information from the "Experiment details" sheet, returns DF formated
    to be added to the database."""

    column_renamer = {
        "Experiment name":"Experiment ID",
        'Analysis name': 'Experiment ID', # old label for experiment name
        'Experiment description (a few sentances)': 'Experiment description',
        'Experiment description (a few sentences)':'Experiment description',
        # Reference information used to be in a field called Citation.
        #   I leave this here for backwards compat, but it could cause issues
        #   if we start specifying actual citation strings in experiment details.
        "Citation":"Reference",
        'Date screen completed (yyyy-mm-dd)':'Date',
        'Date screen completed':'Date',
        'Date published':'Date',
    }

    experiment_details = [xpmet.copy().drop(np.nan, errors='ignore') for xpmet in experiment_details]

    logger.debug(experiment_details)
    for expmet in experiment_details:
        df_rename_columns(expmet, column_renamer, inplace=True, axis='index')

    experiment_details_table = pd.DataFrame(experiment_details)
    experiment_details_table = experiment_details_table.reset_index(drop=True)

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
        s = '\n\t'.join(bad_cite)
        logger.warning("The following experiments have bad citation info:\n"+s)

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

    experiment_details_table.DOI = experiment_details_table.DOI.apply(doi_to_link).values

    experiment_details_table = experiment_details_table.infer_objects()

    return experiment_details_table


def tabulate_comparisons(analysis_wb:AnalysisWorkbook) -> pd.DataFrame:
    """Create DF of comparison information for a single experiment"""
    logger.debug(f"Tabulating comps of {analysis_wb.experiment_details['Experiment name']}")
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

        contrast = comprow['Contrast']
        if len(contrast) == 0:
            contrast = get_treatment_str(
                analysis_wb.wb['Sample details'],
                ctrl, treat
            )
        comparison_info['Contrast'] = contrast

        treat_row = analysis_wb.samples.loc[treat]
        comparison_info['Treatment'] = treat_row['Treatment']
        comparison_info['KO'] = treat_row['KO']
        comparison_info['Dose'] = treat_row['Dose']
        comparison_info['Growth inhibition %'] = treat_row['Growth inhibition %']
        comparison_info['Days grown'] = treat_row['Days grown']
        comparison_info['Cell line'] = treat_row['Cell line']
        comparison_info['Notes'] = treat_row['Notes']

        ctrl_row = analysis_wb.samples.loc[ctrl]
        comparison_info['ControlTreatment'] = ctrl_row['Treatment']
        comparison_info['ControlKO'] = ctrl_row['KO']
        comparison_info['ControlDose'] = ctrl_row['Dose']
        comparison_info['ControlGrowth inhibition %'] = ctrl_row['Growth inhibition %']
        comparison_info['ControlDays grown'] = ctrl_row['Days grown']
        comparison_info['ControlCell line'] = ctrl_row['Cell line']
        comparison_info['ControlNotes'] = ctrl_row['Notes']

        exp_id = analysis_wb.expd['experiment_id']
        comparison_info['Experiment ID'] = exp_id
        comparison_info['Library'] = analysis_wb.experiment_details['Library']
        comparison_info['Comparison ID'] = f"{exp_id}.{ctrl}-{treat}"

        # format strings
        for k in ('Dose', 'ControlDose'):
            if not pd.isna(comparison_info[k]):
                comparison_info[k] = str(comparison_info[k]).replace('uM', 'μM')
        for k in ('KO', 'ControlKO'):        
            if pd.isna(comparison_info[k]) or (comparison_info[k] == '') or (comparison_info[k] == 'WT'):
                comparison_info[k] = 'Wildtype'

        comparison_info['Timepoint'] = comprow['Group'].split('_')[0]

        comparisons_metadata.append(comparison_info)

    comparisons = pd.DataFrame(comparisons_metadata)

    # this is sometimes put in wrong...
    m = comparisons['Timepoint'] == 'endpoint'
    comparisons.loc[m, 'Timepoint'] = 'endpoints'

    # fill in blank treatments
    comparisons.loc[comparisons.Treatment.isna(), 'Treatment'] = 'No treatment'
    # these cols could be blank and aren't essential to have values
    for col in ['Cell line', 'Library', 'Source']:
        if col not in comparisons.columns:
            comparisons.loc[:, col] = 'Unspecified'
        comparisons.loc[comparisons[col].isna(), col] = 'Unspecified'

    # Replace e.g. "endpoints" with "Matched time-point"
    for old, new in timepoint_labels.items():
        comparisons.loc[comparisons.Timepoint == old, 'Timepoint'] = new

    return comparisons


def get_gene_symbols_in_db(session:Session) -> set[str]:
    """All distinct gene symbols in the database."""
    gns = set([g[0] for g in session.query(GeneTable.symbol).distinct().all()])
    return gns


def add_genes_from_symbols(
        symbols: list[str],
        organism: Literal['Human'] | Literal['Mouse'],
        session:Session,
):
    """Adds records for genes in results that are not currently in the table. """

    symbols_to_add = set(symbols).difference(get_gene_symbols_in_db(session))

    logger.debug('Num ID-less genes being added: '+str(len(symbols_to_add)))
    empty_records = [dict(id=s, symbol=s, organism=organism) for s in symbols_to_add]
    insert_records(GeneTable, empty_records, session)


def tabulate_statistics(info:AnalysisInfo) -> pd.DataFrame:
    """Create table of analysis statistics (normZ score, FDR etc) for given
    experiment."""
    experiment_id = info.experiment_id
    tables = []
    for analysis_type, fn in info.results_paths.items():
        logger.debug(f"tabulating from {fn}")
        try:
            stats_table = load_stats_csv(fn)
            stats_table.index.name = 'gene_id'
            for cmp in stats_table.columns.levels[0]:
                table = stats_table[cmp].reset_index()
                df_rename_columns(table, {'lfc': 'score', 'normZ': 'score',
                                          'fdr_log10': 'fdr10'}, inplace=True)
                table.loc[:, 'comparison_id'] = f"{experiment_id}.{cmp}"
                table.loc[:, 'analysis_type_id'] = analysis_type.id
                table.loc[:, 'experiment_id'] = experiment_id
                if analysis_type.id == 1 or analysis_type.id == 2:
                    table = table.loc[:, ['gene_id', 'score', 'fdr', 'fdr10', 'pos_p', 'neg_p',
                                          'comparison_id', 'analysis_type_id', 'experiment_id']]
                # special case when it's chronos
                elif analysis_type.id == 3:
                    table = table.loc[:, ['gene_id', 'chronos_score',
                                          'comparison_id', 'analysis_type_id', 'experiment_id']]
                    table['score'] = table['chronos_score']
                    table['fdr'] = table['chronos_score']
                    table['fdr10'] = table['chronos_score']
                    table['pos_p'] = table['chronos_score']
                    table['neg_p'] = table['chronos_score']
                    table = table.drop('chronos_score', axis = 1)

                # special case when it's manual
                elif analysis_type.id == 4:
                    table = table.loc[:, ['gene_id', 'score', 'pval',
                                          'comparison_id', 'analysis_type_id', 'experiment_id']]
                    table['fdr'] = table['pval']
                    table['fdr10'] = table['pval']
                    table['pos_p'] = table['pval']
                    table['neg_p'] = table['pval']
                    table = table.drop(['score', 'pval'], axis = 1)

                tables.append(table)
        # Don't fail when we encounter a ParserError on an empty Chronos table
        except (pd.errors.ParserError, AttributeError):
            continue

    return pd.concat(tables)


def add_statistics(analysesinfo:list[AnalysisInfo], session:Session):
    """Take a table, as output by crispr_pipeline, add rows to StatTable."""
    for info in analysesinfo:
        table = tabulate_statistics(info)

        add_genes_from_symbols(
            table.gene_id,
            info.analysis_workbook.experiment_details.Organism,
            session,
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

    experiments_metadata = tabulate_experiments_metadata(
            [d.analysis_workbook.experiment_details for d in analysesinfo]
    )

    comparisons_metadata.set_index('Comparison ID', drop=False, inplace=True)
    experiments_metadata.set_index('Experiment ID', drop=False, inplace=True)

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
    """Remove rows from DB and metadata tables associated with the given
    experiment IDs.

    This functions used as part of updating a database,
    use remove_experiments_fromm_db to alter DB on disk if that's the
    only thing you're doing."""

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
) -> None:
    """Remove experiment from a database on disk. Commits and writes changes."""
    db_dir = Path(db_dir)
    metadata = MetadataTables.from_files(db_dir)
    engine = create_engine(get_db_url(db_dir))
    with Session(engine) as session:
        mod_metadata = remove_experiments(metadata, experiments_to_remove, session)
        session.commit()

    mod_metadata.to_files(db_dir)


def update_database(
        db_dir:str|Path,
        analysis_infos: list[AnalysisInfo],
        update_experiments=False,
) -> None:
    """Update an existing database, inserting new experiments.

    Args:
        db_dir: Path to database directory
        analysis_infos: Info about experiments to be added
        update_experiments: If true, existing experiments (determined by experiment ID)
            will be removed from DB and replaced with new version."""
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
            modified_metadata = metadata
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
        outdir:str|Path,
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
            debug_messages=True,
        )

        app.run_server(debug=False, host='0.0.0.0', port=port, )


def run_from_cli(args: typing.Sequence[str]):
    from crispr_screen_viewer.cli import parse_update_database_args
    pargs = parse_update_database_args(args)

    for d in pargs.details_xlsx:
        if not os.path.isfile(d):
            raise ValueError(f"Input directory {d} does not exist")
    analysis_infos = get_paths(
        details_xlsx=pargs.details_xlsx,
        results_dir=pargs.results_dir,
        count_dir=pargs.counts_dir,
        analysis_filename_prefix=pargs.filename_prefix
    )

    set_loguru_level(logger, pargs.verbosity_str())

    if pargs.new_db:
        create_database(
            pargs.out_dir,
            analysis_infos,
            ask_before_deleting=not pargs.force_overwrite,
        )
        return 0

    # else we're adding to existing DB
    update_database(
        pargs.out_dir,
        analysis_infos,
        update_experiments=pargs.update_existing,
    )


if __name__ == '__main__':
    import sys
    run_from_cli(sys.argv[1:])



