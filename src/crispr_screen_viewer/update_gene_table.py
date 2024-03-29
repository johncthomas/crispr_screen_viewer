from pathlib import Path
import os
import pandas as pd

from sqlalchemy.orm import Session

from crispr_screen_viewer.database import GeneTable
from crispr_screen_viewer.update_database import upsert_records, create_engine_with_schema

from loguru import logger

EXORCISEFN = 'exorcise.tsv.gz'

def locate_exorcise_files(d:str|Path) -> list[Path]:
    d = Path(d)

    if os.path.isfile(d):
        return [d]

    if os.path.isfile(d/EXORCISEFN):
        return [d/EXORCISEFN]

    paths = []

    exor_libraries = set(os.listdir(d))

    for lib in exor_libraries:
        if lib.startswith('.'):
            continue
        fn = d / lib / 'exorcise.tsv.gz'
        if not os.path.isfile(fn):
            continue
        paths.append(fn)

    return paths

def tabulate_genes_from_exorcise(exo_fn:list[str|Path], hgnc_fn: str|Path):
    extern_id_col = 'inherit.externalId'
    table_rows = []

    ids_added = set()

    # iterate through file names
    for fn in exo_fn:
        logger.info(f"Parsing from {fn}")
        tab = pd.read_csv(fn, sep='\t', dtype=str)
        isctrl = tab.exo_symbol.str.contains('Non-targ') | tab.exo_symbol.str.contains('Cutting')
        tab = tab.loc[~isctrl]

        if tab.shape[0] == 0:
            print(f'{fn}, All the genes are non-targeting or file is empty. Moving on...')
            continue

        # determine organism
        ishuman = tab[extern_id_col].str.startswith('HGNC:').any()
        ismouse = tab[extern_id_col].str.startswith('MGI:').any()

        if ishuman == ismouse:
            raise ValueError(f"Table contains neither (or both) HGNC or MGI IDs, hgnc={ishuman} mouse={ismouse}")

        if ishuman:
            organism = 'Human'
        else:
            organism = 'Mouse'

        tab:pd.DataFrame = tab.drop_duplicates('exo_symbol').loc[:, ['exo_symbol', extern_id_col]]
        tab.columns = ['symbol', 'official_id']
        tab.loc[:, 'organism'] = organism
        tab.loc[:, 'id'] = tab.symbol

        tab = tab.loc[~tab.symbol.isin(ids_added)]
        ids_added.update(tab.symbol.values)

        tab.official_id = tab.official_id.str.replace('X', '')

        table_rows.append(tab)

        logger.info(f'{tab.shape[0]} new genes added.')

    gene_table = pd.concat(table_rows, ignore_index=True, ).astype(str)



    assert not gene_table.symbol.duplicated().any()

    gene_table = get_symbol_ids_from_hgnc(gene_table, hgnc_fn)

    ismgi = gene_table.official_id.str.startswith('MGI:')
    gene_table.loc[ismgi, 'symbol_with_ids'] = gene_table.loc[ismgi, 'official_id'].map(lambda s: f"({s})")

    return gene_table


def get_symbol_ids_from_hgnc(gene_table:pd.DataFrame, hgnc_fn:str|Path) -> pd.DataFrame:
    """Add symbol_with_ids column for every row where row.official_id contains 'HGNC:".

    Modifies in place, and returns gene_table."""
    hgnc = pd.read_csv(hgnc_fn, dtype=str, sep='\t')
    hgnc.set_index('hgnc_id', inplace=True, drop=False)
    hgnc.fillna('', inplace=True)

    ishgnc = gene_table.official_id.str.startswith('HGNC:')
    hids = gene_table.loc[ishgnc, 'official_id'].values

    hgnc = hgnc.loc[hids]

    hgnc.entrez_id = hgnc.entrez_id.map(lambda i: f"NCBI:{i}")

    def hgnc_row_to_long(row):
        prev_id = ', '.join([x for x in (row.prev_symbol, row.hgnc_id, row.entrez_id) if x])
        if prev_id:
            prev_id = f"  ({prev_id})"
        return row.symbol + prev_id

    gene_table.loc[ishgnc.values, 'symbol_with_ids'] = hgnc.apply(hgnc_row_to_long, axis=1).values
    return gene_table


def upsert_genes(records:list[dict], session:Session):
    """Update GeneTable. Adds new records when ID not found, updates otherwise"""
    upsert_records(
        records,
        session,
        GeneTable,
        primary_key='id'
    )

def run_from_commandline(args):
    from argparse import ArgumentParser

    parser = ArgumentParser(
        'CrSV gene table updater',
        description='Update gene table `symbol_with_ids` column, used for disambiguation when searching for a gene.'
    )

    parser.add_argument(
        '--exorcise-tables', '-x',
        help=f'Location of output of Exorcise. Either a single file, a directory containing "{EXORCISEFN}" '
             f'or a directory containing directories at least some of which contain a file called "{EXORCISEFN}',
        required=True,
    )

    parser.add_argument(
        '--database-path', '-d',
        required=True,
    )

    parser.add_argument(
        '--hgnc-table', '-h',
        help='Table from https://www.genenames.org/download/archive/, "tab separated hgnc_complete_set file"',
        required=True,
    )



    return parser.parse_args(args)


def update_gene_table_hgnc(exorcise_location:str|Path,
                           hgnc_path:str|Path,
                           database_path:str|Path):

    from crispr_screen_viewer.dataset import get_db_url
    db_url = get_db_url(database_path)
    engine = create_engine_with_schema(db_url, echo=False)
    files = locate_exorcise_files(exorcise_location)
    table = tabulate_genes_from_exorcise(files, hgnc_path)
    with Session(engine) as session:
        upsert_genes(table.to_dict(orient='records'), session)
        session.commit()


if __name__ == '__main__':

    import sys
    args = run_from_commandline(sys.argv[1:])
    update_gene_table_hgnc(args.exorcise_tables, args.hgnc_table, args.database_path)