import os
from glob import glob
from pathlib import Path

from loguru import logger
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from crispr_screen_viewer.functions_etc import get_resource_path
from crispr_screen_viewer.update_database import (
    get_paths_branched_structure_v1,
    get_paths,
    create_database,
    create_engine_with_schema
)
from crispr_screen_viewer.dataset import get_db_url
from crispr_screen_viewer.update_gene_table import (
    locate_exorcise_files,
    tabulate_genes_with_ids_from_excorcise,
    upsert_genes,
)

from loguru import logger


INFOS_branched = get_paths_branched_structure_v1(
    glob(get_resource_path('tests/test_data/branched_style/*'))
)

TEST_DB_DIR = Path(get_resource_path('tests/test_data/test_db'))
test_db_url =  get_db_url(TEST_DB_DIR)
DATA_DIR = Path(get_resource_path('tests/test_data/flat_style'))
INFOS = get_paths(
    details_xlsx=glob(f'{get_resource_path(str(DATA_DIR/'details'))}/*'),
    results_dir=DATA_DIR/'results',
    count_dir=DATA_DIR/'counts',
    analysis_filename_prefix='result',
)

def delete_test_files():
    fn = glob(str(TEST_DB_DIR / '*'))
    for f in fn:
        os.remove(f)

def create_test_database(**kwargs):
    """Create a database from ./tests/test_data, output to ./data/test_db
    by default. **kwargs passed to create_database()"""

    #g = glob(f'{get_resource_path("tests/test_data/branched_style")}/*')

    #infos = get_paths_exorcise_structure_v1(g)
    logger.debug(INFOS)

    default_kwargs = dict(
        outdir=TEST_DB_DIR,
        analysis_infos=INFOS,
        ask_before_deleting=False)
    create_database(**default_kwargs | kwargs)


def run_test_server(port=8050):
    datadir = get_resource_path('tests/test_data/test_db')
    from crispr_screen_viewer.launch import init_app
    app = init_app(
        datadir,
        debug_messages=True
    )
    app.run_server(debug=False, host='0.0.0.0', port=port, )


def test_update_gene_table():
    engine = create_engine_with_schema(test_db_url)
    files = locate_exorcise_files(get_resource_path('tests/test_data/exorcise_libs_A'))
    table = tabulate_genes_with_ids_from_excorcise(
        files,
        get_resource_path('tests/test_data/hgnc_table_cutdown.tsv.gz')
    )
    with Session(engine) as session:
        upsert_genes(
            table.to_dict(orient='records'),
            session
        )
        session.commit()


if __name__ == '__main__':
    create_test_database()
    test_update_gene_table()
    run_test_server(8054)


