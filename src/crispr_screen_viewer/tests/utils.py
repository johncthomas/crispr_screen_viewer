import os
from glob import glob
from pathlib import Path

from loguru import logger
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from crispr_screen_viewer.functions_etc import get_resource_path
from crispr_screen_viewer.update_database import get_paths_exorcise_structure_v1, create_database, \
    create_engine_with_schema
from crispr_screen_viewer.dataset import get_db_url
from crispr_screen_viewer.update_gene_table import locate_exorcise_files, tabulate_genes_from_exorcise, upsert_genes


test_db_dir = Path(get_resource_path('tests/test_data/test_db'))
test_db_url =  get_db_url(test_db_dir)
def delete_test_files():
    fn = glob(str(test_db_dir/'*'))
    for f in fn:
        os.remove(f)

def create_test_database(**kwargs):
    """Create a database from ./tests/test_data, output to ./data/test_db
    by default. **kwargs passed to create_database()"""

    g = glob(f'{get_resource_path("tests/test_data/exorcise_style")}/*')
    logger.debug(g)
    infos = get_paths_exorcise_structure_v1(g)
    logger.debug(infos)

    default_kwargs = dict(
        outdir=test_db_dir,
        analysis_infos=infos,
        ask_before_deleting=False)
    create_database(**default_kwargs | kwargs)


def run_test_server(port=8050):
    datadir = get_resource_path('tests/test_data/test_db')
    from crispr_screen_viewer.launch import init_app
    app = init_app(
        datadir,
        debug_messages=False
    )
    app.run_server(debug=False, host='0.0.0.0', port=port, )


def test_update_gene_table():
    engine = create_engine_with_schema(test_db_url)
    files = locate_exorcise_files('/Users/thomas03/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/ddrcs/test_data/exorcise_libs')
    table = tabulate_genes_from_exorcise(
        files,
        '/Users/thomas03/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/ddrcs/hgnc_complete_set.txt'
    )
    with Session(engine) as session:
        upsert_genes(table.to_dict(orient='records'), session)




if __name__ == '__main__':
    delete_test_files()
    test_update_gene_table()
    create_test_database()
    run_test_server(8054)


