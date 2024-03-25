from glob import glob

from loguru import logger

from crispr_screen_viewer.functions_etc import get_resource_path
from crispr_screen_viewer.update_database import get_paths_exorcise_structure_v1, create_database


def create_test_database(**kwargs):
    """Create a database from ./tests/test_data, output to ./data/test_db
    by default. **kwargs passed to create_database()"""
    outdir = get_resource_path('tests/test_data/test_db')
    g = glob(f'{get_resource_path("tests/test_data/exorcise_style")}/*')
    logger.debug(g)
    infos = get_paths_exorcise_structure_v1(g)
    logger.debug(infos)

    default_kwargs = dict(
        outdir=outdir,
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


if __name__ == '__main__':
    create_test_database()
    run_test_server(8054)