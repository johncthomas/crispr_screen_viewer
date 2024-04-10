#!/usr/bin/env python
import sys
from argparse import Namespace, ArgumentParser
from dataclasses import dataclass
from typing import Collection, Sequence

EXORCISEFN = 'exorcise.tsv.gz'

# NOTE DO NOT import anything from crispr_screen_viewer at the top level,
#   the aim here is to keep help options for the command line snappy,
#   and it's also likely to cause circular import issues.

@dataclass
class ArgsUpdateDB:
    details_xlsx: list[str]
    out_dir: str
    results_dir: str
    counts_dir: str
    filename_prefix: str
    new_db: bool
    update_existing: bool
    force_overwrite: bool
    verbosity: int

    def verbosity_str(self) -> str:
        return ['WARNING', 'INFO', 'DEBUG'][self.verbosity]

def parse_update_database_args(args:Sequence[str]) -> ArgsUpdateDB:
    from argparse import ArgumentParser

    parser = ArgumentParser(
        'CrSV database updater',
        description='Write database files for launching the CRISPR Screen Explorer.'
    )
    parser.add_argument(
        'details_xlsx', nargs='+', metavar='DETAILS_XLSX',
    )
    parser.add_argument(
        '--out-dir', '-o', metavar='DIRECTORY',
        help='Location to which output files will be written',
        required=True
    )
    parser.add_argument(
        '--results-dir', '-r', metavar='DIRECTORY',
        help='Directory containing output results',
        default='results',
    )
    parser.add_argument(
        '--counts-dir', '-c', metavar='DIRECTORY',
        help='Directory containing counts files',
        default='counts'
    )
    parser.add_argument(
        '--filename-prefix', '-p',
        help='Default "result". String prepended to input files.',
        default='result',
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
        help='Set verbosity level: 0 = warnings only, 1 = info (default), 2 = debug'
    )

    args = ArgsUpdateDB(**vars(parser.parse_args(args)))

    return args


def parse_gene_table_args(args:Sequence[str]) -> Namespace:
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
        '--hgnc-table', '-H',
        help='Table from https://www.genenames.org/download/archive/, "tab separated hgnc_complete_set file"',
        required=True,
    )

    args =  parser.parse_args(args)
    return args

def parse_launch_args(args:Sequence[str]) -> Namespace:
    launcher_parser = ArgumentParser(add_help=False)

    launcher_parser.add_argument(
        '-p', '--port', metavar='PORT',
        help='Port used to serve the app',
        required=True,
    )

    launcher_parser.add_argument(
        '-d', '--data-path',
        dest='data_path',
        help="Name of the directory that contains the screens' data.",
        required=True, default=None,
    )
    launcher_parser.add_argument(
        '--app-debug', action='store_true',
        help='Launch Dash app in debug mode, tends to break things, but allows you to look at the'
             ' callback graph and JS debug messages.'
    )
    launcher_parser.add_argument(
        '--debug-messages', action='store_true', default=False,
        help='Set log level to debug – print messages describing the internal state of the app. '
             'Also hide Werkzeug messages'
    )
    launcher_parser.add_argument(
        '--public-version', action='store_true',
        help="Don't hide the data-source and analysis-type selectors."
             " In the future analysis-type might have its own option."
    )
    launcher_parser.add_argument(
        '--url-pathname', default="/",
        help="URL base pathname. Needs to end in a /."
    )
    launcher_parser.add_argument(
        '-u', '--database-url',
        help='A SQL database URL, as understood by SQL Alchemy. '
             '\nSee: https://docs.sqlalchemy.org/en/20/core/engines.html#database-urls.'
             '\nE.g.: sqlite:////home/user/some/dir/databasefile.db',
        default=None, required=False
    )

    parser = ArgumentParser(parents=[launcher_parser],
                            description="Dash app for exploring CRISPR screen data.",
                            add_help=True,)

    args = parser.parse_args(args)
    return args


def run():
    args = sys.argv[1:]

    # I don't like how argparse subparsers work, handling it myself.

    # print help
    if (not args) or (args[0] in ('-h', '--help')):
        print('''
    CRISPR Screen Viewer
    
    Use "crispr-screen-viewer COMMAND".
    
    Commands:
        database: Create or update results database
        genes: Update genes from exorcise directory
        launch: Launch the viewer
        remove: Remove experiments from a database.
        test: Run test server, default port 8054.
    
    See "crispr-screen-viewer COMMAND --help" for more information on each command.
    ''')

        exit(0)

    command = args[0]
    cmd_args = args[1:]
    del args # avoid mixups

    # **************
    # ** Reminder **
    # **************
    # Before you add a new command, update the help string.

    if command == 'database':
        parse_update_database_args(cmd_args)

        from crispr_screen_viewer.update_database import run_from_cli
        run_from_cli(cmd_args)

    elif command == 'genes':
        parse_gene_table_args(cmd_args) # print help only

        from crispr_screen_viewer.update_gene_table import run_from_cli
        run_from_cli(cmd_args)

    elif command == 'launch':
        parse_launch_args(cmd_args) # print help only

        from crispr_screen_viewer.launch import from_cli
        from_cli(cmd_args)

    elif command == 'remove':
        if (not cmd_args) or (cmd_args[0] in ('-h', '--help')):
            print('''
    crispr-screen-viewer DATABASE_DIRECTORY EXP_ID [EXP_ID]
    
    Provide a path to the database you wish to alter, then a list of Experiment IDs to remove. 
    If you don't know the IDs, they can be found in experiments_metadata.csv.gz.
''')
        else:
            from crispr_screen_viewer.update_database import remove_experiments_from_db
            remove_experiments_from_db(cmd_args[0], cmd_args[1:])

    elif command == 'test':
        try:
            port = int(cmd_args[0])
        except:
            port = 8054
        from crispr_screen_viewer.tests.utilities import create_test_database, run_test_server
        from loguru import logger
        from crispr_screen_viewer.functions_etc import set_loguru_level
        set_loguru_level(logger, 'DEBUG')
        create_test_database()
        run_test_server(port)

    else:
        print(f'Command {command} not recognised.')

if __name__ == '__main__':
    run()

