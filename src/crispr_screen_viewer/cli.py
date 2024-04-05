#!/usr/bin/env python
import sys

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
        gene_db: Update genes from exorcise directory
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
        from crispr_screen_viewer.update_database import run_from_cli
        run_from_cli(cmd_args)

    elif command == 'gene_db':
        from crispr_screen_viewer.update_gene_table import run_from_cli
        run_from_cli(cmd_args)

    elif command == 'launch':
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
        create_test_database()
        run_test_server(port)

    else:
        print(f'Command {command} not recognised.')

if __name__ == '__main__':
    run()