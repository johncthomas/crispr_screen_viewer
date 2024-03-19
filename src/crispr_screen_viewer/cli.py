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
        database: Create or manage databases
        launch: Launch the viewer
    
    See "crispr-screen-viewer COMMAND --help" for more information on each command.
    ''')

        exit(0)

    command = args[0]
    cmd_args = args[1:]

    # **************
    # ** Reminder **
    # **************
    # If you add a new command, don't forget to update the help string.

    if command == 'database':
        from crispr_screen_viewer.update_database import parse_cli_args
        parse_cli_args(cmd_args)
    if command == 'launch':
        from crispr_screen_viewer.launch import from_cli
        from_cli(cmd_args)

if __name__ == '__main__':
    run()