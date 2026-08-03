import sys

import pyrematching


def cli_argv():
    pyrematching.cli(command_line_args=sys.argv[1:])
