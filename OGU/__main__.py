#!/usr/bin/python3

from sys import argv
from .cli import cli_main
from .utils import get_all_third_party


def main():
    if argv[-1] in ('init', '-init', '--init'):
        get_all_third_party()
    elif argv[-1] in ('-h', '--help'):
        cli_main()
    elif len(argv) > 1:
        cli_main()
    else:
        try:
            from OGU.ui import ui_main
            ui_main()
        except Exception:
            cli_main()


if __name__ == '__main__':
    main()
