#!/usr/bin/python3

import re
import logging

# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
log = logging.getLogger('barcodefinder')
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass


def plastid_rename():
    """
    Use name database.
    """
    pass


def safe_average(x):
    """
    Safe average.
    """
    if len(x) == 0:
        return 0
    else:
        return sum(x) / len(x)


def safe_path(old):
    """
    Remove illegal character in file path or name.
    """
    return re.sub(r'\W', '_', old)
