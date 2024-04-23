#!/usr/bin/python3
import logging
import coloredlogs

# try to handle matplotlib's annoying logs about font
logging.getLogger('matplotlib').setLevel(logging.WARNING)
# project name
name = 'OGU'
# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
log = logging.getLogger(name)
coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
# store global values here
global_dict = {}