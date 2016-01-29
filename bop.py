#!/usr/local/bin/python3

#This program is to add filename into fasta sequence's id, and merge same gene into one file.

from glob import glob
from Bio import SeqIO
import argparse
import re

parser = argparse.ArgumentParser(description='This program will add filename
into fasta sequence's id, and merge same gene files into single one.')
