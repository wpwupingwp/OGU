#!/usr/bin/python3

import os
import xlrd

def read_xls(filename):
    book=xlrd.open_workbook(filename)
#Use first sheet
    sheet=book.sheet_by_index(0)


def get_filename():

