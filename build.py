#!/usr/bin/python3

from distutils.core import setup
import py2exe

filename=input('filename:\n')
option={'py2exe':
         {
             'compressed':1,
             'optimize':2,
             'ascii':1,
             'includes':includes,
             'bundle_files':1
         }
         }
setup(
    options=option,
    zipfile=None,
    console=[{'script':filename}]
    )
