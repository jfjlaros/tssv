"""
TSSV: Targeted characterisation of short structural variation.


The library file consists of four tab-separated columns:
- name of the marker pair
- marker 1
- marker 2
- space-separated description of the expected pattern

If the -d option is used, a folder will be created containing the library
table and per marker a subfolder containing the new alleles and split FASTA
files.

Copyright (c) 2014 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2014 Jeroen F.J. Laros <j.f.j.laros@lumc.nl>
Copyright (c) 2012 Jaap W.F. van der Heijden

Licensed under the MIT license, see the LICENSE file.
"""

import argparse
import os

__version_info__ = ('0', '2', '4')

__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Jeroen F.J. Laros'
__contact__ = 'j.f.j.laros@lumc.nl'
__homepage__ = 'https://humgenprojects.lumc.nl/trac/tssv'

usage = __doc__.split("\n\n\n")

class ProtectedFileType(argparse.FileType):
    def __call__(self, string):
        if 'w' in self._mode and os.path.exists(string):
            raise IOError('failed to create "{}": file exists.'.format(string))
        return super(ProtectedFileType, self).__call__(string)

def doc_split(func):
    return func.__doc__.split("\n\n")[0]

def version(name):
    return "%s version %s\n\nAuthor   : %s <%s>\nHomepage : %s" % (name,
        __version__, __author__, __contact__, __homepage__)
