"""
Tests for the annotation CLI.
"""
from __future__ import (
    absolute_import, division, print_function, unicode_literals)
from future.builtins import str, zip

from hashlib import md5
from io import StringIO

from fake_open import md5_check

from tssv import annotate


class TestAnnotation(object):
    """
    Test the annotation CLI.
    """
    def setup(self):
        self._input = open('data/m1_newalleles.csv')
        self._output = StringIO()

    def test_with_dna(self):
        annotate.annotate(self._input, 'TCCGTCCCATGCATGC', self._output, 0)
        assert md5_check(
            self._output.getvalue(), '596241e84c6d20a4155034236c234c51')
