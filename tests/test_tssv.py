"""Tests for the annotation CLI."""
from hashlib import md5
from io import StringIO

from fake_open import md5_check

from tssv.tssv import parse_library


class TestTSSV(object):
    """Test the annotation CLI."""
    def setup(self):
        self._output = StringIO()

    def test_parse_library_with_pattern(self):
        library = parse_library(open('data/library.csv'), 0)
        assert len(library) == 4
        assert library['m3']['flanks'] == ['TTATTATCTCTC', 'CTATCGAGAGAGAT']
        assert library['m3']['reg_exp'].pattern == '^(TTTAT){1,1}(GGGA){0,1}$'

    def test_parse_library_without_pattern(self):
        library = parse_library(open('data/library_lite.csv'), 0)
        assert library['m3']['reg_exp'].pattern == '(?!x)x'
