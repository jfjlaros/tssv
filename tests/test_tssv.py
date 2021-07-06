"""Tests for the annotation CLI."""
from io import StringIO

from tssv.tssv import parse_library
from tssv.tssv import make_json_report


class TestTSSV(object):
    """Test the annotation CLI."""
    def setup(self):
        self._output = StringIO()

    def test_parse_library(self):
        library = parse_library(open('data/library.csv'), 0.1)
        assert len(library) == 4
        assert library['m3']['flanks'] == ['TTATTATCTCTC', 'CTATCGAGAGAGAT']
        assert library['m3']['reg_exp'].pattern == '^(TTTAT){1,1}(GGGA){0,1}$'
        assert library['m3']['thresholds'] == [2, 2]

    def test_parse_library_without_pattern(self):
        library = parse_library(open('data/library_lite.csv'), 0)
        assert library['m3']['reg_exp'].pattern == '(?!x)x'

    def test_parse_library_with_mismatches(self):
        library = parse_library(open('data/library.csv'), 0.1, 1)
        assert library['m3']['thresholds'] == [1, 1]

    def test_clean_allele(self):
        raw = ['A', 'AA', 'A(1.0)', 'AA(2.4)']
        expected = ['A', 'AA', 'A', 'AA']

        out = [make_json_report.clean_allele(x) for x in raw]

        assert out == expected
