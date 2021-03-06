"""Tests for the annotation CLI."""
from io import StringIO

from fake_open import md5_check

from tssv_extras import annotate


class TestAnnotation(object):
    """Test the annotation CLI."""
    def setup(self):
        self._input = open('data/m1_newalleles.csv')
        self._output = StringIO()

    def test_annotate(self):
        annotate.annotate(self._input, 'TCCGTCCCATGCATGC', self._output, 0)
        assert md5_check(
            self._output.getvalue(), '596241e84c6d20a4155034236c234c51')
