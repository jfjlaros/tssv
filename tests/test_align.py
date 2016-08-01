"""
Tests for the alignment modules.
"""
from Bio import Seq

from tssv import align_pair, sg_align


class TestAlign(object):
    def setup(self):
        self._reference = 'GACTGTCGTGGGCTCTTACGCACATATTATAACTTTCATAAGTTTGTCAG'
        self._reference_rc = Seq.reverse_complement(self._reference)

    def test_align_pair_perfect(self):
        result = align_pair.align_pair(
            self._reference, self._reference_rc, ('TTATGAAAGT', 'CGTAAGAGC'))
        assert result == ((3, 34), (0, 11))

    def test_align_pair_subst(self):
        result = align_pair.align_pair(
            self._reference, self._reference_rc, ('TTATGTAAGT', 'CGTATGAGC'))
        assert result == ((3, 34), (1, 11))

    def test_align_pair_del(self):
        result = align_pair.align_pair(
            self._reference, self._reference_rc, ('TTATGAAGT', 'CGTAGAGC'))
        assert result == ((2, 34), (1, 11))

    def test_align_pair_ins(self):
        result = align_pair.align_pair(
            self._reference, self._reference_rc, ('TTATGACAAGT', 'CGTAACGAGC'))
        assert result == ((3, 43), (1, 11))

    def test_align_pair_del_punish_indel(self):
        result = align_pair.align_pair(
            self._reference, self._reference_rc, ('TTATGAAGT', 'CGTAGAGC'), 3)
        assert result == ((3, 35), (3, 11))

    def test_align_pair_ins_punish_indel(self):
        result = align_pair.align_pair(
            self._reference, self._reference_rc, ('TTATGACAAGT', 'CGTAACGAGC'),
            3)
        assert result == ((4, 34), (3, 11))
