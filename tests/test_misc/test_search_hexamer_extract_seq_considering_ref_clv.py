import unittest
from unittest.mock import MagicMock

import kleat.misc.settings as S
from kleat.misc.search_hexamer import extract_seq


def test_extract_seq_for_plus_strand_clv_supported_by_suffix():
    """
             AA  <-tail of suffix contig
       ATCGAC┘   <-suffix contig
       012345    <-contig coord
    ...789012... <-genome coord
            ^ref_clv
    """
    clv = 1
    strand = '+'
    ref_fa = MagicMock()
    contig = MagicMock()
    contig.query_sequence = 'ATCGACAA'
    contig.cigartuples = ((S.BAM_CMATCH, 6), (S.BAM_CSOFT_CLIP, 2))

    args = contig, strand, clv, ref_fa
    assert extract_seq(*args) == 'ATCGAC'
    assert extract_seq(*args, window=1) ==    'C'
    assert extract_seq(*args, window=2) ==   'AC'
    assert extract_seq(*args, window=3) ==  'GAC'
    assert extract_seq(*args, window=4) == 'CGAC'


def test_extract_seq_for_minus_strand_clv_supported_by_suffix():
    """
    TTT         <-tail of suffix contig
      └ACATC    <-suffix contig
       01234    <-contig coord
    ...89012... <-genome coord
       ^ref_clv
    """
    clv = 1
    strand = '-'
    ref_fa = MagicMock()
    contig = MagicMock()
    contig.query_sequence = 'TTTACATCG'
    contig.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 6))

    args = contig, strand, clv, ref_fa
    assert extract_seq(*args) == 'ACATCG'
    assert extract_seq(*args, window=1) == 'A'
    assert extract_seq(*args, window=2) == 'AC'
    assert extract_seq(*args, window=3) == 'ACA'
    assert extract_seq(*args, window=4) == 'ACAT'


# TODO
# def test_extract_seq_where_ref_clv_is_not_at_the_end_of_the_contig_for_minus_strand_clv():
#     """
#     TT
#     |└AC TTTATT   <-bridge read
#     AACGG┘||||└CG <-bridge contig
#     0123456789012 <-contig coord
#     7890123456789 <-genome coordinate
#       ^ref_clv
#     """
#     c = MagicMock()
#     c.query_sequence = 'AACGGTTTATT'
#     c.cigartuples = ((S.BAM_CMATCH, 13))
#     assert extract_seq(c) == 'CGGTTTATT'


class TestExtractSeqForSoftClippedSeq(unittest.TestCase):
    def setUp(self):
        # these won't't affect the results
        self.clv = 0
        self.ref_fa = MagicMock()

    def test_extract_seq_with_starting_softclip(self):
        strand = '-'
        contig = MagicMock()
        contig.query_sequence = 'TTCCA'
        contig.cigartuples = ((S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 4))
        assert extract_seq(contig, strand, self.clv, self.ref_fa) == 'CCA'

    def test_extract_seq_with_ending_softclip(self):
        strand = '+'
        contig = MagicMock()
        contig.query_sequence = 'GGGAA'
        contig.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))
        assert extract_seq(contig, strand, self.clv, self.ref_fa) == 'GGG'

    def test_extract_seq_with_both_ends_clipped(self):
        strand = '-'
        contig = MagicMock()
        contig.query_sequence = 'TTTGGGAA'
        contig.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))
        assert extract_seq(contig, strand, self.clv, self.ref_fa) == 'GGG'


class TestExtractSeqForHardClippedSeq(unittest.TestCase):
    """The same to the above test, but replaced BAM_CSOFT_CLIP with BAM_CHARD_CLIP"""
    def setUp(self):
        # these won't't affect the results
        self.clv = 0
        self.ref_fa = MagicMock()

    def test_extract_seq_with_starting_softclip(self):
        strand = '-'
        contig = MagicMock()
        contig.query_sequence = 'TTCCA'
        contig.cigartuples = ((S.BAM_CHARD_CLIP, 2), (S.BAM_CMATCH, 3))
        assert extract_seq(contig, strand, self.clv, self.ref_fa) == 'CCA'

    def test_extract_seq_with_ending_softclip(self):
        strand = '+'
        contig = MagicMock()
        contig.query_sequence = 'GGGAA'
        contig.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CHARD_CLIP, 2))
        assert extract_seq(contig, strand, self.clv, self.ref_fa) == 'GGG'

    def test_extract_seq_with_both_ends_clipped(self):
        strand = '-'
        contig = MagicMock()
        contig.query_sequence = 'TTGGGCCAAA'
        contig.cigartuples = ((S.BAM_CHARD_CLIP, 2), (S.BAM_CMATCH, 5), (S.BAM_CHARD_CLIP, 3))
        assert extract_seq(contig, strand, self.clv, self.ref_fa) == 'GGGCC'
