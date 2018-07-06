import pytest

import unittest
from unittest.mock import MagicMock

import kleat.misc.settings as S
from kleat.misc.search_hexamer import extract_seq


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
    def test_extract_seq_with_starting_softclip(self):
        c = MagicMock()
        c.query_sequence = 'TTCCA'
        c.cigartuples = ((S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 3))
        assert extract_seq(c) == 'CCA'

    def test_extract_seq_with_ending_softclip(self):
        c = MagicMock()
        c.query_sequence = 'GGGAA'
        c.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))
        assert extract_seq(c) == 'GGG'

    def test_extract_seq_with_both_ends_clipped(self):
        c = MagicMock()
        c.query_sequence = 'TTTGGGAA'
        c.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))
        assert extract_seq(c) == 'GGG'


class TestExtractSeqForHardClippedSeq(unittest.TestCase):
    """The same to the above test, but replaced BAM_CSOFT_CLIP with BAM_CHARD_CLIP"""
    def test_extract_seq_with_starting_softclip(self):
        c = MagicMock()
        c.query_sequence = 'TTCCA'
        c.cigartuples = ((S.BAM_CHARD_CLIP, 2), (S.BAM_CMATCH, 3))
        assert extract_seq(c) == 'CCA'

    def test_extract_seq_with_ending_softclip(self):
        c = MagicMock()
        c.query_sequence = 'GGGAA'
        c.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CHARD_CLIP, 2))
        assert extract_seq(c) == 'GGG'

    def test_extract_seq_with_both_ends_clipped(self):
        c = MagicMock()
        c.query_sequence = 'TTTGGGAA'
        c.cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 3), (S.BAM_CHARD_CLIP, 2))
        assert extract_seq(c) == 'GGG'
