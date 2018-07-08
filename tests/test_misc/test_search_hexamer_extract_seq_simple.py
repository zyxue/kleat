import unittest
from unittest.mock import MagicMock

import kleat.misc.settings as S
from kleat.misc.search_hexamer import extract_seq


def test_extract_seq_for_plus_strand_clv_supported_by_link():
    """
       ATCGAC    <-link contig
       012345    <-contig coord
    ...789012... <-genome coord
            ^ref_clv
    """
    clv = 99999                 # its value doesn't matter
    strand = '+'
    ref_fa = MagicMock()
    contig = MagicMock()
    contig.query_sequence = 'ATCGAC'
    contig.cigartuples = ((S.BAM_CMATCH, 6),)

    args = contig, strand, clv, ref_fa
    assert extract_seq(*args) == 'ATCGAC'
    assert extract_seq(*args, window=1) ==    'C'
    assert extract_seq(*args, window=2) ==   'AC'
    assert extract_seq(*args, window=3) ==  'GAC'
    assert extract_seq(*args, window=4) == 'CGAC'


def test_extract_seq_for_minus_strand_clv_supported_by_link():
    """
       ACATC    <-link contig
       01234    <-contig coord
    ...89012... <-genome coord
       ^ref_clv
    """
    clv = 99999
    strand = '-'
    ref_fa = MagicMock()
    contig = MagicMock()
    contig.query_sequence = 'ACATCG'
    contig.cigartuples = ((S.BAM_CMATCH, 6),)

    args = contig, strand, clv, ref_fa
    assert extract_seq(*args) == 'ACATCG'
    assert extract_seq(*args, window=1) == 'A'
    assert extract_seq(*args, window=2) == 'AC'
    assert extract_seq(*args, window=3) == 'ACA'
    assert extract_seq(*args, window=4) == 'ACAT'


def test_extract_seq_for_plus_strand_clv_supported_by_suffix():
    """
             AA  <-tail of suffix contig
       ATCGAC┘   <-suffix contig
       012345    <-contig coord
    ...789012... <-genome coord
            ^ref_clv
    """
    clv = 99999
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
    clv = 99999
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
