import unittest
from unittest.mock import MagicMock
from unittest.mock import patch

import kleat.misc.settings as S
from kleat.hexamer.search import extract_seq


def test_extract_seq_for_plus_strand_clv_supported_by_link():
    """
       ATCGAC    <-link contig
       012345    <-contig coord
            ^ctg_clv
    ...789012... <-genome coord
          1 ^ref_clv
    """
    strand = '+'
    ref_clv = 12
    ctg_clv = 5
    ref_fa = MagicMock()
    contig = MagicMock()
    contig.query_sequence = 'ATCGAC'
    contig.cigartuples = ((S.BAM_CMATCH, 6),)

    args = contig, strand, ref_clv, ref_fa, ctg_clv
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
    strand = '-'
    ref_clv = 8
    ctg_clv = 0
    ref_fa = MagicMock()
    contig = MagicMock()
    contig.query_sequence = 'ACATCG'
    contig.cigartuples = ((S.BAM_CMATCH, 6),)

    args = contig, strand, ref_clv, ref_fa, ctg_clv
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
    strand = '+'
    ctg_clv = 5
    ref_clv = 12
    ref_fa = MagicMock()
    contig = MagicMock()
    contig.query_sequence = 'ATCGACAA'
    contig.cigartuples = ((S.BAM_CMATCH, 6), (S.BAM_CSOFT_CLIP, 2))

    args = contig, strand, ref_clv, ref_fa, ctg_clv
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
    strand = '-'
    ref_clv = 8
    ctg_clv = 0
    ref_fa = MagicMock()
    contig = MagicMock()
    contig.query_sequence = 'TTTACATCG'
    contig.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 6))

    args = contig, strand, ref_clv, ref_fa, ctg_clv
    assert extract_seq(*args) == 'ACATCG'
    assert extract_seq(*args, window=1) == 'A'
    assert extract_seq(*args, window=2) == 'AC'
    assert extract_seq(*args, window=3) == 'ACA'
    assert extract_seq(*args, window=4) == 'ACAT'


class TestExtractSeqForSoftClippedSeq(unittest.TestCase):
    def test_extract_seq_with_starting_softclip(self):
        """
        TT
         └CAA
        012345 <-contig coord
        678901 <-genome coord
        """
        strand = '-'
        contig = MagicMock()
        contig.query_sequence = 'TTCCA'
        contig.cigartuples = ((S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 4))
        assert extract_seq(contig, strand, ref_clv=8, ref_fa=MagicMock(), ctg_clv=2) == 'CCA'

    def test_extract_seq_with_ending_softclip(self):
        """
           AA
        GGG┘|
        012345   <-contig coord
        234567 <-genome coord
        """
        contig = MagicMock()
        contig.query_sequence = 'GGGAA'
        contig.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))
        assert extract_seq(contig, strand='+', ref_clv=7, ref_fa=MagicMock(), ctg_clv=2) == 'GGG'

    def test_extract_seq_with_both_ends_clipped(self):
        """
        TTT   AA
        ||└CCC┘|
        01234567
        90123456
         1 |
           ^clv
        """
        contig = MagicMock()
        contig.query_sequence = 'TTTCCCAA'
        contig.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))
        assert extract_seq(contig, strand='-', ref_clv=12, ref_fa=MagicMock(), ctg_clv=3) == 'CCC'


class TestExtractSeqForHardClippedSeq(unittest.TestCase):
    """The same to the above test, but replaced BAM_CSOFT_CLIP with BAM_CHARD_CLIP"""
    @patch('kleat.hexamer.search.apautils')
    def test_extract_seq_with_starting_softclip(self, mock_apautils):
        """
        TT
         └CAA
        012345 <-contig coord
        678901 <-genome coord
        """
        strand = '-'
        contig = MagicMock()
        mock_apautils.infer_query_sequence.return_value = 'TTCCA'
        contig.cigartuples = ((S.BAM_CHARD_CLIP, 2), (S.BAM_CMATCH, 4))
        assert extract_seq(contig, strand, ref_clv=8, ref_fa=MagicMock(), ctg_clv=2) == 'CCA'

    @patch('kleat.hexamer.search.apautils')
    def test_extract_seq_with_ending_softclip(self, mock_apautils):
        """
           AA
        GGG┘|
        012345   <-contig coord
        234567 <-genome coord
        """
        contig = MagicMock()
        mock_apautils.infer_query_sequence.return_value = 'GGGAA'
        contig.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CHARD_CLIP, 2))
        assert extract_seq(contig, strand='+', ref_clv=7, ref_fa=MagicMock(), ctg_clv=2) == 'GGG'

    @patch('kleat.hexamer.search.apautils')
    def test_extract_seq_with_both_ends_clipped(self, mock_apautils):
        """
        TTT   AA
        ||└CCC┘|
        01234567
        90123456
         1 |
           ^clv
        """
        contig = MagicMock()
        mock_apautils.infer_query_sequence.return_value = 'TTTCCCAA'
        contig.cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 3), (S.BAM_CHARD_CLIP, 2))
        assert extract_seq(contig, strand='-', ref_clv=12, ref_fa=MagicMock(), ctg_clv=3) == 'CCC'
