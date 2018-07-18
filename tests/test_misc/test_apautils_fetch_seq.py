import os
import io

import pysam

from kleat.misc.apautils import fetch_seq

REF_FA_FILE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), 'test_hexamer', 'mock_seq.fa')
REF_FA = pysam.FastaFile(REF_FA_FILE)


def test_fetch_seq_for_circular_DNA_with_positive_beginning_and_ending_less_than_contig_length():
    """
    AATTCCGG
    012345678
    """
    assert fetch_seq(REF_FA, seqname='chrM', beg=0, end=8) == 'AATTCCGG'
    assert fetch_seq(REF_FA, seqname='chrM', beg=2, end=6) == 'TTCC'


def test_fetch_seq_for_circular_DNA_when_beginning_is_negative():
    """
    AATTCCGG
    012345678
    """
    assert fetch_seq(REF_FA, seqname='chrM', beg=-1, end=1) == 'GA'
    assert fetch_seq(REF_FA, seqname='chrM', beg=-2, end=3) == 'GGAAT'


def test_fetch_seq_for_circular_DNA_when_ending_is_beyond_contig_length():
    """
    AATTCCGG
    012345678
    """
    assert fetch_seq(REF_FA, seqname='chrM', beg=6, end=11) == 'GGAAT'
    assert fetch_seq(REF_FA, seqname='chrM', beg=7, end=9) == 'GA'


def test_fetch_seq_for_circular_DNA_when_ending_beg_is_larger_than_end():
    """
    AATTCCGG
    01234567890
    """
    assert fetch_seq(REF_FA, seqname='chrM', beg=8, end=1) == 'A'
    assert fetch_seq(REF_FA, seqname='chrM', beg=8, end=2) == 'AA'
    assert fetch_seq(REF_FA, seqname='chrM', beg=8, end=7) == 'AATTCCG'

    assert fetch_seq(REF_FA, seqname='chrM', beg=9, end=1) == ''
    assert fetch_seq(REF_FA, seqname='chrM', beg=9, end=2) == 'A'
    assert fetch_seq(REF_FA, seqname='chrM', beg=9, end=3) == 'AT'
    assert fetch_seq(REF_FA, seqname='chrM', beg=9, end=3) == 'AT'


def test_fetch_seq_for_circular_DNA_when_ending_beg_is_larger_than_end_set_seq_len_is_ten_for_convenience():
    """
              ┬          ┬
    AATTCCGGAC AATTCCGGAC AATTCCGGAC
    0123456789 0123456789 0123456789
    """
    assert fetch_seq(REF_FA, seqname='MT', beg=8, end=1) == 'ACA'
    assert fetch_seq(REF_FA, seqname='MT', beg=8, end=2) == 'ACAA'

    assert fetch_seq(REF_FA, seqname='MT', beg=10, end=1) == 'A'
    assert fetch_seq(REF_FA, seqname='MT', beg=10, end=2) == 'AA'
    assert fetch_seq(REF_FA, seqname='MT', beg=10, end=3) == 'AAT'
    assert fetch_seq(REF_FA, seqname='MT', beg=10, end=4) == 'AATT'

    assert fetch_seq(REF_FA, seqname='MT', beg=11, end=1) == ''
    assert fetch_seq(REF_FA, seqname='MT', beg=11, end=2) == 'A'
    assert fetch_seq(REF_FA, seqname='MT', beg=11, end=3) == 'AT'

    assert fetch_seq(REF_FA, seqname='MT', beg=21, end=1) == ''
    assert fetch_seq(REF_FA, seqname='MT', beg=21, end=3) == 'AT'


def test_fetch_seq_for_linear_DNA_with_positive_beginning_and_ending_less_than_contig_length():
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=0, end=8) == 'AATTCCGG'
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=2, end=6) == 'TTCC'


def test_fetch_seq_for_linear_DNA_when_beginning_is_negative():
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=-1, end=1) == 'A'
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=-2, end=3) == 'AAT'


def test_fetch_seq_for_linear_DNA_when_ending_is_beyond_contig_length():
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=6, end=11) == 'GG'
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=7, end=9) == 'G'
