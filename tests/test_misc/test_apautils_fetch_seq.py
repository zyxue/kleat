import os
import io

import pysam

from kleat.misc.apautils import fetch_seq

REF_FA_FILE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), 'test_hexamer', 'mock_seq.fa')
REF_FA = pysam.FastaFile(REF_FA_FILE)


def test_fetch_seq_for_circular_DNA_with_positive_beginning_and_ending_less_than_contig_length():
    assert fetch_seq(REF_FA, seqname='chrM', beg=0, end=8) == 'AATTCCGG'
    assert fetch_seq(REF_FA, seqname='chrM', beg=2, end=6) == 'TTCC'


def test_fetch_seq_for_circular_DNA_when_beginning_is_negative():
    assert fetch_seq(REF_FA, seqname='chrM', beg=-1, end=1) == 'GA'
    assert fetch_seq(REF_FA, seqname='chrM', beg=-2, end=3) == 'GGAAT'


def test_fetch_seq_for_circular_DNA_when_ending_is_beyond_contig_length():
    assert fetch_seq(REF_FA, seqname='chrM', beg=6, end=11) == 'GGAAT'
    assert fetch_seq(REF_FA, seqname='chrM', beg=7, end=9) == 'GA'


def test_fetch_seq_for_linear_DNA_with_positive_beginning_and_ending_less_than_contig_length():
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=0, end=8) == 'AATTCCGG'
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=2, end=6) == 'TTCC'


def test_fetch_seq_for_linear_DNA_when_beginning_is_negative():
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=-1, end=1) == 'A'
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=-2, end=3) == 'AAT'


def test_fetch_seq_for_linear_DNA_when_ending_is_beyond_contig_length():
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=6, end=11) == 'GG'
    assert fetch_seq(REF_FA, seqname='chr_mock1', beg=7, end=9) == 'G'
