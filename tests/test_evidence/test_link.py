from collections import defaultdict
from unittest.mock import MagicMock

import pytest

from kleat.evidence import link
import kleat.misc.settings as S


def test_analyze_forward_link_for_polyA_read():
    """
      CCG-------AAA   <-right-tail read
    XXCCGX            <-contig
    0123456           <-genome coord
    """
    r = MagicMock()
    r.query_sequence = 'AAA'

    c = MagicMock()
    c.reference_start = 0
    c.reference_end = 6

    strand, ref_clv = link.analyze_forward_link(c, r) 
    assert strand == '+'
    assert ref_clv == 5


def test_analyze_forward_link_for_polyT_read():
    """
    TTT-------CCG    <-right-tail read
             XCCGXX  <-contig
            01234567 <-genome coord
    """
    r = MagicMock()
    r.query_sequence = 'TTT'

    c = MagicMock()
    c.reference_start = 1
    c.reference_end = 7

    strand, ref_clv = link.analyze_forward_link(c, r)
    assert strand == '-'
    assert ref_clv == 1


def test_analyze_reverse_link_for_polyA_read():
    """
    polyA read is polyT before reverse

    TTT-------CCG    <-right-tail read
             XCCGXX  <-contig
            01234567 <-genome coord
    """
    r = MagicMock()
    r.query_sequence = 'AAA'

    c = MagicMock()
    c.reference_start = 1
    c.reference_end = 7

    strand, ref_clv = link.analyze_reverse_link(c, r)
    assert strand == '-'
    assert ref_clv == 1


def test_analyze_reverse_link_for_polyT_reads():
    """
    polyT read is polyA before reverse

       CCG-------AAA   <-right-tail read
      XCCGXX           <-contig
    012345678          <-genome coord
    """
    r = MagicMock()
    r.query_sequence = 'TTT'

    c = MagicMock()
    c.reference_start = 2
    c.reference_end = 8

    strand, ref_clv = link.analyze_reverse_link(c, r)
    assert strand == '+'
    assert ref_clv == 7
