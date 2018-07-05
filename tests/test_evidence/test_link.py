from collections import defaultdict
from unittest.mock import MagicMock

import pytest

from kleat.evidence import link
import kleat.misc.settings as S


def test_analyze_forward_link_for_polyA_read():
    """
      CCG-------AAA   <-right-tail read
     XCCGXX           <-contig
    0123456           <-genome coord
    """
    r = MagicMock()
    r.query_sequence = 'AAA'

    c = MagicMock()
    c.reference_start = 0
    c.reference_end = 6

    assert link.analyze_forward_link(c, r) == ('+', 6)


def test_analyze_forward_link_for_polyT_read():
    """
    TTT-------CCG   <-right-tail read
             XCCGXX <-contig
            0123456 <-genome coord
    """
    r = MagicMock()
    r.query_sequence = 'TTT'

    c = MagicMock()
    c.reference_start = 0
    c.reference_end = 7

    assert link.analyze_forward_link(c, r) == ('-', 1)

