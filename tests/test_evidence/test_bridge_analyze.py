from unittest.mock import MagicMock

from kleat.evidence import bridge
import kleat.misc.settings as S


def test_analyze_left_tail_bridge_read_aligned_to_a_forward_contig():
    """
     TTT
       └ACG   <-left-tail read
      XXACGX  <-contig
      123456  <-contig coord
     01234567 <-contig coord
    """
    r = MagicMock()
    r.reference_start = 2
    r.reference_end = 5
    r.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))
    r.query_sequence = 'TTTACG'

    c = MagicMock()
    c.reference_name = 'chr1'
    c.is_reverse = False
    c.reference_start = 0
    c.reference_end = 7
    c.cigartuples = ((S.BAM_CMATCH, 7),)

    assert bridge.analyze_bridge(c, r) == ('chr1', '-', 3, 3)


def test_analyze_right_tail_bridge_read_aligned_to_a_forward_contig():
    """
          AA
       CCG┘   <-right-tail read
      XCCGXX  <-contig
     0123456  <-contig coord
     01234567 <-contig coord
    """
    r = MagicMock()
    r.reference_start = 1
    r.reference_end = 4
    r.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))
    r.query_sequence = 'CCGAA'

    c = MagicMock()
    c.reference_name = 'chr1'
    c.is_reverse = False
    c.reference_start = 0
    c.reference_end = 7
    c.cigartuples = ((S.BAM_CMATCH, 7),)

    assert bridge.analyze_bridge(c, r) == ('chr1', '+', 4, 2)
