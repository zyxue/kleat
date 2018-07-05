from unittest.mock import MagicMock

from kleat.evidence import bridge
import kleat.misc.settings as S


def get_mock_read(reference_start, reference_end, cigartuples):
    r = MagicMock()
    r.reference_start = reference_start
    r.reference_end = reference_end
    r.cigartuples = cigartuples
    return r


def test_do_forwad_contig_left_tail_brdige_read():
    """
    e.g. TTTACG, reference_start at 2 points to the position of the first T
    (based on IGV)

     TTT
       └ACG  <-left-tail read
      XXACGX <-contig
     0123456 <-contig coord
    """
    mock_read = get_mock_read(2, 5, [(S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3)])
    ctg_offset = 3          # 2 + 1
    tail_len = 3
    assert bridge.do_fwd_ctg_lt_bdg(mock_read) == ('-', ctg_offset, tail_len)


def test_do_forwad_contig_left_tail_brdige_read_2():
    """
    e.g. TTAATTCCGG
    """
    mock_read = get_mock_read(10, 18, [(S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 8)])
    ctg_offset = 11         # 10 + 1
    tail_len = 2
    assert bridge.do_fwd_ctg_lt_bdg(mock_read) == ('-', ctg_offset, tail_len)


def test_do_forwad_contig_right_tail_brdige_read():
    """
    e.g. CCGGAA, reference_end at 4 points to the position of G (based on IGV)

          AA
       CCG┘  <-right-tail read
      XCCGXX <-contig
     0123456 <-contig coord
    """
    mock_read = get_mock_read(1, 4, [(S.BAM_CMATCH, 4), (S.BAM_CSOFT_CLIP, 2)])
    ctg_offset = 4
    tail_len = 2
    assert bridge.do_fwd_ctg_rt_bdg(mock_read) == ('+', ctg_offset, tail_len)


def test_do_reverse_contig_left_tail_brdige_read():
    """
    e.g. TTTACG
     TTT
       └ACG  <-left-tail read
      XXACGX <-contig
     0123456 <-contig coord
      6543210<-contig coord after reverse
    """
    mock_read = get_mock_read(2, 5, [(S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3)])
    # in genome coordinates, it's reversed, the the clv points to the position
    # of A, while position 0 point to the position after G.
    contig_len = 6
    ctg_offset = 4
    tail_len = 3
    assert bridge.do_rev_ctg_lt_bdg(mock_read, contig_len) == ('+', ctg_offset, tail_len)


def test_do_reverse_contig_right_tail_brdige_read():
    """
    e.g. CCGGAA, reference_end at 4 points to the position of G (based on IGV)

          AA
       CCG┘  <-right-tail read
      XCCGXX <-contig
     0123456 <-contig coord
      6543210<-contig coord after reverse
    """
    mock_read = get_mock_read(1, 4, [(S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 2)])
    # in genome coordinates, it's reversed, the the clv points to the position
    # of A, while position 0 point to the position after G.
    contig_len = 6
    ctg_offset = 3
    tail_len = 2
    assert bridge.do_rev_ctg_rt_bdg(mock_read, contig_len) == ('-', ctg_offset, tail_len)


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

