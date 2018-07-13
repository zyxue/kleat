from collections import defaultdict
from unittest.mock import MagicMock

from kleat.evidence import bridge
from kleat.misc.apautils import gen_clv_key_tuple
import kleat.misc.settings as S


def test_bridge_init_evidence_holder():
    assert bridge.init_evidence_holder() == {
        'num_reads': defaultdict(int),
        'max_tail_len': defaultdict(int),
        'hexamer_tuple': defaultdict(lambda: None),
    }


def get_mock_read(ref_beg, ref_end, cigartuples):
    r = MagicMock()
    r.reference_start = ref_beg
    r.reference_end = ref_end
    r.cigartuples = cigartuples
    return r


def test_do_forwad_contig_left_tail_bridge_read():
    """e.g. TTTACG, reference_start at pos 2 (0-based) with the first three Ts
    soft-clipped

    compared to igv visualization, which is 1-based, the contig coord are
    shifted to the right by one-base, a quick comparison is available at
    http://zyxue.github.io/2018/06/21/coordinates-in-bioinformatics.html

    TTT
     |└ACG   <-left-tail read
     XXACGXX  <-contig
     0123456 <-contig coord
       ^ctg_offset
    """
    mock_read = get_mock_read(
        ref_beg=2, ref_end=5, cigartuples=[(S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3)])
    # basically the clv wst. to the contig coordinate when in forward contig
    ctg_offset = 2
    tail_len = 3
    assert bridge.do_fwd_ctg_lt_bdg(mock_read, contig=MagicMock()) == ('-', ctg_offset, tail_len)


def test_do_forwad_contig_left_tail_bridge_read_2():
    """
    TT
    |└AATTCCGG   <-left-tail read
    XXAATTCCGGXX <-contig
    890123456789 <-contig coord
      1
      ^ctg_offset
    """
    mock_read = get_mock_read(
        ref_beg=10, ref_end=18, cigartuples=[(S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 8)])
    ctg_offset = 10
    tail_len = 2
    assert bridge.do_fwd_ctg_lt_bdg(mock_read, contig=MagicMock()) == ('-', ctg_offset, tail_len)


def test_do_forwad_contig_right_tail_bridge_read():
    """
        AA
     CCG┘| <-right-tail read
    XXCCGXX <-contig
    0123456 <-contig coord
       ^ctg_offset
    """
    mock_read = get_mock_read(
        ref_beg=1, ref_end=4, cigartuples=[(S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2)])
    ctg_offset = 3
    tail_len = 2
    assert bridge.do_fwd_ctg_rt_bdg(mock_read, contig=MagicMock()) == ('+', ctg_offset, tail_len)


def test_do_reverse_contig_left_tail_bridge_read():
    """
           TTT
            |└ACG   <-left-tail read
            XXXACGX <-contig
            0123456 <-contig coord
            6543210 <-reversed contig coord, i.e. gnm offset from right to left
    ctg_offset^
    """
    mock_read = get_mock_read(
        ref_beg=2, ref_end=5, cigartuples=[(S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3)])
    # in genome coordinates, it's reversed, the the clv points to the position
    # of A, while position 0 point to the position after G.
    contig = MagicMock()
    contig.infer_query_length.return_value = 7
    ctg_offset = 4
    tail_len = 3
    assert bridge.do_rev_ctg_lt_bdg(mock_read, contig=contig) == ('+', ctg_offset, tail_len)


def test_do_reverse_contig_right_tail_bridge_read():
    """
               AA
            CCG┘|  <-right-tail read
           XXCCGXX <-contig
           0123456 <-contig coord
           6543210 <-reversed contig coord
    ctg_offset^
    """
    mock_read = get_mock_read(
        ref_beg=1, ref_end=4, cigartuples=[(S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 2)])
    # in genome coordinates, it's reversed, the the clv points to the position
    # of A, while position 0 point to the position after G.
    contig = MagicMock()
    contig.infer_query_length.return_value = 7
    ctg_offset = 3
    tail_len = 2
    assert bridge.do_rev_ctg_rt_bdg(mock_read, contig=contig) == ('-', ctg_offset, tail_len)


def test_analyze_left_tail_bridge_read_aligned_to_a_forward_contig():
    """
    TTT
     |└ACG   <-left-tail read
     XXACGXX <-contig
     0123456 <-contig coord
       ^ctg_offset
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

    ref_clv = 2                 # =ctg_offset because c.reference_end = 0
    tail_len = 3

    # just for the sake of fullfilling API requirement
    mock_dd_bridge = bridge.init_evidence_holder()
    clv_key = gen_clv_key_tuple('chr1', '-', ref_clv)
    mock_dd_bridge['hexamer_tuple'][clv_key] = ('NA', -1, -1)

    ref_fa = MagicMock()
    assert bridge.analyze_bridge(c, r, ref_fa, mock_dd_bridge) == ('chr1', '-', ref_clv, tail_len, None)


def test_analyze_right_tail_bridge_read_aligned_to_a_forward_contig():
    """
         AA
      CCG┘|   <-right-tail read
     XXCCGXX  <-contig
     0123456  <-contig coord
        ^ctg_offset
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

    ref_clv = 3                 # =ctg_offset because c.reference_end = 0
    tail_len = 2

    # just for the sake of fullfilling API requirement
    mock_dd_bridge = bridge.init_evidence_holder()
    clv_key = gen_clv_key_tuple('chr1', '+', ref_clv)
    mock_dd_bridge['hexamer_tuple'][clv_key] = ('NA', -1, -1)

    ref_fa = MagicMock()

    assert bridge.analyze_bridge(c, r, ref_fa, mock_dd_bridge) == ('chr1', '+', ref_clv, tail_len, None)
