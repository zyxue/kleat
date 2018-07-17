from unittest.mock import MagicMock

from kleat.evidence import bridge
from kleat.misc.apautils import gen_clv_key_tuple
import kleat.misc.settings as S


def test_analyze_left_tail_bridge_read_aligned_to_a_forward_contig():
    """
    TTT
     |└ACG   <-left-tail read
     XXACGXX <-contig
     0123456 <-contig coord
       ^ctg_clv
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
    c.infer_query_length.return_value = 7
    c.cigartuples = ((S.BAM_CMATCH, 7),)

    ref_clv = 2                 # =ctg_clv because c.reference_end = 0
    tail_len = 3

    # just for the sake of fullfilling API requirement
    mock_dd_bridge = bridge.init_evidence_holder()
    clv_key = gen_clv_key_tuple('chr1', '-', ref_clv)
    mock_dd_bridge['hexamer_tuple'][clv_key] = ('NA', -1, -1)

    ref_fa = MagicMock()
    assert bridge.analyze_bridge(c, r, ref_fa, mock_dd_bridge, bridge_skip_check_size=0) == ('chr1', '-', ref_clv, tail_len, None)


def test_analyze_right_tail_bridge_read_aligned_to_a_forward_contig():
    """
         AA
      CCG┘|   <-right-tail read
     XXCCGXX  <-contig
     0123456  <-contig coord
        ^ctg_clv
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
    c.infer_query_length.return_value = 7
    c.cigartuples = ((S.BAM_CMATCH, 7),)

    ref_clv = 3                 # ctg_clv because c.reference_end = 0
    tail_len = 2

    # just for the sake of fullfilling API requirement
    mock_dd_bridge = bridge.init_evidence_holder()
    clv_key = gen_clv_key_tuple('chr1', '+', ref_clv)
    mock_dd_bridge['hexamer_tuple'][clv_key] = ('NA', -1, -1)

    ref_fa = MagicMock()

    assert bridge.analyze_bridge(c, r, ref_fa, mock_dd_bridge, bridge_skip_check_size=0) == ('chr1', '+', ref_clv, tail_len, None)
