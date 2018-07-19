from unittest.mock import MagicMock, patch

from kleat.evidence import bridge
from kleat.misc.apautils import gen_clv_key_tuple_with_ctg_clv
import kleat.misc.settings as S


# @patch('kleat.hexamer.hexamer.apautils')
def test_analyze_left_tail_bridge_read_aligned_to_a_forward_contig():
    """
    TTT
     |└ACG      <-left-tail read
     CGACGTA    <-contig
     01234567   <-contig coord
       ^ctg_clv
     34567890   <-genome coord
       ^ref_clv
    """
    read = MagicMock()
    read.reference_start = 2
    read.reference_end = 5
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))
    read.query_sequence = 'TTTACG'

    contig = MagicMock()
    contig.reference_name = 'chr1'
    contig.is_reverse = False
    contig.reference_start = 3
    contig.reference_end = 10
    contig.query_sequence = 'CGACGTA'
    contig.infer_query_length.return_value = 7
    contig.cigartuples = ((S.BAM_CMATCH, 7),)

    ref_clv = 5                 # =ctg_clv because contig.reference_end = 0
    ctg_clv = 2
    tail_len = 3

    # just for the sake of fullfilling API requirement
    mock_dd_bridge = bridge.init_evidence_holder()
    clv_key = gen_clv_key_tuple_with_ctg_clv('chr1', '-', ref_clv, ctg_clv)
    mock_dd_bridge['hexamer_tuple'][clv_key] = ('NA', -1, -1)

    ref_fa = MagicMock()
    expected = ('chr1', '-', ref_clv, ctg_clv, tail_len, None)
    assert bridge.analyze_bridge(contig, read, ref_fa, mock_dd_bridge, bridge_skip_check_size=0) == expected


def test_analyze_right_tail_bridge_read_aligned_to_a_forward_contig():
    """
         AA
      CCG┘|      <-right-tail read
     01234567    <-contig coord
     |  ^ctg_clv
     34567890    <-genome coord
        ^ref_clv
    """
    read = MagicMock()
    read.reference_start = 1
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))
    read.query_sequence = 'CCGAA'

    contig = MagicMock()
    contig.reference_name = 'chr1'
    contig.is_reverse = False
    contig.reference_start = 3
    contig.reference_end = 10
    contig.infer_query_length.return_value = 7
    contig.cigartuples = ((S.BAM_CMATCH, 7),)

    ref_clv = 6
    ctg_clv = 3
    tail_len = 2

    # just for the sake of fullfilling API requirement
    mock_dd_bridge = bridge.init_evidence_holder()
    clv_key = gen_clv_key_tuple_with_ctg_clv('chr1', '+', ref_clv, ctg_clv)
    mock_dd_bridge['hexamer_tuple'][clv_key] = ('NA', -1, -1)

    ref_fa = MagicMock()
    expected = ('chr1', '+', ref_clv, ctg_clv, tail_len, None)
    assert bridge.analyze_bridge(contig, read, ref_fa, mock_dd_bridge, bridge_skip_check_size=0) == expected
