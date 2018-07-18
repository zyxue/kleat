from unittest.mock import MagicMock, patch

import kleat.misc.settings as S
from kleat.hexamer.hexamer import extract_seq


"""
cc: ctg_clv; ice: init_clv_end
rc: ref_clv; ire: init_ref_end
"""


@patch('kleat.hexamer.hexamer.apautils')
def test_hardclip_plust_strand(mock_apautils):
    """
           AAA
    CGCACCG┘ |       <-suffix contig with hardclip
    \\\|  |  |       <-hardclip mask
    01234567890      <-contig coord
       |cc^   ^ice
 ...XXXACCGTCG...    <-genome
    234567890123     <-genome coord
          | 1 |
        rc^   ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'CGCACCGAAA'
    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 3),
        (S.BAM_CMATCH, 4),
        (S.BAM_CSOFT_CLIP, 3),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=8, ref_fa=ref_fa, ctg_clv=6)
    assert extract_seq(**kw) == 'CGCACCG'
    assert extract_seq(window=1, **kw) == 'G'
    assert extract_seq(window=3, **kw) == 'CCG'


@patch('kleat.hexamer.hexamer.apautils')
def test_hardclip_minus_strand(mock_apautils):
    """
      TT
      |└CGCACCG       <-suffix contig with hardclip
      | |   ///       <-hardclip mask
      01234567890      <-contig coord
   icb^ ^cc     1
      XXCGCACCG...    <-genome
      3456789012      <-genome coord
      | |    1
   irb^ ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'TTCGCACCG'
    ctg.cigartuples = (
        (S.BAM_CSOFT_CLIP, 2),
        (S.BAM_CMATCH, 4),
        (S.BAM_CHARD_CLIP, 3),
    )
    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=5, ref_fa=ref_fa, ctg_clv=2)
    assert extract_seq(**kw) == 'CGCACCG'
    assert extract_seq(window=1, **kw) == 'C'
    assert extract_seq(window=3, **kw) == 'CGC'

