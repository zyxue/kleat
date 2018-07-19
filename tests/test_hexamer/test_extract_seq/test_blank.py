from unittest.mock import MagicMock, patch

import kleat.misc.settings as S
from kleat.hexamer.hexamer import extract_seq


"""
cc: ctg_clv; ice: init_clv_end
rc: ref_clv; ire: init_ref_end
"""

@patch('kleat.hexamer.hexamer.apautils')
def test_for_blank_contig_with_hardclip_plus_strand(mock_apautils):
    """
    derived_from_a_real_case

       ATGACGT      <-blank contig
       \\    |      <-hardclip mask
       01234567     <-contig offset coord
             ^ctg_clv
       56789012     <-genome offset coord
             ^ref_clv
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    mock_apautils.infer_query_sequence.return_value = 'ATGACGT'

    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 2),
        (S.BAM_CMATCH, 5),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=11, ref_fa=ref_fa, ctg_clv=6)
    assert extract_seq(**kw) == 'ATGACGT'


@patch('kleat.hexamer.hexamer.apautils')
def test_for_blank_contig_with_hardclip_minus_strand(mock_apautils):
    """
    derived_from_a_real_case

       ATGACGT      <-blank contig
       |    //      <-hardclip mask
       01234567     <-contig offset coord
       ^ctg_clv
       56789012     <-genome offset coord
       ^ref_clv
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    mock_apautils.infer_query_sequence.return_value = 'ATGACGT'

    ctg.cigartuples = (
        (S.BAM_CMATCH, 5),
        (S.BAM_CHARD_CLIP, 2),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=5, ref_fa=ref_fa, ctg_clv=0)
    assert extract_seq(**kw) == 'ATGACGT'
