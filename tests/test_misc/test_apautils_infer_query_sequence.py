from unittest.mock import MagicMock

import kleat.misc.settings as S
from kleat.misc.apautils import infer_query_sequence


def test_infer_query_sequence_for_with_left_hard_clip():
    """
    CGATT    <-contig
    \\       <-hardclip mask
    012345   <-contig coordinate
    """

    ctg = MagicMock()
    ctg.query_sequence = 'ATT'
    ctg.get_tag.return_value = 'CG'
    ctg.cigartuples = [
        (S.BAM_CHARD_CLIP, 2),
        (S.BAM_CMATCH, 3),
    ]
    assert infer_query_sequence(ctg) == ('ATT')
    assert infer_query_sequence(ctg, always=True) == ('CGATT')


def test_infer_query_sequence_for_with_right_hard_clip():
    """
    CGATT    <-contig
       //    <-hardclip mask
    012345   <-contig coordinate
    """

    ctg = MagicMock()
    ctg.query_sequence = 'CGA'
    ctg.get_tag.return_value = 'TT'
    ctg.cigartuples = [
        (S.BAM_CMATCH, 3),
        (S.BAM_CHARD_CLIP, 2),
    ]
    assert infer_query_sequence(ctg) == ('CGA')
    assert infer_query_sequence(ctg, always=True) == ('CGATT')
