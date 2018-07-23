from unittest.mock import MagicMock

from kleat.visaln.visaln import calc_xlim
import kleat.misc.settings as S


# def test_calc_xlim_simple():
#     """
#     ATC---GATC///
#     01234567890123  <-genome coordinate
#               1
#     """

#     contig = MagicMock()
#     contig.reference_start = 0
#     contig.reference_end = 10
#     contig.query_length = 13
#     contig.cigartuples = (
#         (S.BAM_CMATCH, 3),
#         (S.BAM_CREF_SKIP, 3),
#         (S.BAM_CMATCH, 4),
#         (S.BAM_CHARD_CLIP, 3),
#     )

#     assert calc_xlim(contig, 0, 1, predicted_clv=4, padding=1) == (-1, 3)
#     assert calc_xlim(contig, 1, 2, predicted_clv=4, padding=1) == (3, 13)


def test_calc_xlim_with_two_skips():
    """
    ATC--GATC-----CGC///
    012345678901234567890  <-Genome coordinate
              1         2
    """

    contig = MagicMock()
    contig.reference_start = 0
    contig.reference_end = 17
    contig.query_length = 20

    contig.cigartuples = (
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 4),
        (S.BAM_CREF_SKIP, 5),
        (S.BAM_CMATCH, 3),
        (S.BAM_CHARD_CLIP, 3),
    )

    assert calc_xlim(contig, 0, 1, clvs=[4], padding=1) == (-1, 3)
    assert calc_xlim(contig, 1, 2, clvs=[4], padding=1) == (3, 9)
    assert calc_xlim(contig, 2, 3, clvs=[4], padding=1) == (13, 20)

    assert calc_xlim(contig, 1, 2, clvs=[4], padding=2) == (3, 10)
    assert calc_xlim(contig, 2, 3, clvs=[12], padding=2) == (10, 21)
