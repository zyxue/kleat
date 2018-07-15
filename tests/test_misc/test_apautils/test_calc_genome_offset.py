import pytest

from kleat.misc.apautils import calc_genome_offset
import kleat.misc.settings as S

"""
When calculating offset, reversed contig should have already been flipped to
match the direction of the reference genome as illustrated in the ascii
drawings in the test cases here

The first example in pytest.mark.parametrize is often illustrated with ascii
drawing
"""


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((S.BAM_CMATCH, 4),), 1, 1],

    [((S.BAM_CMATCH, 10),), 2, 2],
    [((S.BAM_CMATCH, 20),), 5, 5],
])
def test_for_nonskipped_contig(
        ctg_cigartuples, ctg_offset, gnm_offset):
    """
     TT
      └AC     <-bridge read
      AACG    <-bridge contig
      01234   <-contig offset coord: different from "contig coord", it doesn't consider clipped regions
       ^ctg_offset
      01234   <-genome offset coord
       ^gnm_offset
    """
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((S.BAM_CMATCH, 6), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 4)), 2, 2],

    [((S.BAM_CMATCH, 10), (S.BAM_CREF_SKIP, 5), (S.BAM_CMATCH, 10)), 2, 2],
])
def test_for_skipped_with_clv_before_a_skip(
        ctg_cigartuples, ctg_offset, gnm_offset):
    """
     TT
      └AC            <-bridge read
      AACGTA--ATCG    <-bridge contig
      012345  67890   <-contig offset coord
        ^
      0123456789012   <-genome offset coord
        ^ctf/gnm_offset
    """
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((S.BAM_CMATCH, 2), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 6)), 4, 6],

    [((S.BAM_CMATCH, 10), (S.BAM_CREF_SKIP, 5), (S.BAM_CMATCH, 10)), 12, 17],
])
def test_for_skipped_contig_with_clv_after_a_skip(
        ctg_cigartuples, ctg_offset, gnm_offset):
    """
          TT
          |└AC  <-bridge read
      CG--ATCGAT    <-bridge contig
      01  2345678   <-contig offset coord
            ^ctg_offset
      01234567890   <-genome offset coord
            ^gnm_offset
    """
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((S.BAM_CMATCH, 31), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 44)), 5, 5],  # case1

    [((S.BAM_CMATCH, 31), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 44)), 30, 30],  # case2

    [((S.BAM_CMATCH, 31), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 44)), 31, 33],  # case3
    [((S.BAM_CMATCH, 31), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 44)), 32, 34],
    [((S.BAM_CMATCH, 31), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 44)), 33, 35],
])
def test_for_a_long_contig_with_deletion(ctg_cigartuples, ctg_offset, gnm_offset):
    """
       TT                       TT TT
       |└TC                      └C └TC                                            <-bridge read
    ATCGATCGATCGATCGATCGATCGATCGATC__ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG  <-bridge contig
    0123456789012345678901234567890  123456789012345678901234567890123456789012345 <-bridge offset coord
         |    1         2         3  |        4         5         6         7
         ^ctg_offset       (case2)^  ^ctg_offset(case3)
    012345678901234567890123456789012345678901234567890123456789012345678901234567 <-genome offset coord
         |    1         2         3  |      4         5         6         7
         ^gnm_offset(case1)          ^gnm_offset(case3)
    """
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [3, 0],                   # overlap with clv extracted from suffix evidence
    [4, 1],                   # bridge tail is a bit after the contig tail
])
def test_for_contig_with_softclip(ctg_offset, expected_gnm_offset):
    """
    TTT
    012       <-contig offset coord for tail
      └XXX    <-contig
       345    <-contig offset coord
       012    <-genome offset coord

    see parameters in the decorator for various ctg_offset
    """
    ctg_cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [2, 0],                   # overlap with clv extracted from suffix evidence
    [3, 1],                   # bridge tail is a bit after the contig tail
])
def test_for_contig_with_hardclip(ctg_offset, expected_gnm_offset):
    """
    The calculation with soft-clipped contig is the same except the CIGAR

     TT
     01       <-contig offset coord for tail
      └XXX    <-contig
       234    <-contig offset coord
       012    <-genome offset coord

    """
    ctg_cigartuples = ((S.BAM_CHARD_CLIP, 2), (S.BAM_CMATCH, 3))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset
