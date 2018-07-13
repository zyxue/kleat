import pytest

from kleat.misc.apautils import calc_genome_offset
import kleat.misc.settings as S

"""When calculating offset, reversed contig should have already been flipped to
match the strand of the genome as illustrated in the ascii drawings in the test
cases here"""


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((S.BAM_CMATCH, 4),), 1, 1],  # with ascii drawing below

    [((S.BAM_CMATCH, 10),), 2, 2],
    [((S.BAM_CMATCH, 20),), 5, 5],
])
def test_calc_genome_offset_for_nonskipped_contig(
        ctg_cigartuples, ctg_offset, gnm_offset):
    """
     TT
      └AC  <-bridge read
      AACG <-bridge contig
      0123 <-contig coord
      0123 <-genome offset
       ^ctf/gnm_offset
    """
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((S.BAM_CMATCH, 6), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 4)), 2, 2],  # with ascii drawing below
    [((S.BAM_CMATCH, 10), (S.BAM_CREF_SKIP, 5), (S.BAM_CMATCH, 10)), 2, 2],
])
def test_calc_genome_offset_for_skipped_contig_when_ctg_offset_is_before_skipping_happens(
        ctg_cigartuples, ctg_offset, gnm_offset):
    """
      TT
      |└AC  <-bridge read
      AACGTA--ATCG <-bridge contig
      012345  6789 <-contig coord
      012345678901 <-genome offset
        ^ctf/gnm_offset
    """
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((S.BAM_CMATCH, 2), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 6)), 4, 6],  # with ascii drawing below

    [((S.BAM_CMATCH, 10), (S.BAM_CREF_SKIP, 5), (S.BAM_CMATCH, 10)), 12, 17],
])
def test_calc_genome_offset_for_skipped_contig_when_ctg_offset_is_after_skipping_happens(
        ctg_cigartuples, ctg_offset, gnm_offset):
    """
          TT
          |└AC  <-bridge read
      CG--ATCGAT <-bridge contig
      01  234567 <-contig coord
      0123456789 <-genome offset
            ^ctf/gnm_offset
    """
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((S.BAM_CMATCH, 31), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 44)), 5, 5],
    [((S.BAM_CMATCH, 31), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 44)), 31, 31],
    [((S.BAM_CMATCH, 31), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 44)), 32, 34],
    [((S.BAM_CMATCH, 31), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 44)), 33, 35],
])
def test_calc_genome_offset_for_contig_with_deletion(ctg_cigartuples, ctg_offset, gnm_offset):
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    # before the insertion, see example in the docstring
    [2, 2],
    [3, 3],

    # inside the insertion, here only left-tail case is tested, see
    # ./test_bridge_clv_inside_insertion.py for more comprehensive test cases
    [4, 3],
    [5, 3],
    [6, 3],

    # after the insertion
    [7, 4],
    [8, 5],
])
def test_calc_genome_offset_for_contig_with_three_base_insertion(ctg_offset, expected_gnm_offset):
    """
       AGC  <-inserted sequence
       456  <-contig coord for inserted sequence
        ┬
    XXXX XX <-contig
    0123 78 <-contig coord
    0123 45 <-genome offset

    see parameters in the decorator for various ctg_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CINS, 3), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    # before the insertion, see example in the docstring
    [2, 2],
    [3, 3],

    # in the insertion
    [4, 3],

    # after the insertion
    [5, 4],
    [6, 5],
])
def test_calc_genome_offset_for_contig_with_one_base_insertion(ctg_offset, expected_gnm_offset):
    # thought one-base case might be more suitable for testing edgecases
    """
        G   <-inserted sequence
        4   <-contig coord for inserted sequence
        ┬
    XXXX XX <-contig
    0123 56 <-contig coord
    0123 45 <-genome offset

    see parameters in the decorator for various ctg_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CINS, 1), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [3, 0],                   # overlap with clv extracted from suffix evidence
    [4, 1],                   # bridge tail is a bit after the contig tail
])
def test_calc_genome_offset_for_contig_with_softclip(ctg_offset, expected_gnm_offset):
    """
    TTT
    012    <-contig coord for tail
      └XXX <-contig
       345 <-contig coord
       012 <-genome offset

    see parameters in the decorator for various ctg_offset
    """
    ctg_cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [2, 0],                   # overlap with clv extracted from suffix evidence
    [3, 1],                   # bridge tail is a bit after the contig tail
])
def test_calc_genome_offset_for_contig_with_hardclip(ctg_offset, expected_gnm_offset):
    """
    The calculation with soft-clipped contig is the same except the CIGAR

     TT
     01    <-contig coord for tail
      └XXX <-contig
       234 <-contig coord
       012 <-genome offset

    """
    ctg_cigartuples = ((S.BAM_CHARD_CLIP, 2), (S.BAM_CMATCH, 3))
    assert calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset
