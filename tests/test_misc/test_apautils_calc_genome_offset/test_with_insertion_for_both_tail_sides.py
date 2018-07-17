import pytest

import kleat.misc.settings as S
from kleat.misc.apautils import calc_genome_offset


"""
test cases here validates the calculation of clv when it happens to be
inside an insertion in a tail-side sensitive way, i.e.

when it's left-tailed (always in the forward sense, if contig is reversed,
reverse it again), pick the left closest genome coordinate to represent the clv
location.

when it's right-tailed, pick the right closest genome coordinate to represent
the clv location.
"""


@pytest.mark.parametrize("ctg_clv, tail_side, expected_gnm_offset", [
    [2, 'left',  2],
    [2, 'right', 2],

    [3, 'left',  3],
    [3, 'right', 3],

])
def test_for_contig_with_clv_before_insertion_so_tail_side_does_not_matter(ctg_clv, tail_side, expected_gnm_offset):
    """
    TT AA      <-bridge read tail
     └C┘       <-bridge read
      |        # blank line to separate the bridge read the insertion
      |AGC     <-inserted sequence
      |456     <-contig offset coord for inserted sequence
      | ┬
    ATCG GT    <-contig
    0123 789   <-contig offset coord
      ^ctg_clv
    0123 456   <-genome offset coord
      ^gnm_offset
    """
    ctg_cigartuples = (
        (S.BAM_CMATCH, 4),
        (S.BAM_CINS, 3),
        (S.BAM_CMATCH, 2)
    )
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, tail_side) == expected_gnm_offset


@pytest.mark.parametrize("ctg_clv, tail_side, expected_gnm_offset", [
    # inside the insertion, tail_side matters
    [4, 'left',  3],
    [4, 'right', 4],

    [5, 'left',  3],
    [5, 'right', 4],

    [6, 'left',  3],
    [6, 'right', 4],
])
def test_for_contig_with_clv_inside_insertion_so_tail_side_matters(ctg_clv, tail_side, expected_gnm_offset):
    """
            TT AA      <-bridge read tail
             └A┘       <-bridge read, for visual convenience two cases for different tail sides are merged with only one base shown
              |        # blank line to separate the bridge read the insertion
              AGC      <-inserted sequence
              456      <-contig offset coord for inserted sequence
       ctg_clv^|
              |┬
           ATCG GT    <-contig
           0123 78    <-contig offset coord
           0123 45    <-genome offset coord
              ^gnm_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 3), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, tail_side) == expected_gnm_offset


@pytest.mark.parametrize("ctg_clv, tail_side, expected_gnm_offset", [
    [7, 'left',  4],
    [7, 'right', 4],

    [8, 'left',  5],
    [8, 'right', 5],
])
def test_for_contig_with_clv_after_insertion_so_tail_side_does_not_matter(ctg_clv, tail_side, expected_gnm_offset):
    """
            TT AA      <-bridge read tail
             └A┘       <-bridge read, for visual convenience two cases for different tail sides are merged with only one base shown
              |        # blank line to separate the bridge read the insertion
              AGC      <-inserted sequence
              456      <-contig offset coord for inserted sequence
               ┬
           ATCG GT    <-contig
           0123 789   <-contig offset coord
                ^ctg_clv
           0123 456   <-genome offset coord
                ^gnm_offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 3), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, tail_side) == expected_gnm_offset


#################
# a 1-base case #
#################

@pytest.mark.parametrize("ctg_clv, tail_side, expected_gnm_offset", [
    [2, 'left',  2],
    [2, 'right', 2],

    [3, 'left',  3],
    [3, 'right', 3],
])
def test_for_contig_with_one_base_insertion_with_clv_before_insertion(ctg_clv, tail_side, expected_gnm_offset):
    """
    TT AA
     └C┘       <-bridge read, for visual convenience two cases for different tail sides are merged with only one base shown
      |        # blank line to separate the bridge read the insertion
      | G      <-inserted sequence
      | 4      <-contig offset coord for inserted sequence
      | ┬
    ATCG AC    <-contig
    0123 56    <-contig offset coord
      ^ctg_clv
    0123 45    <-genome offset coord
      ^gnm_offset
    see parameters in the decorator for various ctg_clv
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 1), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, tail_side) == expected_gnm_offset


@pytest.mark.parametrize("ctg_clv, tail_side, expected_gnm_offset", [
    [4, 'left', 3],
    [4, 'right', 4],
])
def test_for_contig_with_1_base_insertion_with_clv_inside_the_clv(ctg_clv, tail_side, expected_gnm_offset):
    """
      TT
       └GC     <-bridge read
        |      # blank line to separate the bridge read the insertion
        G      <-inserted sequence
        4      <-contig offset coord for inserted sequence
        ^ctg_clv
        ┬
    ATCG AC    <-contig
    0123 56    <-contig offset coord
    0123 45    <-genome offset coord
       ^gnm_offset
    see parameters in the decorator for various ctg_clv
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 1), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, tail_side) == expected_gnm_offset


@pytest.mark.parametrize("ctg_clv, tail_side, expected_gnm_offset", [
    [5, 'left',  4],
    [5, 'right', 4],

    [6, 'left',  5],
    [6, 'right', 5],
])
def test_for_contig_with_1_base_insertion_with_clv_after_insertion(ctg_clv, tail_side, expected_gnm_offset):
    """
       TT
        └AC    <-bridge read
         |     # blank line to separate the bridge read the insertion
        G|     <-inserted sequence
        4|     <-contig offset coord for inserted sequence
        ┬|
    ATCG AC    <-contig
    0123 56    <-contig offset coord
         ^ctg_clv
    0123 45    <-genome offset coord
         ^gnm_offset
    see parameters in the decorator for various ctg_clv
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CINS, 1), (S.BAM_CMATCH, 2))
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, tail_side) == expected_gnm_offset
