import pytest

import kleat.misc.settings as S
from kleat.misc.apautils import calc_genome_offset


@pytest.mark.parametrize("ctg_clv, tail_side, skip_check_size, expected_gnm_offset", [
    # [1, 'left',  0, 1],
    # [1, 'right', 0, 1],
    # [2, 'left',  0, 3],
    # [2, 'right', 0, 3],

    [1, 'left',  1, 3],
    [1, 'right', 1, 1],
    [2, 'left',  1, 3],
    [2, 'right', 1, 1],
])
def test_clv_before_insertion(ctg_clv, tail_side, skip_check_size, expected_gnm_offset):
    """
   TT AA        <-bridge read tail
    └C┘         <-bridge read, for visual convenience two cases for different tail sides are merged with only one base shown
     |          # blank line to separate the bridge read the insertion
     | AGC      <-inserted sequence
     | 345      <-contig offset coord for inserted sequence
     |  ┬
    AT-G GT     <-contig
    01 2 678    <-contig offset coord
     ^ctg_clv
    0123 456    <-genome offset coord
     ^gnm_offset
    """
    ctg_cigartuples = (
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 1),
        (S.BAM_CINS, 3),
        (S.BAM_CMATCH, 2)
    )
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, tail_side, skip_check_size) == expected_gnm_offset



@pytest.mark.parametrize("ctg_clv, tail_side, skip_check_size, expected_gnm_offset", [
    [3, 'left',  0, 3],
    [3, 'right', 0, 4],

    [5, 'left',  0, 3],
    [5, 'right', 0, 4],

    [3, 'left',  1, 3],
    [3, 'right', 1, 4],

    [5, 'left',  1, 3],
    [5, 'right', 1, 4],

    [3, 'left',  2, 3],
    # this case is considered fine for now
    [3, 'right', 2, 4],

    [2, 'left',  2, 3],
    [2, 'right', 2, 1],

    [5, 'left',  2, 3],
    [5, 'right', 2, 4],
])
def test_clv_inside_insertion(ctg_clv, tail_side, skip_check_size, expected_gnm_offset):
    """
            TT AA      <-bridge read tail
             └A┘       <-bridge read, for visual convenience two cases for different tail sides are merged with only one base shown
              |        # blank line to separate the bridge read the insertion
              AGC      <-inserted sequence
              345      <-contig offset coord for inserted sequence
       ctg_clv^|
              |┬
           AT-G GT     <-contig
           01 2 678    <-contig offset coord
           0123 456    <-genome offset coord
              ^gnm_offset
    """
    ctg_cigartuples = (
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 1),
        (S.BAM_CINS, 3),
        (S.BAM_CMATCH, 2)
    )
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, tail_side, skip_check_size) == expected_gnm_offset


@pytest.mark.parametrize("ctg_clv, tail_side, skip_check_size, expected_gnm_offset", [
    [6, 'left',  0, 4],
    [6, 'right', 0, 4],

    [7, 'left',  0, 5],
    [7, 'right', 0, 5],

    [6, 'left',  2, 4],
    [6, 'right', 2, 4],

    [6, 'left',  3, 4],
    [6, 'right', 3, 4],

    [6, 'left',  4, 4],
    [6, 'right', 4, 4],

    [6, 'left',  5, 4],
    [6, 'right', 5, 4],

    [6, 'left',  6, 4],
    [6, 'right', 6, 4],

    # varying skip_check_size won't change the result the current
    # implementation only considers the last skip. Here the insertion is
    # between skip and match when tail_side is right
])
def test_clv_after_insertion(ctg_clv, tail_side, skip_check_size, expected_gnm_offset):
    """
            TT AA      <-bridge read tail
             └A┘       <-bridge read, for visual convenience two cases for different tail sides are merged with only one base shown
              |        # blank line to separate the bridge read the insertion
              AGC      <-inserted sequence
              345      <-contig offset coord for inserted sequence
               ┬
           AT-G GT     <-contig
           01 2 678    <-contig offset coord
                ^ctg_clv
           0123 456    <-genome offset coord
                ^gnm_offset
    """
    ctg_cigartuples = (
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 1),
        (S.BAM_CINS, 3),
        (S.BAM_CMATCH, 2)
    )
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, tail_side, skip_check_size) == expected_gnm_offset

