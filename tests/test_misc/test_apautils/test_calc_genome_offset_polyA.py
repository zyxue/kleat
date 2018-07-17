import pytest

from kleat.misc.apautils import calc_genome_offset
import kleat.misc.settings as S


def test_for_skipped_contig_with_clv_right_before_the_skip():
    """
          AA
        AT┘         <-bridge read
      CGAT--CGAT    <-bridge contig
      0123  45678   <-contig offset coord
         ^ctg_clv
      01234567890   <-genome offset coord
         ^gnm_offset
    """
    cigartuples = [
        (S.BAM_CMATCH, 4),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 4)
    ]
    gnm_offset = 3
    assert calc_genome_offset(cigartuples, ctg_clv=3, tail_side='right') == gnm_offset


def test_for_contig_with_two_skips_with_clv_right_before_the_skip():
    """
             AA
           TC┘         <-bridge read
      CG--ATC--GAT    <-bridge contig
      01  234  5678   <-contig offset coord
            ^ctg_clv
      0123456789012   <-genome offset coord
            ^gnm_offset
    """
    cigartuples = [
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
    ]
    gnm_offset = 6
    assert calc_genome_offset(cigartuples, ctg_clv=4, tail_side='right') == gnm_offset


def test_for_contig_with_two_skips_and_soft_clip_with_clv_before_the_skip():
    """
            AA
          AT┘    AAAA    <-bridge read
      CG--ATC--GAT┘   <-bridge contig
      01  234  5678901   <-contig offset coord
           ^ctg_clv
      0123456789012   <-genome offset coord
           ^gnm_offset
    """
    cigartuples = [
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CSOFT_CLIP, 4),
    ]
    gnm_offset = 5
    assert calc_genome_offset(cigartuples, ctg_clv=3, tail_side='right') == gnm_offset


def test_for_contig_with_two_skips_and_soft_clip_edge_case_1():
    """
             AA
           TC┘    AAAA    <-bridge read
      CG--ATC--GAT┘   <-bridge contig
      01  234  56789012   <-contig offset coord
            ^ctg_clv
      0123456789012   <-genome offset coord
            ^gnm_offset
    """
    cigartuples = [
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CSOFT_CLIP, 4),
    ]
    gnm_offset = 6
    assert calc_genome_offset(cigartuples, ctg_clv=4, tail_side='right') == gnm_offset


# TODO: think if such is a realistic treatment of the clv, maybe it's better to
# assign clv to C(4), this would be caused by c2g alignment, the tool may
# wander around where the skip boundry should be. This would fix the real case below it
def test_for_contig_with_two_skips_and_soft_clip_edge_case_2():
    """
                AA
            C--G┘    AAAA    <-bridge read
      CG--ATC--GAT┘   <-bridge contig
      01  234  56789012   <-contig offset coord
               ^ctg_clv
      0123456789012   <-genome offset coord
               ^gnm_offset
    """
    cigartuples = [
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CSOFT_CLIP, 4),
    ]
    gnm_offset = 9
    assert calc_genome_offset(cigartuples, ctg_clv=5, tail_side='right') == gnm_offset


# def test_for_contig_real_case_E1_L_4362_chr16_plus_strand():
#     cigartuples = [
#         (S.BAM_CMATCH, 501),
#         (S.BAM_CREF_SKIP, 743),
#         (S.BAM_CMATCH, 169),
#         (S.BAM_CREF_SKIP, 661),
#         (S.BAM_CMATCH, 145),
#         (S.BAM_CREF_SKIP, 1005),
#         (S.BAM_CMATCH, 94),
#         (S.BAM_CREF_SKIP, 1712),
#         (S.BAM_CMATCH, 220),
#         (S.BAM_CREF_SKIP, 3379),
#         (S.BAM_CMATCH, 156),
#         (S.BAM_CREF_SKIP, 1043),
#         (S.BAM_CMATCH, 78),
#         (S.BAM_CREF_SKIP, 492),
#         (S.BAM_CMATCH, 155),
#         (S.BAM_CREF_SKIP, 101),
#         (S.BAM_CMATCH, 183),
#         (S.BAM_CREF_SKIP, 183),
#         (S.BAM_CMATCH, 508),
#         (S.BAM_CSOFT_CLIP, 13)
#     ]
#     gnm_offset = 5250
#     assert calc_genome_offset(cigartuples, ctg_clv=1129, tail_side='right') == gnm_offset
