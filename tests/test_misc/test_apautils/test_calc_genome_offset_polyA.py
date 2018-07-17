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


def test_for_contig_with_two_skips_and_soft_clip_with_clv_right_before_the_skip():
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
