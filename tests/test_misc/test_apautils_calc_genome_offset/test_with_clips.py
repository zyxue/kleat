import pytest

from kleat.misc.calc_genome_offset import calc_genome_offset
import kleat.misc.settings as S


@pytest.mark.parametrize("ctg_clv, expected_gnm_offset", [
    [3, 0],                   # overlap with clv extracted from suffix evidence
    [4, 1],                   # bridge tail is a bit after the contig tail
])
def test_for_contig_with_softclip(ctg_clv, expected_gnm_offset):
    """
    TTT
    012       <-contig offset coord for tail
      └XXX    <-contig
       345    <-contig offset coord
       ^ctg_clv
       012    <-genome offset coord
       ^gnm_offset
    """
    ctg_cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_clv, expected_gnm_offset", [
    [2, 0],
    [3, 1],
    [4, 2],
])
def test_for_contig_with_hardclip(ctg_clv, expected_gnm_offset):
    """
    TT
     └ACG
    \\ACGT     <-contig
    0123456    <-contig offset coord
      ^ctg_clv
      0123     <-genome offset coord
      ^gnm_offset
    """
    ctg_cigartuples = ((S.BAM_CHARD_CLIP, 2), (S.BAM_CMATCH, 3))
    assert calc_genome_offset(ctg_cigartuples, ctg_clv, 'left') == expected_gnm_offset


def test_for_bridge_support_on_suffix_contig_edgecase():
    """
                 AA
           CGTACT┘|      <-bridge read
           012345678
           |||AAA  |
       ATGACGT┘ |  |     <-suffix contig
       0123456789012     <-contig offset coord
           |    ^ctg_clv
       0123456789012     <-genome offset coord
                ^gnm_offset
    """
    ctg_cigartuples = (
        (S.BAM_CMATCH, 7),
        (S.BAM_CSOFT_CLIP, 3)
    )

    expected_gnm_offset = 9
    assert calc_genome_offset(ctg_cigartuples, ctg_clv=9, tail_side='right') == expected_gnm_offset


# not sure if such case would ever happen yet, but test it anyway
def test_for_contig_with_clv_in_hardclip_spanning_clv():
    """
    TT
     └ACG
    CGACGTA      <-contig
    \\\   |      <-hardclip mask
    01234567     <-contig coord
   76543210      <-rev contig coord
      ^ctg_clv
      |01234     <-genome offset coord
   76543210      <-rev genome offset coord
      ^gnm_offset
    """
    ctg_cigartuples = (
        (S.BAM_CHARD_CLIP, 3),
        (S.BAM_CMATCH, 4),
    )
    expected_gnm_offset = -1
    assert calc_genome_offset(ctg_cigartuples, ctg_clv=2, tail_side='left') == expected_gnm_offset

    # just shift ctg_clv to the left a bit more
    assert calc_genome_offset(ctg_cigartuples, ctg_clv=1, tail_side='left') == -2
