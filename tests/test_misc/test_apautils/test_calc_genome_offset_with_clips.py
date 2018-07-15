import pytest

from kleat.misc.apautils import calc_genome_offset
import kleat.misc.settings as S


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
