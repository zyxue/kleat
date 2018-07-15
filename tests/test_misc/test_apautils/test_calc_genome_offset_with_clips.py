import pytest

from kleat.misc.apautils import calc_genome_offset
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


# In such case, there won't be a gnm_offset, not sure if would ever happen yet,
# think about it later

# @pytest.mark.parametrize("ctg_clv, expected_gnm_offset", [
#     [2, 0],
#     [3, 1],
#     [4, 2],
# ])
# def test_for_contig_with_clv_in_hardclip(ctg_clv, expected_gnm_offset):
#     """
#     TT
#      └ACG
#     \\\CGT     <-contig
#     0123456    <-contig offset coord
#       ^ctg_clv
#       0123     <-genome offset coord
#       ^gnm_offset
#     """
#     ctg_cigartuples = ((S.BAM_CHARD_CLIP, 2), (S.BAM_CMATCH, 3))
#     assert calc_genome_offset(ctg_cigartuples, ctg_clv, 'left') == expected_gnm_offset
