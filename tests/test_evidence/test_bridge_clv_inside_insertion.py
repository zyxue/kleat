import pytest

from kleat.evidence import bridge
import kleat.misc.settings as S


"""in addition to bridge.py, test cases here validates the calculation of clv
when it happens to be inside an insertion in a tail-direction away way, i.e.

when it's left-tailed (always in the forward sense, if contig is reversed,
reverse it again), pick the left closest genome coordinate to represent the clv
location.

when it's right-tailed, pick the right closest genome coordinate to represent
the clv location.
"""

@pytest.mark.parametrize("ctg_offset_cutoff, tail_direction, expected_gnm_offset", [
    # before or after the insertion, tail_direction doesn't matter, see example
    # in the docstring

    # before the insertion
    [2, 'left', 2],
    [3, 'left', 3],
    [2, 'right', 2],
    [3, 'right', 3],

    # after the insertion
    [7, 'left', 4],
    [8, 'left', 5],
    [7, 'right', 4],
    [8, 'right', 5],

    # inside the insertion, tail_direction matters
    [4, 'left', 3],
    [5, 'left', 3],
    [6, 'left', 3],
    [4, 'right', 4],
    [5, 'right', 4],
    [6, 'right', 4],
])
def test_calc_genome_offset_for_contig_with_three_base_insertion_with_varying_tail_directions(ctg_offset_cutoff, tail_direction, expected_gnm_offset):
    """
       AGC  <-inserted sequence
       456  <-contig coord for inserted sequence
        ┬
     XXX XX <-contig
    0123 78 <-contig coord
    0123 45 <-genome offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CINS, 3), (S.BAM_CMATCH, 2))
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset_cutoff, tail_direction) == expected_gnm_offset



@pytest.mark.parametrize("ctg_offset_cutoff, tail_direction, expected_gnm_offset", [
    # before the insertion, see example in the docstring
    [2, 'left', 2],
    [3, 'left', 3],
    [2, 'right', 2],
    [3, 'right', 3],

    # inside the insertion
    [4, 'left', 3],
    [4, 'right', 4],

    # after the insertion
    [5, 'left', 4],
    [6, 'left', 5],
    [5, 'right', 4],
    [6, 'right', 5],
])
def test_calc_genome_offset_for_contig_with_one_base_insertion_with_varying_tail_directions(ctg_offset_cutoff, tail_direction, expected_gnm_offset):
    # thought one-base case would be more suitable for testing edgecases for if
    #
    # cur_ctg_ofs >= ctg_offset
    # or
    # cur_ctg_ofs > ctg_offset
    #
    # but turns out it doesn't matter, as the increase of the contig coordinate
    # already takes into consideration the >, = wouldn't happen.
    """
        G   <-inserted sequence
        4   <-contig coord for inserted sequence
        ┬
     XXX XX <-contig
    0123 56 <-contig coord
    0123 45 <-genome offset
    """
    ctg_cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CINS, 1), (S.BAM_CMATCH, 2))
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset_cutoff, tail_direction) == expected_gnm_offset
