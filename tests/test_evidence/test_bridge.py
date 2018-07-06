from collections import defaultdict
from unittest.mock import MagicMock

import pytest

from kleat.evidence import bridge
from kleat.misc.settings import (
    BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP,
    BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP)


def test_bridge_init_evidence_holder():
    assert bridge.init_evidence_holder() == {
        'num_reads': defaultdict(int),
        'max_tail_len': defaultdict(int),
    }


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((BAM_CMATCH, 4),), 1, 1],  # with ascii drawing below

    [((BAM_CMATCH, 10),), 2, 2],
    [((BAM_CMATCH, 20),), 5, 5],
])
def test_calc_genome_offset_for_nonskipped_forward_contig(
        ctg_cigartuples, ctg_offset, gnm_offset):
    """
     TT
      └AC  <-bridge read
      AACG <-bridge contig
      0123 <-contig coord
      0123 <-genome offset
    """
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    # ctg_offset before skip happens
    [((BAM_CMATCH, 10), (BAM_CREF_SKIP, 5), (BAM_CMATCH, 10)), 2, 2],
    # ctg_offset after skip happens
    [((BAM_CMATCH, 10), (BAM_CREF_SKIP, 5), (BAM_CMATCH, 10)), 12, 17],
])
def test_calc_genome_offset_for_skipped_contig(
        ctg_cigartuples, ctg_offset, gnm_offset):
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset, gnm_offset", [
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CMATCH, 44)), 5, 5],
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CMATCH, 44)), 31, 31],
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CMATCH, 44)), 32, 34],
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CMATCH, 44)), 33, 35],
])
def test_calc_genome_offset_for_contig_with_deletion(ctg_cigartuples, ctg_offset, gnm_offset):
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == gnm_offset


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
     XXX XX <-contig
    0123 78 <-contig coord
    0123 45 <-genome offset
    """
    ctg_cigartuples = ((BAM_CMATCH, 3), (BAM_CINS, 3), (BAM_CMATCH, 2))
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


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
     XXX XX <-contig
    0123 56 <-contig coord
    0123 45 <-genome offset
    """
    ctg_cigartuples = ((BAM_CMATCH, 3), (BAM_CINS, 1), (BAM_CMATCH, 2))
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset




@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [3, 1],                   # overlap with clv extracted from suffix evidence
    [4, 2],                   # bridge tail is a bit after the contig tail
])
def test_calc_genome_offset_for_contig_with_softclip(ctg_offset, expected_gnm_offset):
    """
     TT
    012
      └XXX <-contig
       345 <-contig coord
      0123 <-genome offset
    """
    ctg_cigartuples = ((BAM_CSOFT_CLIP, 2), (BAM_CMATCH, 3))
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset


@pytest.mark.parametrize("ctg_offset, expected_gnm_offset", [
    [3, 1],                   # overlap with clv extracted from suffix evidence
    [4, 2],                   # bridge tail is a bit after the contig tail
])
def test_calc_genome_offset_for_contig_with_hardclip(ctg_offset, expected_gnm_offset):
    """
    The calculation with soft-clipped contig is the same except the CIGAR

     TT
    012
      └XXX <-contig
       345 <-contig coord
      0123 <-genome offset

    """
    ctg_cigartuples = ((BAM_CHARD_CLIP, 2), (BAM_CMATCH, 3))
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset, 'left') == expected_gnm_offset
