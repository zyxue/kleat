from collections import defaultdict

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


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset_cutoff, gnm_offset", [
    [((BAM_CMATCH, 10),), 2, 2],
    [((BAM_CMATCH, 20),), 5, 5],
])
def test_calc_genome_offset_for_nonskipped_contig(
        ctg_cigartuples, ctg_offset_cutoff, gnm_offset):
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset_cutoff) == gnm_offset


@pytest.mark.parametrize("ctg_cigartuples, ctg_offset_cutoff, gnm_offset", [
    # ctg_offset_cutoff before skip happens
    [((BAM_CMATCH, 10), (BAM_CREF_SKIP, 5), (BAM_CMATCH, 10)), 2, 2],
    # ctg_offset_cutoff after skip happens
    [((BAM_CMATCH, 10), (BAM_CREF_SKIP, 5), (BAM_CMATCH, 10)), 12, 17],
])
def test_calc_genome_offset_for_skipped_contig(
        ctg_cigartuples, ctg_offset_cutoff, gnm_offset):
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset_cutoff) == gnm_offset

@pytest.mark.parametrize("ctg_cigartuples, ctg_offset_cutoff, gnm_offset", [
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CMATCH, 44)), 5, 5],
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CMATCH, 44)), 31, 31],
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CMATCH, 44)), 32, 34],
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CMATCH, 44)), 33, 35],

    # insertion, softclip or hardclip shouldn't have an effect
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CINS, 100), (BAM_CMATCH, 44)), 5, 5],
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CSOFT_CLIP, 100), (BAM_CMATCH, 44)), 31, 31],
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CHARD_CLIP, 200), (BAM_CMATCH, 44)), 32, 34],
    [((BAM_CMATCH, 31), (BAM_CDEL, 2), (BAM_CINS, 300), (BAM_CMATCH, 44)), 33, 35],
])
def test_calc_genome_offset_for_skipped_contig_with_deletion(
        ctg_cigartuples, ctg_offset_cutoff, gnm_offset):
    assert bridge.calc_genome_offset(ctg_cigartuples, ctg_offset_cutoff) == gnm_offset
