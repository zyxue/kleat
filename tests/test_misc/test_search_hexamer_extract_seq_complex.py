import unittest
from unittest.mock import MagicMock

import kleat.misc.settings as S
from kleat.misc.search_hexamer import extract_seq


def test_extract_seq_with_skipped_region_for_plus_strand_clv():
    """
               AAA             <-tail of suffix contig
       ACGG--GC┘||             <-suffix contig with skip
       0123  456789            <-contig coord
              |  1
       ctg_clv^  ^init_ctg_idx <-contig coord
    ...ACGGTTGCGGT...          <-genome
       789012345678            <-genome coord
          1   |  |
       ref_clv^  ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    ctg.query_sequence = 'ACGGGCAAA'  # len 9
    ctg.cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 2), (S.BAM_CSOFT_CLIP, 3))

    ref_fa = MagicMock()
    ref_fa.fetch = MagicMock(return_value='TT')
    assert extract_seq(contig=ctg, strand='+', ref_clv=14, ref_fa=ref_fa, ctg_clv=5) == 'ACGGTTGC'
    ref_fa.fetch.assert_called_with('chr1', 11, 13)


def test_extract_seq_with_skipped_region_for_minus_strand_clv():
    """
    cc: ctg_clv; ici: init_clv_id
    rc: ref_clv; iri: init_ref_id

                TT          <-tail of suffix contig
                |└GT--C     <-suffix contig with skip
                01234  56   <-contig coord
        init_ctg_clv^ ^ctg_clv  <-contig coord
             ...ACGGTTGC... <-genome
                567890123   <-genome coord
                | | 1
    init_ref_idx^ ^ref_clv
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'TTGTC'
    ctg.cigartuples = ((S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 2), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 1))

    ref_fa = MagicMock()
    ref_fa.fetch = MagicMock(return_value='TG')
    assert extract_seq(contig=ctg, strand='-', ref_clv=8, ref_fa=ref_fa, ctg_clv=3) == 'GTTGC'
    ref_fa.fetch.assert_called_with('chr2', 10, 12)
