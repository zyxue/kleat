import unittest
from unittest.mock import MagicMock, call

import kleat.misc.settings as S
from kleat.hexamer.search import extract_seq


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


# def test_extract_seq_where_for_plus_strand_clv_supported_by_suffix_with_specified_window_size():
#     """
#                       AA   <-tail of suffix contig
#        ACGGAATTCCGCGGCC┘   <-suffix contig
#        0123456789012345    <-contig coord
#               1      |
#     ...7890123456789012... <-genome coord
#           1         2|
#                      ^ref_clv
#     """
#     clv = 1
#     strand = '+'
#     ref_fa = MagicMock()
#     window = 3
#     contig = MagicMock()
#     contig.query_sequence = 'ACGGAATTCCCGGCCAA'
#     contig.cigartuples = ((S.BAM_CMATCH, 15), (S.BAM_CSOFT_CLIP, 2))

    # assert extract_seq(contig, strand, clv, ref_fa) == 'ACGGAATTCCCGGCC'


# def test_extract_seq_where_for_minus_strand_clv_supported_by_suffix():
#     """
#      TTT                <-tail of suffix contig
#        └ACATCGATCGGC    <-suffix contig
#         012345678901    <-contig coord
#         |         1
#      ...890123456789... <-genome coord
#         | 1
#         ^ref_clv
#     """
#     clv = 11
#     strand = '+'
#     contig = MagicMock()
#     contig.query_sequence = 'TTACATCGATCGC'
#     contig.cigartuples = ((S.BAM_CMATCH, 15), (S.BAM_CSOFT_CLIP, 2))
#     ref_fa = MagicMock()
#     assert extract_seq(contig, strand, clv, ref_fa) == 'ACATCGATCGC'


# def test_extract_seq_where_contig_has_skipped_cigars_for_plus_strand():
