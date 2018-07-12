from unittest.mock import MagicMock, call

import kleat.misc.settings as S
from kleat.hexamer.search import extract_seq


"""
cc: ctg_clv; ici: init_clv_idx
rc: ref_clv; iri: init_ref_idx
"""


def test_extract_seq_with_skipped_region():
    """
                TTT             <-tail of suffix contig
                ||└GT--C        <-suffix contig with skip
                01234  56      <-contig coord
    init_ctg_clv^  ^ctg_clv     <-contig coord
                     | |
             ...XXXGTTGC...    <-genome
                5678901234      <-genome coord
                   | 1
                   ^ref_clv/iri
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'TTTGTC'
    ctg.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 2), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 1))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch = MagicMock(return_value='TG')
    kw = dict(contig=ctg, strand='-', ref_clv=8, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'GTTGC'
    ref_fa.fetch.assert_called_with('chr2', 10, 12)
    # **kw needs go after window for py34 syntax
    assert extract_seq(window=3, **kw) == 'GTT'


def test_extract_seq_with_two_skipped_regions_and_a_mismatch():
    """
                TTT              <-tail of suffix contig
                ||└GT--CAG-AC    <-suffix contig with skip
                01234  567 890   <-contig coord
    init_ctg_clv^  ^cc  x        <-contig coord
             ...XXXGTTGCGGCAC... <-genome
                56789012345678   <-genome coord
                   | 1
                   ^ref_clv/iri
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'TTTGTCAGAC'
    ctg.cigartuples = (
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 2),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch.side_effect = ['TG', 'C']
    kw = dict(contig=ctg, strand='-', ref_clv=8, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'GTTGCAGCAC'
    assert ref_fa.fetch.call_count == 2
    ref_fa.fetch.assert_has_calls([call('chr2', 10, 12), call('chr2', 15, 16)])


def test_extract_seq_with_2_base_insertion():
    """
              GA
         TTT  ┬       <-tail of suffix contig
         ||└AC TCG    <-suffix contig
         01234 567    <-contig coord
      ici^  ^ctg_clv
         XXXXX XXX    <-genome
         456789012... <-genome coord
            |  1
            ^ref_clv/iri
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    ctg.query_sequence = 'TTTACGATCG'
    ctg.cigartuples = (
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 2), 
        (S.BAM_CINS, 2),
        (S.BAM_CMATCH, 3),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=7, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'ACGATCG'
    assert extract_seq(window=3, **kw) == 'ACG'


def test_extract_seq_with_skipped_region_and_insertion_and_mismatches():
    """
               GA
         TTT   ┬       <-tail of suffix contig
         ||└A-C TAG    <-suffix contig
         0123 4 7890   <-contig coord
      ici^  ^ctc x
         XXXXGX XTX... <-genome
         456789 0123   <-genome coord
            |   1
            ^ref_clv/iri
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    ctg.query_sequence = 'TTTACGATAG'
    ctg.cigartuples = (
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 1),
        (S.BAM_CINS, 2),
        (S.BAM_CMATCH, 3),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch.return_value = 'G'
    kw = dict(contig=ctg, strand='-', ref_clv=7, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'AGCGATAG'
    ref_fa.fetch.assert_called_with('chr1', 8, 9)
    assert extract_seq(window=1, **kw) == 'A'
    assert extract_seq(window=5, **kw) == 'AGCGA'


def test_extract_seq_with_skipped_region_and_indels_and_mismatches():
    """
               GA
         TTT   ┬           <-tail of suffix contig
         ||└A-C TAG__GT    <-suffix contig
         0123 4 78901234   <-contig coord
      ici^  ^ctc x 1
         XXXXGX XTXTAXX... <-genome
         456789 01234567   <-genome coord
            |   1
            ^ref_clv/iri
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    ctg.query_sequence = 'TTTACGATAGGTA'
    ctg.cigartuples = (
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 1),
        (S.BAM_CINS, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CDEL, 2),
        (S.BAM_CMATCH, 2),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch.return_value = 'G'
    kw = dict(contig=ctg, strand='-', ref_clv=7, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'AGCGATAGGT'
    ref_fa.fetch.assert_called_with('chr1', 8, 9)
    assert extract_seq(window=9, **kw) == 'AGCGATAGG'
