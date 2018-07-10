from unittest.mock import MagicMock, call

import kleat.misc.settings as S
from kleat.hexamer.search import extract_seq


"""
cc: ctg_clv; ice: init_clv_e
rc: ref_clv; ire: init_ref_end
"""


def test_extract_seq():
    """
             AA
           GT┘      <-bread read
       GACGGTTGC    <-bridge contig
       0123456789   <-contig coord
     ctg_clv^   ^ice  <-contig coord
    ...GACGGTTGC... <-genome
       5678901234    <-genome coord
            1   |
     ref_clv^   ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GACGGTTGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 9),)

    ref_fa = MagicMock()
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=5)
    assert extract_seq(**kw) == 'GACGGT'


def test_extract_seq_with_skip_after_ctg_clv():
    """
             AA
           GT┘      <-bread read
       GACGGT-GC    <-bridge contig
       012345 678   <-contig coord
     ctg_clv^   ^ice  <-contig coord
    ...GACGGTTGC... <-genome
       5678901234    <-genome coord
            1   |
     ref_clv^   ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GACGGTGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 6), (S.BAM_CREF_SKIP, 1), (S.BAM_CMATCH, 2))

    ref_fa = MagicMock()
    ref_fa.fetch.return_value = 'T'
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=5)
    assert extract_seq(**kw) == 'GACGGT'
    assert ref_fa.call_count == 0


def test_extract_seq_with_skip_before_ctg_clv():
    """
             AA
           GT┘      <-bread read
       G--AGTTGC    <-bridge contig
       0  1234567   <-contig coord
     ctg_clv^   ^ice  <-contig coord
    ...GACAGTTGC... <-genome
       5678901234    <-genome coord
            1   |
     ref_clv^   ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GAGTTGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 1), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 6))

    ref_fa = MagicMock()
    ref_fa.fetch.return_value = 'AC'
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'GACAGT'
    ref_fa.fetch.assert_called_once_with('chr2', 6, 8)


def test_extract_seq_with_skip_before_and_after_ctg_clv():
    """
             AA
           GT┘      <-bread read
       G--AGT-GC    <-bridge contig
       0  123 456   <-contig coord
     ctg_clv^   ^ice  <-contig coord
    ...GACAGTTGC... <-genome
       5678901234    <-genome coord
            1   |
     ref_clv^   ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GAGTGC'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 2),
    )

    ref_fa = MagicMock()
    ref_fa.fetch.return_value = 'AC'
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'GACAGT'
    ref_fa.fetch.assert_called_once_with('chr2', 6, 8)


def test_extract_seq_with_skip_before_and_after_ctg_clv_and_a_mismatch():
    """
    The same example as test_extract_seq_with_skip_before_and_after_ctg_clv,
    but with a mismatch introduced

             AA
           GT┘      <-bread read
       G--AGT-GC    <-bridge contig
       0  x23 456   <-contig coord
     ctg_clv^   ^ice  <-contig coord
    ...GACAGTTGC... <-genome
       5678901234    <-genome coord
            1   |
     ref_clv^   ^ire

    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GAGTGC'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 2),
    )

    ref_fa = MagicMock()
    ref_fa.fetch.return_value = 'AC'
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'GACAGT'
    ref_fa.fetch.assert_called_once_with('chr2', 6, 8)
