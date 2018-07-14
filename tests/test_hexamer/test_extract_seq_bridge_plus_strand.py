from unittest.mock import MagicMock, patch

import kleat.misc.settings as S
from kleat.hexamer.search import extract_seq


"""
cc: ctg_clv; ice: init_clv_e
rc: ref_clv; ire: init_ref_end
"""


def test_extract_seq():
    """
             AA
           GT┘       <-bridge read
       GACGGTTGC     <-bridge contig
       0123456789    <-contig coord
     ctg_clv^   ^ice
    ...GACGGTTGC...  <-genome
       5678901234    <-genome coord
            1   |
     ref_clv^   ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GACGGTTGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 9),)

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=5)
    assert extract_seq(**kw) == 'GACGGT'


def test_extract_seq_with_skip_after_ctg_clv():
    """
             AA
           GT┘       <-bridge read
       GACGGT-GC     <-bridge contig
       012345 678    <-contig coord
     ctg_clv^   ^ice
    ...GACGGTTGC...  <-genome
       5678901234    <-genome coord
            1   |
     ref_clv^   ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GACGGTGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 6), (S.BAM_CREF_SKIP, 1), (S.BAM_CMATCH, 2))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch.return_value = 'T'
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=5)
    assert extract_seq(**kw) == 'GACGGT'
    assert ref_fa.call_count == 0


def test_extract_seq_with_skip_before_ctg_clv():
    """
             AA
           GT┘       <-bridge read
       G--AGTTGC     <-bridge contig
       0  1234567    <-contig coord
     ctg_clv^   ^ice
    ...GACAGTTGC...  <-genome
       5678901234    <-genome coord
            1   |
     ref_clv^   ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GAGTTGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 1), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 6))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch.return_value = 'AC'
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'GACAGT'
    ref_fa.fetch.assert_called_once_with('chr2', 6, 8)


def test_extract_seq_with_skip_before_and_after_ctg_clv():
    """
             AA
           GT┘       <-bridge read
       G--AGT-GC     <-bridge contig
       0  123 456    <-contig coord
     ctg_clv^   ^ice
    ...GACAGTTGC...  <-genome
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
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch.return_value = 'AC'
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'GACAGT'
    ref_fa.fetch.assert_called_once_with('chr2', 6, 8)


def test_extract_seq_with_skip_before_and_after_ctg_clv_and_a_mismatch():
    """
    The same example as test_extract_seq_with_skip_before_and_after_ctg_clv,
    but with a mismatch introduced

             AA
           GT┘        <-bridge read
       G--AGT-GC      <-bridge contig
       0  x23 456     <-contig coord
     ctg_clv^   ^ice
    ...GACAGTTGC...   <-genome
       5678901234     <-genome coord
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
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch.return_value = 'AC'
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'GACAGT'
    ref_fa.fetch.assert_called_once_with('chr2', 6, 8)


def test_extract_seq_for_bridge_with_multiple_skips_before_clv():
    """
               AA
             TA┘       <-bridge read
       G-C--CTAGC      <-bridge contig
       0 1  234567     <-contig coord
        ||| x ^cc^ice
    ...GACTGGTAGC...   <-genome
       56789012345     <-genome coord
            1 |  |
            rc^  ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GCCTAGC'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 5)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch.side_effect = list(reversed(['A', 'TG']))  # from Right => Left
    kw = dict(contig=ctg, strand='+', ref_clv=12, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'GACTGCTA'


def test_extract_seq_for_bridge_with_deletion():
    """
               AA
             CG┘      <-bridge read
       GAC__TCGTC     <-bridge contig
       012  345678    <-contig coord
          | ||x
          | |^cc ^ice
    ...GACGGTCCTC...  <-genome
       56789012345    <-genome coord
            1|   |
             ^rc ^rce
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GACTCGTC'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 3),
        (S.BAM_CDEL, 2),
        (S.BAM_CMATCH, 5),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=15, ref_fa=ref_fa, ctg_clv=5)
    assert extract_seq(**kw) == 'GACTCG'


def test_extract_seq_for_bridge_with_insertion():
    """
         AG   AA      <-inserted bases
         ┬  GT┘       <-bread read
       GA CGGTCGC     <-bridge contig
       01 45678901    <-contig coord
        x   1|   |
        x  cc^   ^ice
    ...GT CGGTCGC...  <-genome
       56 78901234    <-genome coord
             |   |
           rc^   ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GAAGCGGTCGC'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 2),
        (S.BAM_CINS, 2),
        (S.BAM_CMATCH, 7)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=10, ref_fa=ref_fa, ctg_clv=7)
    assert extract_seq(**kw) == 'GAAGCGGT'


@patch('kleat.hexamer.search.apautils')
def test_extract_seq_with_hardclipped_region(mock_apautils):
    """
           AA
         TC┘|      <-bridge read
    \\\ATTCGT      <-bridge contig (hardcipped, could be chimeric https://www.biostars.org/p/109333/)
       0123456     <-contig coord
        cc^  ^ice
    ...ATTCGXXX... <-genome
       567890123   <-genome coord
          | 1|
        rc^  ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGT'
    ctg.cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 6))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=8, ref_fa=ref_fa, ctg_clv=3)
    assert extract_seq(**kw) == 'ATTC'
