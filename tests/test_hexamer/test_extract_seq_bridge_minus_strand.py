from unittest.mock import MagicMock, call

import kleat.misc.settings as S
from kleat.hexamer.search import extract_seq


"""
cc: ctg_clv; ici: init_clv_idx
rc: ref_clv; iri: init_ref_idx
"""


def test_extract_seq_for_bridge():
    """
        T
        └CG         <-bread read
       GACGGTTGC    <-bridge contig
       0123456789   <-contig coord
    ici^ ^ctg_clv   <-contig coord
    ...GACGGTTGC... <-genome
       567890123    <-genome coord
       | |  1
    iri^ ^ref_clv
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GACGGTTGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 9),)

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=7, ref_fa=ref_fa, ctg_clv=2)
    assert extract_seq(**kw) == 'CGGTTGC'


def test_extract_seq_for_bridge_with_skip():
    """
        T
        └GT         <-bread read
       GACGGT-GC    <-bridge contig
       012345 678   <-contig coord
    ici^ ^ctg_clv   <-contig coord
    ...GACGGTAGC... <-genome
       5678901234   <-genome coord
       | |  1
    iri^ ^ref_clv
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GACGGTGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 6), (S.BAM_CREF_SKIP, 1), (S.BAM_CMATCH, 2))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch.return_value = 'A'

    kw = dict(contig=ctg, strand='-', ref_clv=7, ref_fa=ref_fa, ctg_clv=2)
    assert extract_seq(**kw) == 'CGGTAGC'
    ref_fa.fetch.assert_called_with('chr2', 11, 12)


def test_extract_seq_for_bridge_with_skip_before_clv():
    """
             T
             └AG       <-bridge read
       GA--GGTAGC      <-bridge contig
       01  2345678     <-contig coord
    ici^ | |  ^ctg_clv <-contig coord
    ...GACTGGTAGC...   <-genome
       5678901234      <-genome coord
       |    1 |
    iri^      ^ref_clv
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GAGGTAGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 2), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 6))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=12, ref_fa=ref_fa, ctg_clv=5)
    assert extract_seq(**kw) == 'AGC'


def test_extract_seq_for_bridge_with_multiple_skips_before_clv():
    """
             T
             └AG       <-bridge read
       G-C--GTAGC      <-bridge contig
       0 1  234567     <-contig coord
    ici^||| | ^ctg_clv <-contig coord
    ...GACTGGTAGC...   <-genome
       5678901234      <-genome coord
       |    1 |
    iri^      ^ref_clv
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GCGTAGC'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 5)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=12, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'AGC'


def test_extract_seq_for_bridge_with_deletion():
    """
        T
        └GT         <-bread read
       GACGGT_CGC    <-bridge contig
       012345 678   <-contig coord
    ici^ ^cc | x 1   <-contig coord
    ...GACGGTCCTC... <-genome
       56789012345    <-genome coord
       | |  1
    iri^ ^ref_clv
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GACGGTCGC'
    ctg.cigartuples = ((S.BAM_CMATCH, 6), (S.BAM_CDEL, 1), (S.BAM_CMATCH, 3))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100

    kw = dict(contig=ctg, strand='-', ref_clv=7, ref_fa=ref_fa, ctg_clv=2)
    assert extract_seq(**kw) == 'CGGTCGC'


def test_extract_seq_for_bridge_with_insertion():
    """
        T    AG      <-inserted bases
        └CG  ┬       <-bread read
       GACGGT CTC    <-bridge contig
       012345 8901   <-contig coord
    ici^ ^cc   x1    <-contig coord
    ...GACGGT CGC... <-genome
       567890 1234    <-genome coord
       | |  1
    iri^ ^ref_clv
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'GACGGTAGCTC'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 6),
        (S.BAM_CINS, 2),
        (S.BAM_CMATCH, 3)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100

    kw = dict(contig=ctg, strand='-', ref_clv=7, ref_fa=ref_fa, ctg_clv=2)
    assert extract_seq(**kw) == 'CGGTAGCTC'


def test_extract_seq_with_hardclipped_region():
    """
             TT
             |└GXXX  <-bridge read
           ATTCG///  <-bridge contig (hardcipped), chimeric
           012345678 <-contig coord
               ^ctg_clv/ici
        ...ATTCGXXX...
           567890123 <-genome coordinate
               |1
               ^ref_clv/iri
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    ctg.query_sequence = 'ATTCG'
    ctg.cigartuples = ((S.BAM_CMATCH, 5), (S.BAM_CHARD_CLIP, 3))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'G'
