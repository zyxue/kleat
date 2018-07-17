from unittest.mock import MagicMock, patch

import kleat.misc.settings as S
from kleat.hexamer.hexamer import extract_seq


"""
cc: ctg_clv; icb: init_clv_beg
rc: ref_clv; irb: init_ref_beg
"""


@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_after_clv(mock_apautils):
    """
             TT
             |└GTGA     <-bridge read
           ATTCGTGA     <-bridge contig (hardcipped), chimeric
                 //     <-hardclip mask
           012345678    <-contig coord
        icb^   ^cc
        ...ATTCGTGA...
           567890123    <-genome coordinate
           |   |1
        irb^   ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGTGA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 6),
        (S.BAM_CHARD_CLIP, 2),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'GTGA'
    assert extract_seq(window=5, **kw) == 'GTGA'


@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_before_clv(mock_apautils):
    """
             TT
             |└GTGA     <-bridge read
           ATTCGTGA     <-bridge contig (hardcipped), chimeric
           //           <-hardclip mask
           012345678    <-contig coord
        icb^   ^cc
        ...ATTCGTGA...
           567890123    <-genome coordinate
           |   |1
        irb^   ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGTGA'
    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 2),
        (S.BAM_CMATCH, 6),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'GTGA'


@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_spanning_clv_from_before(mock_apautils):
    """
             TT
             |└GTGA     <-bridge read
           ATTCGTGA     <-bridge contig (hardcipped), chimeric
           //////       <-hardclip mask
           012345678    <-contig coord
        icb^   ^cc
        ...ATTCGTGA...
           567890123    <-genome coordinate
           |   |1
        irb^   ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGTGA'
    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 6),
        (S.BAM_CMATCH, 2),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == ''


@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_spanning_clv_from_before_edgecase_1(mock_apautils):
    """
             TT
             |└GTGA     <-bridge read
           ATTCGTGA     <-bridge contig (hardcipped), chimeric
           /////        <-hardclip mask
           012345678    <-contig coord
        icb^   ^cc
        ...ATTCGTGA...
           567890123    <-genome coordinate
           |   |1
        irb^   ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGTGA'
    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 5),
        (S.BAM_CMATCH, 3),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == ''



@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_not_spanning_clv_from_before_edgecase_2(mock_apautils):
    """
             TT
             |└GTGA     <-bridge read
           ATTCGTGA     <-bridge contig (hardcipped), chimeric
           ////         <-hardclip mask
           012345678    <-contig coord
        icb^   ^cc
        ...ATTCGTGA...
           567890123    <-genome coordinate
           |   |1
        irb^   ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGTGA'
    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 4),
        (S.BAM_CMATCH, 4),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'GTGA'


@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_spanning_clv_from_before_edgecase_3(mock_apautils):
    """
             TT
              └GTGA     <-bridge read
               GTGA     <-bridge contig (hardcipped), chimeric
               /        <-hardclip mask
               012345678    <-contig coord
               ^icb/cc
            ...GTGA...
               0123    <-genome coordinate
               |1
               ^irb/rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'GTGA'
    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 1),
        (S.BAM_CMATCH, 3),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=0, ref_fa=ref_fa, ctg_clv=0)
    assert extract_seq(**kw) == ''


@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_spanning_clv_from_after_edgecase_1(mock_apautils):
    """
             TT
             |└GTGA     <-bridge read
           ATTCGTGA     <-bridge contig (hardcipped), chimeric
               ////     <-hardclip mask
           012345678    <-contig coord
        icb^   ^cc
        ...ATTCGTGA...
           567890123    <-genome coordinate
           |   |1
        irb^   ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGTGA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 4),
        (S.BAM_CHARD_CLIP, 4),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == ''


@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_not_spanning_clv_from_after_edgecase_2(mock_apautils):
    """
             TT
             |└GTGA     <-bridge read
           ATTCGTGA     <-bridge contig (hardcipped), chimeric
                ///     <-hardclip mask
           012345678    <-contig coord
        icb^   ^cc
        ...ATTCGTGA...
           567890123    <-genome coordinate
           |   |1
        irb^   ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGTGA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 5),
        (S.BAM_CHARD_CLIP, 3),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'GTGA'


@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_spanning_clv_from_after_edgecase_3(mock_apautils):
    """
             TT
             |└GTGA     <-bridge read
           ATTCGTGA     <-bridge contig (hardcipped), chimeric
              /////     <-hardclip mask
           012345678    <-contig coord
        icb^   ^cc
        ...ATTCGTGA...
           567890123    <-genome coordinate
           |   |1
        irb^   ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGTGA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 3),
        (S.BAM_CHARD_CLIP, 5),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == ''


@patch('kleat.hexamer.search.apautils')
def test_with_hardclip_spanning_clv_from_after(mock_apautils):
    """
             TT
             |└GTGA     <-bridge read
           ATTCGTGA     <-bridge contig (hardcipped), chimeric
             //////     <-hardclip mask
           012345678    <-contig coord
        icb^   ^cc
        ...ATTCGTGA...
           567890123    <-genome coordinate
           |   |1
        irb^   ^rc
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'ATTCGTGA'
    ctg.cigartuples = ((S.BAM_CMATCH, 2), (S.BAM_CHARD_CLIP, 6))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == ''
