from unittest.mock import MagicMock, patch

import kleat.misc.settings as S
from kleat.hexamer.hexamer import extract_seq


"""
cc: ctg_clv; ice: init_clv_end
rc: ref_clv; ire: init_ref_end
"""


@patch('kleat.hexamer.search.apautils')
def test_hardclip_before_clv(mock_apautils):
    """
           AA
         TC┘         <-bridge read
    CGCATTCGTCG      <-bridge contig (hardcipped, could be chimeric https://www.biostars.org/p/109333/)
    \\\|  |          <-hardclip mask
    012345678901      <-contig coord
       |cc^    ^ice
 ...XXXATTCGTCG...   <-genome
    234567890123     <-genome coord
          | 1  |
        rc^    ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'CGCATTCGTCG'
    ctg.cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 8))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=8, ref_fa=ref_fa, ctg_clv=6)
    assert extract_seq(**kw) == 'CGCATTC'
    assert extract_seq(window=1, **kw) == 'C'
    assert extract_seq(window=2, **kw) == 'TC'
    assert extract_seq(window=3, **kw) == 'TTC'
    assert extract_seq(window=4, **kw) == 'ATTC'
    assert extract_seq(window=5, **kw) == 'CATTC'




@patch('kleat.hexamer.search.apautils')
def test_hardclip_after_clv(mock_apautils):
    """
             AAA
          GTT┘                 <-bridge read
       A-GGTTGCAGA             <-bridge contig
       | |  | |///             <-hardclip mask
       0 1234567890            <-contig coord
     ctg_clv^     ^ice <-contig coord
    ...ACGGTTGCAGA...          <-genome
       789012345678            <-genome coord
          1 |     |
     ref_clv^     ^init_fe
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    mock_apautils.infer_query_sequence.return_value = 'AGGTTGCAGA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 6),
        (S.BAM_CHARD_CLIP, 3)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch = MagicMock(return_value='C')
    kw = dict(contig=ctg, strand='+', ref_clv=12, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'ACGGTT'
    ref_fa.fetch.assert_called_with('chr1', 8, 9)

    assert extract_seq(window=1, **kw) == 'T'
    assert extract_seq(window=2, **kw) == 'TT'
    assert extract_seq(window=3, **kw) == 'GTT'
    assert extract_seq(window=4, **kw) == 'GGTT'
    assert extract_seq(window=5, **kw) == 'CGGTT'


@patch('kleat.hexamer.search.apautils')
def test_hardclip_spanning_clv_from_before_edgecase_1(mock_apautils):
    """
           AA
         TC┘         <-bridge read
      CATTCGT        <-bridge contig (hardcipped, could be chimeric https://www.biostars.org/p/109333/)
      \\\\| |        <-hardclip mask
      01234567       <-contig coord
        cc^ |^ice
   ...XATTCGT...     <-genome
      23456789       <-genome coord
          |  |
        rc^  ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'CATTCGT'
    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 4),
        (S.BAM_CMATCH, 3)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=6, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'CATTC'


@patch('kleat.hexamer.search.apautils')
def test_hardclip_spanning_clv_from_before_edgecase_2(mock_apautils):
    """
           AA
         TC┘         <-bridge read
      CATTCGT        <-bridge contig (hardcipped, could be chimeric https://www.biostars.org/p/109333/)
      \\\\\ |        <-hardclip mask
      01234567       <-contig coord
        cc^  ^ice
   ...XATTCGT...     <-genome
      23456789       <-genome coord
          |  |
        rc^  ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'CATTCGT'
    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 5),
        (S.BAM_CMATCH, 2)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=6, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == ''



@patch('kleat.hexamer.search.apautils')
def test_hardclip_spanning_clv_from_before_edgecase_3(mock_apautils):
    """
           AA
         TC┘         <-bridge read
      CATTCGT        <-bridge contig (hardcipped, could be chimeric https://www.biostars.org/p/109333/)
      \\\\\\|        <-hardclip mask
      01234567       <-contig coord
        cc^  ^ice
   ...XATTCGT...     <-genome
      23456789       <-genome coord
          |  |
        rc^  ^ire
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr2'
    mock_apautils.infer_query_sequence.return_value = 'CATTCGT'
    ctg.cigartuples = (
        (S.BAM_CHARD_CLIP, 6),
        (S.BAM_CMATCH, 1)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='+', ref_clv=6, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == ''



@patch('kleat.hexamer.search.apautils')
def test_hardclip_spanning_clv_from_after_edgecase_1(mock_apautils):
    """
             AAA
          GTT┘           <-bridge read
       A-GGTTGCA         <-bridge contig
       | |  |///         <-hardclip mask
       0 12345678        <-contig coord
          cc^   ^ice
    ...ACGGTTGCA...      <-genome
       7890123456        <-genome coord
          1 |   |
          rc^   ^ie
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    mock_apautils.infer_query_sequence.return_value = 'AGGTTGCA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 4),
        (S.BAM_CHARD_CLIP, 3)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch = MagicMock(return_value='C')
    kw = dict(contig=ctg, strand='+', ref_clv=12, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'ACGGTT'
    ref_fa.fetch.assert_called_with('chr1', 8, 9)

    assert extract_seq(window=1, **kw) == 'T'
    assert extract_seq(window=2, **kw) == 'TT'
    assert extract_seq(window=3, **kw) == 'GTT'
    assert extract_seq(window=4, **kw) == 'GGTT'
    assert extract_seq(window=5, **kw) == 'CGGTT'



@patch('kleat.hexamer.search.apautils')
def test_hardclip_spanning_clv_from_after_edgecase_2(mock_apautils):
    """
             AAA
          GTT┘           <-bridge read
       A-GGTTGCA         <-bridge contig
       | |  ////         <-hardclip mask
       0 12345678        <-contig coord
          cc^   ^ice
    ...ACGGTTGCA...      <-genome
       7890123456        <-genome coord
          1 |   |
          rc^   ^ie
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    mock_apautils.infer_query_sequence.return_value = 'AGGTTGCA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 3),
        (S.BAM_CHARD_CLIP, 4),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch = MagicMock(return_value='C')
    kw = dict(contig=ctg, strand='+', ref_clv=12, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == ''



@patch('kleat.hexamer.search.apautils')
def test_hardclip_spanning_clv_from_after_edgecase_3(mock_apautils):
    """
             AAA
          GTT┘           <-bridge read
       A-GGTTGCA         <-bridge contig
       | | /////         <-hardclip mask
       0 12345678        <-contig coord
          cc^   ^ice
    ...ACGGTTGCA...      <-genome
       7890123456        <-genome coord
          1 |   |
          rc^   ^ie
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    mock_apautils.infer_query_sequence.return_value = 'AGGTTGCA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 2),
        (S.BAM_CHARD_CLIP, 5),
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch = MagicMock(return_value='C')
    kw = dict(contig=ctg, strand='+', ref_clv=12, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == ''
