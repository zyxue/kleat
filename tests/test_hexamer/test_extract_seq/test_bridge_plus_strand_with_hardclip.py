from unittest.mock import MagicMock, patch

import kleat.misc.settings as S
from kleat.hexamer.search import extract_seq


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


@patch('kleat.hexamer.search.apautils')
def test_extract_seq_with_hardclipping(mock_apautils):
    """
             AAA             <-bridge read contig
          GTT┘
       A-GGTTGCAGA             <-suffix contig
       | |  | |///             <-hardclip mask
       012345678901            <-contig coord
            |    1
     ctg_clv^    ^init_ctg_idx <-contig coord
    ...ACGGTTGCAGA...          <-genome
       789012345678            <-genome coord
          1 |    |
     ref_clv^    ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    mock_apautils.infer_query_sequence.return_value = 'AGGGCAGA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 6),
        (S.BAM_CHARD_CLIP, 3)
    )

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    ref_fa.fetch = MagicMock(return_value='C')
    kw = dict(contig=ctg, strand='+', ref_clv=12, ref_fa=ref_fa, ctg_clv=5)
    assert extract_seq(**kw) == 'AGGTT'
    # assert extract_seq(window=3, **kw) == 'TGC'
    # ref_fa.fetch.assert_called_with('chr1', 11, 13)
