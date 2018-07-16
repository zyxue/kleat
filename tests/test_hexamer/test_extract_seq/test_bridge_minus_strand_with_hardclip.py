from unittest.mock import MagicMock, patch

import kleat.misc.settings as S
from kleat.hexamer.search import extract_seq


@patch('kleat.hexamer.search.apautils')
def test_extract_seq_with_hardclipped_region(mock_apautils):
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
    mock_apautils.infer_query_sequence.return_value = 'ATTCG'
    ctg.cigartuples = ((S.BAM_CMATCH, 5), (S.BAM_CHARD_CLIP, 3))

    ref_fa = MagicMock()
    ref_fa.get_reference_length.return_value = 100
    kw = dict(contig=ctg, strand='-', ref_clv=9, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'G'
