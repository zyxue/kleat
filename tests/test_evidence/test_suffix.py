from unittest.mock import MagicMock

from kleat.evidence.suffix import is_a_suffix_read
import kleat.misc.settings as S


def test_is_a_suffix_read_for_forward_suffix_contig():
    """
     TT
      └ATC                      # suffix read
    TTTATCGC                    # suffix contig
    01234567
    """
    mock_read = MagicMock()
    mock_read.cigartuples = [
        (S.BAM_CSOFT_CLIP, 2),
        (S.BAM_CMATCH, 3),
    ]
    mock_read.reference_start = 3

    mock_contig = MagicMock()
    mock_contig.query_sequence = 'TTTATCGC'
    mock_contig.cigartuples = [
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 5),
    ]
    assert is_a_suffix_read(mock_read, mock_contig)


def test_is_a_suffix_read_for_reverse_suffix_contig():
    """
         AA
      CGC┘                      # suffix read
    ATCGCAAA                    # suffix contig
    01234567
    """
    mock_read = MagicMock()
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 3),
        (S.BAM_CSOFT_CLIP, 2),
    ]
    mock_read.reference_end = 4

    mock_contig = MagicMock()
    mock_contig.query_sequence = 'ATCGCAAA'
    mock_contig.infer_sequence_length.return_value = 8
    mock_contig.cigar = [
        (S.BAM_CMATCH, 5),
        (S.BAM_CSOFT_CLIP, 3),
    ]
