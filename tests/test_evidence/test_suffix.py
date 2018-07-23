import pytest
from unittest.mock import MagicMock

from kleat.evidence import suffix
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
    mock_read.query_sequence = 'TTATC'
    mock_read.reference_start = 3
    mock_read.cigartuples = [
        (S.BAM_CSOFT_CLIP, 2),
        (S.BAM_CMATCH, 3),
    ]

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
    mock_read.query_sequence = 'CGCAA'
    mock_read.reference_end = 4
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 3),
        (S.BAM_CSOFT_CLIP, 2),
    ]

    mock_contig = MagicMock()
    mock_contig.query_sequence = 'ATCGCAAA'
    mock_contig.infer_query_length.return_value = 8
    mock_contig.cigartuples = [
        (S.BAM_CMATCH, 5),
        (S.BAM_CSOFT_CLIP, 3),
    ]
    assert is_a_suffix_read(mock_read, mock_contig)


def test_is_not_a_suffix_read_for_forward_suffix_contig():
    """
       ATC                      # suffix read
    TTTATCGC                    # suffix contig
    01234567
    """
    mock_read = MagicMock()
    mock_read.reference_start = 3
    mock_read.query_sequence = 'ATC'
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 3),
    ]

    mock_contig = MagicMock()
    mock_contig.query_sequence = 'TTTATCGC'
    mock_contig.cigartuples = [
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 5),
    ]
    assert not is_a_suffix_read(mock_read, mock_contig)


def test_is_not_a_suffix_read_for_reverse_suffix_contig():
    """
      CGC                       # suffix read
    ATCGCAAA                    # suffix contig
    01234567
    """
    mock_read = MagicMock()
    mock_read.query_sequence = 'CGC'
    mock_read.reference_end = 4
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 3),
    ]

    mock_contig = MagicMock()
    mock_contig.query_sequence = 'ATCGCAAA'
    mock_contig.infer_query_length.return_value = 8
    mock_contig.cigartuples = [
        (S.BAM_CMATCH, 5),
        (S.BAM_CSOFT_CLIP, 3),
    ]
    assert not is_a_suffix_read(mock_read, mock_contig)


# At the momment, these variables are just copied form
# test_apautils_tail_related.py, could modify at convenience later
mock_non_tailed_read = MagicMock()
mock_non_tailed_read.query_sequence = "ATCG"
mock_non_tailed_read.cigartuples = [(S.BAM_CMATCH, 4)]

#  TTT
#  ||└AT
# 789012 <- coord (0-based)
#    1
mock_left_tailed_read = MagicMock()
mock_left_tailed_read.query_sequence = "TTTAT"
mock_left_tailed_read.cigartuples = [(S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 2)]
# start is set to one base from A, consistent with pysam + IGV observation
mock_left_tailed_read.reference_start = 11  # this is artbitray
mock_left_tailed_read.reference_end = 13

#    AAAA
#  AT┘|||
# 0123456 <- coord
# 1
mock_right_tailed_read = MagicMock()
mock_right_tailed_read.query_sequence = "ATAAAA"
mock_right_tailed_read.cigartuples = [(S.BAM_CMATCH, 2), (S.BAM_CSOFT_CLIP, 4)]
mock_right_tailed_read.reference_start = 11  # this is arbitray
mock_right_tailed_read.reference_end = 13


def test_calc_ref_clv():
    f = suffix.calc_ref_clv
    assert f(mock_left_tailed_read) == 11
    assert f(mock_right_tailed_read) == 12
    with pytest.raises(ValueError, match=r'not a suffix segment'):
        f(mock_non_tailed_read)


def test_calc_ref_clv_passing_tail_side_argument():
    f = suffix.calc_ref_clv
    assert f(mock_left_tailed_read, 'left') == 11
    assert f(mock_right_tailed_read, 'right') == 12


def test_calc_strand_passing_tail_side_argument():
    f = suffix.calc_strand
    assert f('left') == '-'
    assert f('right') == '+'
    err_msg = 'tail_side must be "left" or "right", but None passed'
    with pytest.raises(ValueError, match=err_msg):
        f(None)
