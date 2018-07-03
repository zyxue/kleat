import pytest

from unittest.mock import MagicMock

from kleat.misc.settings import BAM_CSOFT_CLIP, BAM_CMATCH
from kleat.misc import apautils as A

mock_non_tailed_read = MagicMock()
mock_non_tailed_read.query_sequence = "ATCG"
mock_non_tailed_read.cigartuples = [(BAM_CMATCH, 4)]

mock_left_tailed_read = MagicMock()
mock_left_tailed_read.query_sequence = "TTTAT"
mock_left_tailed_read.cigartuples = [(BAM_CSOFT_CLIP, 3), (BAM_CMATCH, 2)]
# start is set to one base from A, consistent with pysam + IGV observation
mock_left_tailed_read.reference_start = 1
mock_left_tailed_read.reference_end = 3

mock_right_tailed_read = MagicMock()
mock_right_tailed_read.query_sequence = "ATAAAA"
mock_right_tailed_read.cigartuples = [(BAM_CMATCH, 2), (BAM_CSOFT_CLIP, 4)]
mock_right_tailed_read.reference_start = 1
mock_right_tailed_read.reference_end = 3


def test_left_tail():
    assert A.left_tail(mock_left_tailed_read) is True
    assert A.left_tail(mock_right_tailed_read) is False
    assert A.left_tail(mock_non_tailed_read) is False


def test_right_tail():
    assert A.right_tail(mock_left_tailed_read) is False
    assert A.right_tail(mock_right_tailed_read) is True
    assert A.right_tail(mock_non_tailed_read) is False


def test_calc_ref_clv():
    assert A.calc_ref_clv(mock_left_tailed_read) == 2
    assert A.calc_ref_clv(mock_right_tailed_read) == 3
    with pytest.raises(ValueError, match=r'not a suffix segment'):
        A.calc_ref_clv(mock_non_tailed_read)


def test_calc_ref_clv_passing_tail_side_argument():
    assert A.calc_ref_clv(mock_left_tailed_read, 'left') == 2
    assert A.calc_ref_clv(mock_right_tailed_read, 'right') == 3


def test_calc_tail_length():
    assert A.calc_tail_length(mock_left_tailed_read) == 3
    assert A.calc_tail_length(mock_right_tailed_read) == 4
    with pytest.raises(ValueError, match=r'not a suffix segment'):
        A.calc_ref_clv(mock_non_tailed_read)


def test_calc_tail_length_passing_tail_side_argument():
    assert A.calc_tail_length(mock_left_tailed_read, 'left') == 3
    assert A.calc_tail_length(mock_right_tailed_read, 'right') == 4
    err_msg = r'not be a right tailed segment as its last CIGAR is not BAM_CSOFT_CLIP'
    with pytest.raises(ValueError, match=err_msg):
        A.calc_tail_length(mock_left_tailed_read, 'right')
