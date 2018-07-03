from unittest.mock import MagicMock

from kleat.misc.settings import BAM_CSOFT_CLIP, BAM_CMATCH
from kleat.misc import apautils

mock_non_tailed_read = MagicMock()
mock_non_tailed_read.query_sequence = "ATCGATCG"
mock_non_tailed_read.cigartuples = [(BAM_CMATCH, 8)]

mock_left_tailed_read = MagicMock()
mock_left_tailed_read.query_sequence = "TTTTATCG"
mock_left_tailed_read.cigartuples = [(BAM_CSOFT_CLIP, 4), (BAM_CMATCH, 4)]

mock_right_tailed_read = MagicMock()
mock_right_tailed_read.query_sequence = "ATCGAAAA"
mock_right_tailed_read.cigartuples = [(BAM_CMATCH, 4), (BAM_CSOFT_CLIP, 4)]

assert apautils.left_tail(mock_left_tailed_read) is True
assert apautils.right_tail(mock_left_tailed_read) is False
assert apautils.left_tail(mock_non_tailed_read) is False

assert apautils.right_tail(mock_left_tailed_read) is False
assert apautils.right_tail(mock_right_tailed_read) is True
assert apautils.right_tail(mock_non_tailed_read) is False
