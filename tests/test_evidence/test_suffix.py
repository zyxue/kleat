import pytest
from unittest.mock import MagicMock

from kleat.evidence import suffix
from kleat.evidence.suffix import is_a_suffix_read
import kleat.misc.settings as S


def test_is_a_left_tail_suffix_read_for_forward_suffix_contig():
    """
     TTATC                      # suffix read
    TTT
     |└ATCGC                    # suffix contig
    012345678                   # contig coord
       ^ctg_clv
    234567890                   # genome coord
       ^ref_clv
    """
    mock_read = MagicMock()
    mock_read.query_sequence = 'TTATC'
    mock_read.reference_start = 1
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 5),
    ]

    mock_contig = MagicMock()
    mock_contig.is_reverse = False
    mock_contig.query_sequence = 'TTTATCGC'
    mock_contig.cigartuples = [
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 5),
    ]
    assert is_a_suffix_read(mock_read, mock_contig, ctg_clv=3)


def test_is_a_left_tail_suffix_read_for_reverse_suffix_contig():
    """
    Actually, it makes no difference to is_a_suffix_read whether contig is reversed or not

            TTATC               # suffix read
           TTT
           | └ATCGCTAC          # suffix contig in r2c alignment
           012345678901         # contig coord in r2c alignment
              |                 # flip
              |AAA
       GTAGCGAT┘ |              # suffix contig in c2g alignment
       012345678901             # contig coord in c2g alignment
       ctg_clv^  1
      109876543210              # rev contig coord, corresponding to contig coord in r2c alignment
       1      |ctg_clv
    ...GTAGCGATCGT...            # genome
       234567890123             # genome coord
       ref_clv^1
    """
    mock_read = MagicMock()
    mock_read.query_sequence = 'TTATC'
    mock_read.reference_start = 1
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 5),
    ]

    mock_contig = MagicMock()
    mock_contig.is_reverse = True
    mock_contig.query_sequence = 'GTAGCGGCGATAAA'
    mock_contig.infer_query_length.return_value = 11
    mock_contig.cigartuples = [
        (S.BAM_CMATCH, 5),
        (S.BAM_CSOFT_CLIP, 3),
    ]
    assert is_a_suffix_read(mock_read, mock_contig, ctg_clv=7)


def test_is_a_right_tail_suffix_read_for_forward_suffix_contig():
    """
         CGCAA                      # suffix read
            AAA
       ATCGC┘ |                   # suffix contig
       012345678                   # contig coord
           ^ctg_clv
    ...ATCGCTGC...
       234567890                   # genome coord
    ref_clv^   1
    """
    mock_read = MagicMock()
    mock_read.query_sequence = 'CGCAA'
    mock_read.reference_end = 6
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 5),
    ]

    mock_contig = MagicMock()
    mock_contig.is_reverse = False
    mock_contig.query_sequence = 'ATCGCAAA'
    mock_contig.cigartuples = [
        (S.BAM_CMATCH, 5),
        (S.BAM_CSOFT_CLIP, 3),
    ]
    assert is_a_suffix_read(mock_read, mock_contig, ctg_clv=4)


def test_is_a_right_tail_suffix_read_for_reverse_suffix_contig():
    """
         CGCAA                  # suffix read
            AAA
    GCGATCGC┘ |                 # suffix contig in r2c alignment
    012345678901                # contig coord in r2c alignment
           |  1                 # flip
        TTTGCGTACGC             # suffix contig in c2g alignment
        012345678901            # contig coord in c2g alignment
    ctg_clv^      1
       109876543210             # rev contig coord, corresponding to contig coord in r2c alignment
        1  |
        ...GCGTACGC...          # genome
           23456789             # genome coord
           ^ref_clv
    """
    mock_read = MagicMock()
    mock_read.query_sequence = 'CGCAA'
    mock_read.reference_end = 9
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 5),
    ]

    mock_contig = MagicMock()
    mock_contig.is_reverse = True
    mock_contig.query_sequence = 'TTTGCGTACGC'
    mock_contig.infer_query_length.return_value = 11
    mock_contig.cigartuples = [
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 8),
    ]
    assert is_a_suffix_read(mock_read, mock_contig, ctg_clv=3)


def test_is_not_a_suffix_read_for_forward_suffix_contig():
    """
       ATC                      # suffix read
    TTT|
      └ATCGC                    # suffix contig
    01234567
       ^ctg_clv
    """
    mock_read = MagicMock()
    mock_read.reference_start = 3
    mock_read.query_sequence = 'ATC'
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 3),
    ]

    mock_contig = MagicMock()
    mock_contig.is_reverse = False
    mock_contig.query_sequence = 'TTTATCGC'
    mock_contig.cigartuples = [
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 5),
    ]
    assert not is_a_suffix_read(mock_read, mock_contig, ctg_clv=3)


def test_is_not_a_suffix_read_for_reverse_suffix_contig():
    """
      CGC                     # suffix read
        |AAA
    ATCGC┘                    # suffix contig
    01234567
        ^ctg_clv
    """
    mock_read = MagicMock()
    mock_read.query_sequence = 'CGC'
    mock_read.reference_end = 4
    mock_read.cigartuples = [
        (S.BAM_CMATCH, 3),
    ]

    mock_contig = MagicMock()
    mock_contig.is_reverse = True
    mock_contig.query_sequence = 'ATCGCAAA'
    mock_contig.infer_query_length.return_value = 8
    mock_contig.cigartuples = [
        (S.BAM_CMATCH, 5),
        (S.BAM_CSOFT_CLIP, 3),
    ]
    assert not is_a_suffix_read(mock_read, mock_contig, ctg_clv=4)


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
