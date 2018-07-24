from unittest.mock import MagicMock

from kleat.evidence.suffix import is_a_suffix_read
import kleat.misc.settings as S


"""it's like testing edge cases for is_a_suffix_read in test_suffix.py.

In principle, it's highly unlikely that the suffix_read_tail_len could be
larger than suffix_contig_tail_len as the contig tail would've been extended
during the assembly stage
"""


def test_is_a_left_tail_suffix_read_for_forward_suffix_contig_with_suffix_read_tail_len_equal_to_suffix_contig_tail_len():
    """
    TTTATC                      # suffix read
    TTT
     |└ATCGC                    # suffix contig
    012345678                   # contig coord
       ^ctg_clv
    234567890                   # genome coord
       ^ref_clv
    """
    mock_read = MagicMock()
    mock_read.query_sequence = 'TTTATC'
    mock_read.reference_start = 0
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
    assert is_a_suffix_read(mock_read, mock_contig, ctg_clv=3) == 3


def test_is_a_left_tail_suffix_read_for_reverse_suffix_contig_with_suffix_read_tail_len_equal_to_suffix_contig_tail_len():
    """
           TTTATC               # suffix read
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
    ...GTAGCGATCGT...           # genome
       234567890123             # genome coord
       ref_clv^1
    """
    mock_read = MagicMock()
    mock_read.query_sequence = 'TTTATC'
    mock_read.reference_start = 0
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
    assert is_a_suffix_read(mock_read, mock_contig, ctg_clv=7) == 3


def test_is_a_right_tail_suffix_read_for_forward_suffix_contig_with_suffix_read_tail_len_equal_to_suffix_contig_tail_len():
    """
         CGCAAA                    # suffix read
            AAA
       ATCGC┘ |                    # suffix contig
       012345678                   # contig coord
           ^ctg_clv
    ...ATCGCTGC...
       234567890                   # genome coord
    ref_clv^   1
    """
    mock_read = MagicMock()
    mock_read.query_sequence = 'CGCAAA'
    mock_read.reference_end = 8
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
    assert is_a_suffix_read(mock_read, mock_contig, ctg_clv=4) == 3


def test_is_a_right_tail_suffix_read_for_reverse_suffix_contig_with_suffix_read_tail_len_equal_to_suffix_contig_tail_len():
    """
         CGCAAA                 # suffix read
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
    mock_read.query_sequence = 'CGCAAA'
    mock_read.reference_end = 11
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
    assert is_a_suffix_read(mock_read, mock_contig, ctg_clv=3) == 3


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
    assert is_a_suffix_read(mock_read, mock_contig, ctg_clv=3) is None
