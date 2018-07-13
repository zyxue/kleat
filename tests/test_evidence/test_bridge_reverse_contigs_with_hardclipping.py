from collections import defaultdict
from unittest.mock import MagicMock

from kleat.evidence import bridge
from kleat.misc.apautils import gen_clv_key_tuple
import kleat.misc.settings as S


###################################################
# test different situations for do_rev_ctg_lt_bdg #
###################################################

def test_do_rev_ctg_lt_bdg_with_left_hard_clipping():
    """
    the hardclip won't have an effect in such case

         TT
          └ACG         <-left-tail read
         \\ACGXX       <-contig
         0123456       <-contig coord
           ^ctg_offset
        ...ACGXX...    <-rev ref genome
          543210       <-rev ref genome offset
          876543       <-rev ref coord
    ref_clv^   ^begining of contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 2
    read.reference_end = 5
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 3))

    contig = MagicMock()
    contig.infer_query_length.return_value = 7  # including hardclip
    contig.cigartuples = (
        # note from right to left
        (S.BAM_CMATCH, 5),
        (S.BAM_CHARD_CLIP, 2),
    )

    ctg_offset = 4
    tail_len = 2
    assert bridge.do_rev_ctg_lt_bdg(read, contig) == ('+', ctg_offset, tail_len)


def test_do_rev_ctg_lt_bdg_with_left_hard_clipping_passing_ctg_clv():
    """
    TTT
      └ACG      <-left-tail read
      \\CGX     <-contig, the hardclipped part would appear in another position in genome, although the read is well aligned to it
      012345    <-contig coord
       ^ctg_offset (won't be captured by this read)
     ...CGX...  <-rev ref genome
       3210     <-rev ref genome offset
       6543     <-rev ref genome coord
          ^begining of contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 1
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))

    contig = MagicMock()
    contig.infer_query_length.return_value = 6  # including hardclip
    contig.cigartuples = (
        # note from right to left
        (S.BAM_CMATCH, 3),
        (S.BAM_CHARD_CLIP, 2),
    )

    assert bridge.do_rev_ctg_lt_bdg(read, contig) is None


def test_do_rev_ctg_lt_bdg_with_left_hard_clipping_right_on_ctg_clv():
    """
        TTT
          └ACG      <-left-tail read
          \ACGX     <-contig, the hardclipped part would appear in another position in genome, although the read is well aligned to it
          012345    <-contig coord
           ^ctg_offset (won't be captured by this read)
         ...CGX...  <-rev ref genome
           3210     <-rev ref genome offset
           6543     <-rev ref genome coord
    ref_clv^  ^begining of contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 1
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))

    contig = MagicMock()
    contig.infer_query_length.return_value = 5  # including hardclip
    contig.cigartuples = (
        # note from right to left
        (S.BAM_CMATCH, 4),
        (S.BAM_CHARD_CLIP, 1),
    )

    ctg_offset = 3
    tail_len = 3
    assert bridge.do_rev_ctg_lt_bdg(read, contig) == ('+', ctg_offset, tail_len)


def test_do_rev_ctg_lt_bdg_with_right_hard_clipping_passing_ctg_clv():
    """
       TTT
         └AG      <-left-tail read
       XX\\\\     <-contig, the hardclipped part would appear in another position in genome, although the read is well aligned to it
       0123456    <-contig coord
          ^ctg_offset (won't be captured by this read)
    ...XX...      <-rev ref genome
       10         <-rev ref genome offset
       65         <-rev ref genome coord
        ^begining of contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 3
    read.reference_end = 5
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 2))

    contig = MagicMock()
    contig.infer_query_length.return_value = 6  # including hardclip
    contig.cigartuples = (
        # note from right to left
        (S.BAM_CHARD_CLIP, 4),
        (S.BAM_CMATCH, 2),
    )
    assert sum(_[1] for _ in contig.cigartuples) == contig.infer_query_length()
    assert bridge.do_rev_ctg_lt_bdg(read, contig) is None


def test_do_rev_ctg_lt_bdg_with_right_hard_clipping_right_on_ctg_clv():
    """
       TTT
         └AG      <-left-tail read
       XXXA\\     <-contig, the hardclipped part would appear in another position in genome, although the read is well aligned to it
       0123456    <-contig coord
          ^ctg_offset (won't be captured by this read)
    ...XXXA...    <-rev ref genome
       3210       <-rev ref genome offset
       6543       <-rev ref genome coord
          ^ref_clv/begining of contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 3
    read.reference_end = 5
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 2))

    contig = MagicMock()
    contig.infer_query_length.return_value = 6  # including hardclip
    contig.cigartuples = (
        # note from right to left
        (S.BAM_CHARD_CLIP, 2),
        (S.BAM_CMATCH, 4),
    )

    ctg_offset = 0
    tail_len = 3
    assert bridge.do_rev_ctg_lt_bdg(read, contig) == ('+', ctg_offset, tail_len)
