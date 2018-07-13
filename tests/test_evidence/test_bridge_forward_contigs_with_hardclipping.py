from unittest.mock import MagicMock

from kleat.evidence import bridge
import kleat.misc.settings as S


###################################################
# test different situations for do_fwd_ctg_lt_bdg #
###################################################

def test_do_fwd_ctg_lt_bdg_with_left_hard_clipping():
    """
      TTT
        └ACG      <-left-tail read
     \\XXACGX     <-contig
     01234567     <-contig coord
         ^ctg_offset
    ...XXACGX...  <-reference genome
       456789     <-genome coord
       | ^ref_clv
       ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 4
    read.reference_end = 7
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CHARD_CLIP, 2), (S.BAM_CMATCH, 6))
    contig.infer_query_length.return_value = 8  # including hardclip

    ctg_offset = 2              # 4 -2
    tail_len = 3
    assert bridge.do_fwd_ctg_lt_bdg(read, contig) == ('-', ctg_offset, tail_len)


def test_do_fwd_ctg_lt_bdg_with_left_hard_clipping_right_after_ctg_clv():
    """
      TTT
        └ACG      <-left-tail read
     \\\\\CGX     <-contig, the "\\\\\" part would appear in another position in genome, although the read is well aligned to it
     012345678    <-contig coord
         ^ctg_offset(won't be captured by this read)
       ...CGX...  <-reference genome
         78901    <-genome coord
         ^starting the contig2genome alignment

    ctg_offset would be 4 (ctg_clv) - 5 (hardclip) = -1 < 0, so this read
    won't capture the genome offset of the clv, but its mate potentially will.
    """
    read = MagicMock()
    read.reference_start = 4
    read.reference_end = 7
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CHARD_CLIP, 5), (S.BAM_CMATCH, 3))
    contig.infer_query_length.return_value = 8  # including hardclip

    assert bridge.do_fwd_ctg_lt_bdg(read, contig) is None


def test_do_fwd_ctg_lt_bdg_with_left_hard_clipping_right_before_ctg_clv():
    """
      TTT
        └ACG      <-left-tail read
     \\\\ACGX     <-contig
     012345678    <-contig coord
         ^ctg_offset
      ...ACGX...  <-reference genome
         78901    <-genome coord
         ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 4
    read.reference_end = 7
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CHARD_CLIP, 4), (S.BAM_CMATCH, 4))
    contig.infer_query_length.return_value = 8  # including hardclip

    ctg_offset = 0              # due to hardclipping
    tail_len = 3
    assert bridge.do_fwd_ctg_lt_bdg(read, contig) == ('-', ctg_offset, tail_len)


def test_do_fwd_ctg_lt_bdg_with_left_hard_clipping_1bp_before_ctg_clv():
    """
      TTT
        └ACG      <-left-tail read
     \\\GACGX     <-contig
     01234567     <-contig coord
         ^ctg_offset
     ...GACGX...  <-reference genome
        678901    <-genome coord
         ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 4
    read.reference_end = 7
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 3), (S.BAM_CMATCH, 3))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 5))
    contig.infer_query_length.return_value = 8  # including hardclip

    ctg_offset = 1              # due to hardclipping
    tail_len = 3
    assert bridge.do_fwd_ctg_lt_bdg(read, contig) == ('-', ctg_offset, tail_len)


def test_do_fwd_ctg_lt_bdg_with_right_hard_clipping():
    """
    such right hardclipping (not passing the ctg_clv) won't have an effect in such case

       TT
        └AC        <-left-tail read
      XXXACGXX//   <-contig
      0123456789   <-contig coord
         ^ctg_offset
    ..XXXACGXX...  <-reference genome
      34567890     <-genome coord
      |  ^ref_clv
      ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 3
    read.reference_end = 5
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 2))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CMATCH, 8), (S.BAM_CHARD_CLIP, 2))
    contig.infer_query_length.return_value = 10  # including hardclip

    ctg_offset = 3
    tail_len = 2
    assert bridge.do_fwd_ctg_lt_bdg(read, contig) == ('-', ctg_offset, tail_len)


def test_do_fwd_ctg_lt_bdg_with_right_hard_clipping_passing_ctg_clv():
    """
       TT
        └AC        <-left-tail read
      XX////       <-contig
      012345       <-contig coord
         ^ctg_offset
    ..XX...        <-reference genome
      34           <-genome coord
      |  ^ref_clv
      ^starting the contig2genome alignment

    ctg_offset would be 3 (ctg_clv) - 5 (hardclip) = -1 < 0, so this read
    """
    read = MagicMock()
    read.reference_start = 3
    read.reference_end = 5
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 2))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CMATCH, 2), (S.BAM_CHARD_CLIP, 4))
    contig.infer_query_length.return_value = 6  # including hardclip

    assert bridge.do_fwd_ctg_lt_bdg(read, contig) is None


def test_do_fwd_ctg_lt_bdg_with_right_hard_clipping_right_on_ctg_clv_edgecase():
    """
       TT
        └AC        <-left-tail read
       XX///       <-contig
       012345       <-contig coord
         ^ctg_offset
    ...XX..        <-reference genome
       34           <-genome coord
       |  ^ref_clv
       ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 2
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 2))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CMATCH, 2), (S.BAM_CHARD_CLIP, 3))
    contig.infer_query_length.return_value = 5  # including hardclip

    assert bridge.do_fwd_ctg_lt_bdg(read, contig) is None


def test_do_fwd_ctg_lt_bdg_with_right_hard_clipping_right_after_ctg_clv_edgecase():
    """
       TT
        └AC        <-left-tail read
      XXXA//       <-contig
      012345       <-contig coord
         ^ctg_offset
    ..XX...        <-reference genome
      34           <-genome coord
      |  ^ref_clv
      ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 3
    read.reference_end = 5
    read.cigartuples = ((S.BAM_CSOFT_CLIP, 2), (S.BAM_CMATCH, 2))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CHARD_CLIP, 2))
    contig.infer_query_length.return_value = 6  # including hardclip

    ctg_offset = 3
    tail_len = 2
    assert bridge.do_fwd_ctg_lt_bdg(read, contig) == ('-', ctg_offset, tail_len)


###################################################
# test different situations for do_fwd_ctg_rt_bdg #
###################################################


def test_do_fwd_ctg_rt_bdg_with_left_hardclipping():
    """
              AA
           CCG┘      <-right-tail read
       \\\XCCGXX     <-contig
       0123456789    <-contig coord
          |  ^ctg_offset
    ...XXXXCCGXX...   <-reference genome
          4567890    <-genome coord
          |  ^ref_clv
          ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 4
    read.reference_end = 7
    read.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 7))
    contig.infer_query_length.return_value = 9

    ctg_offset = 3
    tail_len = 2
    assert bridge.do_fwd_ctg_rt_bdg(read, contig) == ('+', ctg_offset, tail_len)


def test_do_fwd_ctg_rt_bdg_with_left_hardclipping_passing_ctg_clv():
    """
            AA
         CCG┘       <-right-tail read
        \\\\\X      <-contig
        0123456     <-contig coord
           ^ctg_offset
          ...X...   <-reference genome
           45678    <-genome coord
    ref_clv^ ^starting the contig2genome alignment

    ref_clv won't be captured by this bridge read
    """
    read = MagicMock()
    read.reference_start = 1
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CHARD_CLIP, 5), (S.BAM_CMATCH, 1))
    contig.infer_query_length.return_value = 6

    assert bridge.do_fwd_ctg_rt_bdg(read, contig) is None


def test_do_fwd_ctg_rt_bdg_with_left_hardclipping_right_on_ctg_clv():
    """
            AA
         CCG┘       <-right-tail read
        \\\\XX      <-contig
        0123456     <-contig coord
           ^ctg_offset
         ...XX...   <-reference genome
           45678    <-genome coord
    ref_clv^ ^starting the contig2genome alignment

    ref_clv won't be captured by this bridge read
    """
    read = MagicMock()
    read.reference_start = 1
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CMATCH, 3), (S.BAM_CSOFT_CLIP, 2))

    contig = MagicMock()
    contig.cigartuples = ((S.BAM_CHARD_CLIP, 4), (S.BAM_CMATCH, 2))
    contig.infer_query_length.return_value = 6

    assert bridge.do_fwd_ctg_rt_bdg(read, contig) is None


def test_do_fwd_ctg_rt_bdg_with_right_hardclipping():
    """
    right hardclipping won't have an effect in such case
           AAAAA
         CG┘       <-right-tail read
       XCCGXX//    <-contig
       012345678   <-contig coord
       |  ^ctg_offset
    ...XCCGXX...   <-reference genome
       345678901   <-genome coord
       |  ^ref_clv
       ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 2
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CMATCH, 2), (S.BAM_CSOFT_CLIP, 5))

    contig = MagicMock()
    contig.infer_query_length.return_value = 8
    contig.cigartuples = (
        (S.BAM_CMATCH, 6),
        (S.BAM_CHARD_CLIP, 2),
    )

    ctg_offset = 3
    tail_len = 5
    assert bridge.do_fwd_ctg_rt_bdg(read, contig) == ('+', ctg_offset, tail_len)


def test_do_fwd_ctg_rt_bdg_with_right_hardclipping_right_before_ctg_clv():
    """
    right hardclipping won't have an effect in such case

           AAAAA
         CG┘       <-right-tail read
       XCCG//      <-contig
       0123456     <-contig coord
       |  ^ctg_offset
    ...XCCG...     <-reference genome
       34567       <-genome coord
       |  ^ref_clv
       ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 2
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CMATCH, 2), (S.BAM_CSOFT_CLIP, 5))

    contig = MagicMock()
    contig.infer_query_length.return_value = 6
    contig.cigartuples = (
        (S.BAM_CMATCH, 4),
        (S.BAM_CHARD_CLIP, 2),
    )

    ctg_offset = 3
    tail_len = 5
    assert bridge.do_fwd_ctg_rt_bdg(read, contig) == ('+', ctg_offset, tail_len)


def test_do_fwd_ctg_rt_bdg_with_right_hardclipping_passing_ctg_clv():
    """
           AAAAA
         CG┘       <-right-tail read
       XC///       <-contig
       012345678   <-contig coord
       |  ^ctg_offset
    ...XC...       <-reference genome
       345678901   <-genome coord
       |  ^ref_clv (won't be captured by this bridge read)
       ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 2
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CMATCH, 2), (S.BAM_CSOFT_CLIP, 5))

    contig = MagicMock()
    contig.infer_query_length.return_value = 5
    contig.cigartuples = (
        (S.BAM_CMATCH, 2),
        (S.BAM_CHARD_CLIP, 3),
    )

    assert bridge.do_fwd_ctg_rt_bdg(read, contig) is None


def test_do_fwd_ctg_rt_bdg_with_right_hardclipping_right_on_ctg_clv():
    """
           AAAAA
         CG┘       <-right-tail read
       XCC///      <-contig
       0123456     <-contig coord
       |  ^ctg_offset
    ...XCC...      <-reference genome
       3456789     <-genome coord
       |  ^ref_clv
       ^starting the contig2genome alignment
    """
    read = MagicMock()
    read.reference_start = 2
    read.reference_end = 4
    read.cigartuples = ((S.BAM_CMATCH, 2), (S.BAM_CSOFT_CLIP, 5))

    contig = MagicMock()
    contig.infer_query_length.return_value = 6
    contig.cigartuples = (
        (S.BAM_CMATCH, 3),
        (S.BAM_CHARD_CLIP, 3),
    )

    assert bridge.do_fwd_ctg_rt_bdg(read, contig) is None
