from collections import defaultdict
from unittest.mock import MagicMock

from kleat.evidence import bridge
from kleat.misc.apautils import gen_clv_key_tuple
import kleat.misc.settings as S


def test_do_do_fwd_ctg_lt_bdg_with_left_hard_clipping():
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

    ctg_offset = 2              # 4 -2
    tail_len = 3
    assert bridge.do_fwd_ctg_lt_bdg(read, contig) == ('-', ctg_offset, tail_len)


def test_do_fwd_ctg_lt_bdg_with_right_hard_clipping():
    """
    right hardclipping won't have an effect in such case

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

    ctg_offset = 3
    tail_len = 2
    assert bridge.do_fwd_ctg_lt_bdg(read, contig) == ('-', ctg_offset, tail_len)


def test_do_fwd_ctg_rt_bdg_with_left_hardclipping():
    """
              AA
           CCG┘      <-right-tail read
       \\\XCCGXX     <-contig
       0123456789    <-contig coord
          |  ^ctg_offset
    ...XXXCCGXX...   <-reference genome
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

    ctg_offset = 3
    tail_len = 2
    assert bridge.do_fwd_ctg_rt_bdg(read, contig) == ('+', ctg_offset, tail_len)
