from unittest.mock import MagicMock, patch

import kleat.misc.settings as S
from kleat.hexamer.xseq_minus import init_ctg_beg, init_ref_beg


"""
cc: ctg_clv; icb: init_clv_beg
rc: ref_clv; irb: init_ref_end
"""


def test_init_begs():
    """
         TT
          └GT        <-bridge read
       GACGGTTGC     <-bridge contig
       0123456789    <-contig coord
    irb^    ^cc
    ...GACGGTTGC...  <-genome
       5678901234    <-genome coord
       |    1
    irc^    ^rc
    """
    ref_clv = 10
    cigartuples = ((S.BAM_CMATCH, 9),)
    ctg_clv = 5
    ctg_seq = 'GACGGTTGC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 5
    assert init_ctg_beg(ctg_seq) == 0


def test_init_begs_with_skip_after_clv():
    """
         TT
          └GT       <-bridge read
       GACGGT-GC     <-bridge contig
       012345 678    <-contig coord
    irc^    ^cc
    ...GACGGTTGC...  <-genome
       5678901234    <-genome coord
       |    1
    irb^    ^rc
    """
    ref_clv = 10
    cigartuples = ((S.BAM_CMATCH, 6), (S.BAM_CREF_SKIP, 1), (S.BAM_CMATCH, 2))
    ctg_clv = 5
    ctg_seq = 'GACGGTGC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 5
    assert init_ctg_beg(ctg_seq) == 0


def test_init_begs_with_skip_before_clv():
    """
         TT
          └GT       <-bridge read
       G--AGTTGC     <-bridge contig
       0  1234567    <-contig coord
    icb^    ^cc
    ...GACAGTTGC...  <-genome
       5678901234    <-genome coord
       |    1
    icb^    ^rc
    """
    ref_clv = 10
    cigartuples = ((S.BAM_CMATCH, 1), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 6))
    ctg_clv = 3
    ctg_seq = 'GAGTTGC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 5
    assert init_ctg_beg(ctg_seq) == 0


def test_init_begs_with_skip_both_before_and_after_clv():
    """
         TT
          └GT       <-bridge read
       G--AGT-GC     <-bridge contig
       0  123 456    <-contig coord
    icb^    ^cc
    ...GACAGTTGC...  <-genome
       5678901234    <-genome coord
       |    1
    ifb^    ^rc
    """
    ref_clv = 10
    cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 2),
    )
    ctg_clv = 3
    ctg_seq = 'GAGTGC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 5
    assert init_ctg_beg(ctg_seq) == 0


def test_init_begs_with_skip_both_before_and_after_ctg_clv_and_a_mismatch():
    """
         TT
          └GT        <-bridge read
       G--TGT-GC      <-bridge contig
       0  123 456     <-contig coord
    irb^  x ^cc
    ...GACAGTTGC...   <-genome
       5678901234     <-genome coord
       |    1
    ifb^    ^rc
    """

    ref_clv = 10
    cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 2),
    )
    ctg_clv = 3
    ctg_seq = 'GAGTGC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 5
    assert init_ctg_beg(ctg_seq) == 0


def test_init_begs_with_multiple_skips_before_clv():
    """
           TT
            └TA        <-bridge read
       G-C--CTAGC      <-bridge contig
       0 1  234567     <-contig coord
    icb^    x ^cc
    ...GACTGGTAGC...   <-genome
       56789012345     <-genome coord
       |    1 |  |
    irb^      ^rc
    """

    ref_clv = 12
    cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 5)
    )
    ctg_clv = 4
    ctg_seq = 'GCCTAGC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 5
    assert init_ctg_beg(ctg_seq) == 0


def test_init_begs_with_deletion():
    """
           TT
            └CG       <-bridge read
       GAC__TCGTC     <-bridge contig
       012  345678    <-contig coord
       |  | | |
    icb^  | | ^cc
    ...GACGGTCCTC...  <-genome
       56789012345    <-genome coord
       |    1 |
    irb^      ^rc
    """

    ref_clv = 11
    cigartuples = (
        (S.BAM_CMATCH, 3),
        (S.BAM_CDEL, 2),
        (S.BAM_CMATCH, 5),
    )
    ctg_clv = 4
    ctg_seq = 'GACTCGTC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 5
    assert init_ctg_beg(ctg_seq) == 0


def test_init_begs_with_insertion():
    """
        AG   TT       <-inserted bases
         ┬    └G      <-bread read
       GA CGGTCGC     <-bridge contig
       01 45678901    <-contig coord
       |x      |1
    icb^x      ^cc
    ...GT CGGTCGC...  <-genome
       56 78901234    <-genome coord
       |       |
    irb^       ^rc
    """

    ref_clv = 12
    cigartuples = (
        (S.BAM_CMATCH, 2),
        (S.BAM_CINS, 2),
        (S.BAM_CMATCH, 7)
    )
    ctg_clv = 9
    ctg_seq = 'GAAGCGGTCGC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 5
    assert init_ctg_beg(ctg_seq) == 0


def test_bridge_init_begs_with_sofclip_before_clv():
    """
          TTT                      <-polyA softclip
            └GTT                   <-bridge read
       CCC   |||                   <-non-polyA softclip
       ||└GA-GGTTGCAGA             <-suffix contig
       || |  |||  ////             <-hardclip mask
       01234 5678901234            <-contig coord
    icb^       ^cc|
       ...ACGGTTGC...              <-genome
       4567890123456789             <-genome coord
       |     1 |
    irb^       ^rc
    """
    ref_clv = 12
    cigartuples = (
        (S.BAM_CSOFT_CLIP, 3),
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 5),
        (S.BAM_CHARD_CLIP, 4)
    )
    ctg_clv = 7
    ctg_seq = 'CCCGAGGTTGCAGA'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 4
    assert init_ctg_beg(ctg_seq) == 0


def test_bridge_init_begs_with_softclip_after_clv():
    """
          TT           <-polyA clip
           └TC         <-bridge read
            ||  CC     <-non-polyA softclip
       \\\ATTCGT┘|     <-bridge contig (hardcipped, could be chimeric https://www.biostars.org/p/109333/)
       01234567890     <-contig coord
    icb^     ^cc
       ...ATTCGTXX...  <-genome
       234567890123    <-genome coord
       |     | 1  |
    irb^     ^rc
    """
    ref_clv = 8
    cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 6), (S.BAM_CSOFT_CLIP, 2))
    ctg_clv = 6
    ctg_seq = 'ATTCGTCC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 2
    assert init_ctg_beg(ctg_seq) == 0


def test_bridge_init_begs_with_hardclip_before_clv():
    """
          TT
           └TC        <-bridge read
       \\\ATTCGT      <-bridge contig (hardcipped, could be chimeric https://www.biostars.org/p/109333/)
       0123456789     <-contig coord
    icb^     ^cc
       ...ATTCGXX...  <-genome
       2345678901     <-genome coord
       |     | 1
    irb^     ^rc
    """

    ref_clv = 8
    cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 6))
    ctg_clv = 6
    ctg_seq = 'ATTCGT'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 2
    assert init_ctg_beg(ctg_seq) == 0


def test_bridge_init_begs_with_hardclip_after_clv():
    """
       TTT
         └GTT                   <-bridge read
       A-GGTTGCAGA              <-suffix contig
       | |  |  ///              <-hardclip mask
       0 1234567890             <-contig coord
    icb^    ^cc
    ...ACGGTTGC...              <-genome
       789012345678             <-genome coord
       |  1 |
    irb^    ^rc
    """

    ref_clv = 12
    cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 6),
        (S.BAM_CHARD_CLIP, 3)
    )
    ctg_clv = 4
    ctg_seq = 'AGGTTGCAGA'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 7
    assert init_ctg_beg(ctg_seq) == 0


def test_link_init_begs():
    """
    TT......AA
          ATCGAC    <-link contig
          0123456   <-contig coord
       icb^ ^ctg_clv
       ...7890123... <-genome coord
       irb^ ^ref_clv
    """
    ref_clv = 9
    cigartuples = ((S.BAM_CMATCH, 6),)
    ctg_clv = 2
    ctg_seq = 'ATCGAC'

    assert init_ref_beg(ref_clv, cigartuples, ctg_clv) == 7
    assert init_ctg_beg(ctg_seq) == 0
