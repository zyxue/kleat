from unittest.mock import MagicMock, patch

import kleat.misc.settings as S
from kleat.hexamer.xseq_plus import init_ctg_end, init_ref_end


"""
cc: ctg_clv; icb: init_clv_beg
rc: ref_clv; irb: init_ref_end
"""


def test_init_ends():
    """
             AA
           GT┘       <-bridge read
       GACGGTTGC     <-bridge contig
       0123456789    <-contig coord
          cc^   ^ice
    ...GACGGTTGC...  <-genome
       5678901234    <-genome coord
            1   |
          rc^   ^ire
    """
    ref_clv = 10
    cigartuples = ((S.BAM_CMATCH, 9),)
    ctg_clv = 5
    ctg_seq = 'GACGGTTGC'

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 14
    assert init_ctg_end(ctg_seq) == 9


def test_init_ends_with_skip_after_clv():
    """
             AA
           GT┘       <-bridge read
       GACGGT-GC     <-bridge contig
       012345 678    <-contig coord
          cc^   ^ice
    ...GACGGTTGC...  <-genome
       5678901234    <-genome coord
            1   |
          rc^   ^ire
    """
    ref_clv = 10
    cigartuples = ((S.BAM_CMATCH, 6), (S.BAM_CREF_SKIP, 1), (S.BAM_CMATCH, 2))
    ctg_clv = 5
    ctg_seq = 'GACGGTGC'

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 14
    assert init_ctg_end(ctg_seq) == 8


def test_init_ends_with_skip_before_clv():
    """
             AA
           GT┘       <-bridge read
       G--AGTTGC     <-bridge contig
       0  1234567    <-contig coord
          cc^   ^ice
    ...GACAGTTGC...  <-genome
       5678901234    <-genome coord
            1   |
          rc^   ^ire
    """
    ref_clv = 10
    cigartuples = ((S.BAM_CMATCH, 1), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 6))
    ctg_clv = 3
    ctg_seq = 'GAGTTGC'

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 14
    assert init_ctg_end(ctg_seq) == 7


def test_init_ends_with_skip_both_before_and_after_clv():
    """
             AA
           GT┘       <-bridge read
       G--AGT-GC     <-bridge contig
       0  123 456    <-contig coord
          cc^   ^ice
    ...GACAGTTGC...  <-genome
       5678901234    <-genome coord
            1   |
          rc^   ^ire
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

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 14
    assert init_ctg_end(ctg_seq) == 6


def test_init_ends_with_skip_both_before_and_after_ctg_clv_and_a_mismatch():
    """
             AA
           GT┘        <-bridge read
       G--AGT-GC      <-bridge contig
       0  x23 456     <-contig coord
          cc^   ^ice
    ...GACAGTTGC...   <-genome
       5678901234     <-genome coord
            1   |
          rc^   ^ire
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

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 14
    assert init_ctg_end(ctg_seq) == 6


def test_init_ends_with_multiple_skips_before_clv():
    """
               AA
             TA┘       <-bridge read
       G-C--CTAGC      <-bridge contig
       0 1  234567     <-contig coord
        ||| x ^cc^ice
    ...GACTGGTAGC...   <-genome
       56789012345     <-genome coord
            1 |  |
            rc^  ^ire
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

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 15
    assert init_ctg_end(ctg_seq) == 7


def test_init_ends_with_deletion():
    """
               AA
             CG┘      <-bridge read
       GAC__TCGTC     <-bridge contig
       012  345678    <-contig coord
          | ||x
          | |^cc ^ice
    ...GACGGTCCTC...  <-genome
       56789012345    <-genome coord
            1|   |
             ^rc ^rce
    """

    ref_clv = 11
    cigartuples = (
        (S.BAM_CMATCH, 3),
        (S.BAM_CDEL, 2),
        (S.BAM_CMATCH, 5),
    )
    ctg_clv = 4
    ctg_seq = 'GACTCGTC'

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 15
    assert init_ctg_end(ctg_seq) == 8


def test_init_ends_with_insertion():
    """
         AG   AA      <-inserted bases
         ┬  GT┘       <-bread read
       GA CGGTCGC     <-bridge contig
       01 45678901    <-contig coord
        x    |  1|
        x  cc^   ^ice
    ...GT CGGTCGC...  <-genome
       56 78901234    <-genome coord
             |   |
           rc^   ^ire
    """

    ref_clv = 10
    cigartuples = (
        (S.BAM_CMATCH, 2),
        (S.BAM_CINS, 2),
        (S.BAM_CMATCH, 7)
    )
    ctg_clv = 7
    ctg_seq = 'GAAGCGGTCGC'

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 14
    assert init_ctg_end(ctg_seq) == 11


def test_bridge_init_ends_with_sofclip_before_clv():
    """
             AAA                <-polyA softclip
          GTT┘                  <-bridge read
    CCC   |||                   <-non-polyA softclip
    01└GA-GGTTGCAGA             <-suffix contig
       |  |||  ////             <-hardclip mask
       34 5678901234            <-contig coord
          cc^  |   ^ice
    ...ACGGTTGC...              <-genome
       7890123456789             <-genome coord
          1 |      |
          rc^      ^init_ref_idx
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

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 19
    assert init_ctg_end(ctg_seq) == 14


def test_bridge_init_ends_with_softclip_after_clv():
    """
           AA       <-polyA clip
         TC┘|       <-bridge read
         ||  CC     <-non-polyA softclip
    \\\ATTCGT┘|     <-bridge contig (hardcipped, could be chimeric https://www.biostars.org/p/109333/)
       012345678    <-contig coord
        cc^    ^ice
    ...ATTCGTXX...  <-genome
       567890123    <-genome coord
          | 1  |
        rc^    ^ire
    """
    ref_clv = 8
    cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 6), (S.BAM_CSOFT_CLIP, 2))
    ctg_clv = 3
    ctg_seq = 'ATTCGTCC'

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 13
    assert init_ctg_end(ctg_seq) == 8


def test_bridge_init_ends_with_hardclip_before_clv():
    """
           AA
         TC┘|      <-bridge read
    \\\ATTCGT      <-bridge contig (hardcipped, could be chimeric https://www.biostars.org/p/109333/)
       0123456     <-contig coord
        cc^  ^ice
    ...ATTCGXX...  <-genome
       5678901     <-genome coord
          | 1|
        rc^  ^ire
    """

    ref_clv = 8
    cigartuples = ((S.BAM_CHARD_CLIP, 3), (S.BAM_CMATCH, 6))
    ctg_clv = 3
    ctg_seq = 'ATTCGT'

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 11
    assert init_ctg_end(ctg_seq) == 6


def test_bridge_init_ends_with_hardclip_after_clv():
    """
             AAA
          GTT┘                  <-bridge read
       A-GGTTGCAGA              <-suffix contig
       | |  |  ///              <-hardclip mask
       0 1234567890             <-contig coord
          cc^     ^ice
    ...ACGGTTGC...              <-genome
       789012345678             <-genome coord
          1 |     |
          rc^     ^init_ref_idx
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

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 18
    assert init_ctg_end(ctg_seq) == 10


def test_link_init_ends():
    """
          TT...AA
       ATCGAC    <-link contig
       0123456   <-contig coord
            ^ctg_clv
    ...7890123... <-genome coord
          1 ^ref_clv
    """
    ref_clv = 12
    cigartuples = ((S.BAM_CMATCH, 6),)
    ctg_clv = 5
    ctg_seq = 'ATCGAC'

    assert init_ref_end(ref_clv, cigartuples, ctg_clv, ctg_seq) == 13
    assert init_ctg_end(ctg_seq) == 6
