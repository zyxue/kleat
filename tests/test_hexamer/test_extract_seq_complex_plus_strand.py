import unittest
from unittest.mock import MagicMock, call

import kleat.misc.settings as S
from kleat.hexamer.search import extract_seq


def test_extract_seq_with_skipped_region_for_plus_strand_clv():
    """
               AAA             <-tail of suffix contig
       ACGG--GC┘||             <-suffix contig with skip
       0123  456789            <-contig coord
              |  1
       ctg_clv^  ^init_ctg_idx <-contig coord
    ...ACGGTTGCGGT...          <-genome
       789012345678            <-genome coord
          1   |  |
       ref_clv^  ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    ctg.query_sequence = 'ACGGGCAAA'
    ctg.cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CREF_SKIP, 2), (S.BAM_CMATCH, 2), (S.BAM_CSOFT_CLIP, 3))

    ref_fa = MagicMock()
    ref_fa.fetch = MagicMock(return_value='TT')
    kw = dict(contig=ctg, strand='+', ref_clv=14, ref_fa=ref_fa, ctg_clv=5)
    assert extract_seq(**kw) == 'ACGGTTGC'
    assert extract_seq(**kw, window=3) == 'TGC'
    ref_fa.fetch.assert_called_with('chr1', 11, 13)


def test_extract_seq_with_1_base_insertion_for_plus_strand_clv():
    """
           T
           ┬  AAA              <-tail of suffix contig
       ACGG GC┘||              <-suffix contig with skip
       0123 56789             <-contig coord
          | ||  |
      ctg_clv^  ^init_ctg_idx <-contig coord
    ...ACGG GCXXX...           <-genome
       7890 1234567            <-genome coord
          1  |  |
      ref_clv^  ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    ctg.query_sequence = 'ACGGTGCAAA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 4), 
        (S.BAM_CINS, 5),
        (S.BAM_CMATCH, 2),
        (S.BAM_CSOFT_CLIP, 3),
    )

    ref_fa = MagicMock()
    kw = dict(contig=ctg, strand='+', ref_clv=12, ref_fa=ref_fa, ctg_clv=6)
    assert extract_seq(**kw) == 'ACGGTGC'
    assert extract_seq(**kw, window=4) == 'GTGC'


def test_extract_seq_with_5_base_inserted_region_for_plus_strand_clv():
    """
         AATCC
           ┬   AA              <-tail of suffix contig
       ACGG GCG┘|              <-suffix contig with skip
       0123 9012            <-contig coord
          | |1| |
       ctg_clv^ ^init_ctg_idx <-contig coord
    ...ACGG GCGXXX...           <-genome
       7890 1234567            <-genome coord
         1    | |
       ref_clv^  ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    ctg.query_sequence = 'ACGGAATCCGCGAA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 4), 
        (S.BAM_CINS, 5),
        (S.BAM_CMATCH, 3),
        (S.BAM_CSOFT_CLIP, 2),
    )

    ref_fa = MagicMock()
    kw = dict(contig=ctg, strand='+', ref_clv=13, ref_fa=ref_fa, ctg_clv=10)
    assert extract_seq(**kw) == 'ACGGAATCCGCG'
    assert extract_seq(**kw, window=5) == 'CCGCG'


def test_extract_seq_with_deleted_region_for_plus_strand_clv():
    """
               AAA             <-tail of suffix contig
       ACGG__GC┘||             <-suffix contig with skip
       0123  456789            <-contig coord
              |  1
       ctg_clv^  ^init_ctg_idx <-contig coord
    ...ACGGTTGCGGT...          <-genome
       789012345678            <-genome coord
          1   |  |
       ref_clv^  ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr1'
    ctg.query_sequence = 'ACGGGCAAA'  # len 9
    ctg.cigartuples = ((S.BAM_CMATCH, 4), (S.BAM_CDEL, 2), (S.BAM_CMATCH, 2), (S.BAM_CSOFT_CLIP, 3))

    ref_fa = MagicMock()
    kw = dict(contig=ctg, strand='+', ref_clv=14, ref_fa=ref_fa, ctg_clv=5)
    assert extract_seq(**kw) == 'ACGGGC'
    assert extract_seq(**kw, window=3) == 'GGC'


def test_extract_seq_with_two_skipped_region_for_plus_strand_clv():
    """
               AAA             <-tail of suffix contig
       A-TT--GC┘||             <-suffix contig with skip
       0 12  345678            <-contig coord
              |
       ctg_clv^  ^init_ctg_idx <-contig coord
    ...ACTTAAGCGGT...          <-genome
       789012345678            <-genome coord
          1   |  |
       ref_clv^  ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr3'
    ctg.query_sequence = 'ATTGCAAA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 2),
        (S.BAM_CSOFT_CLIP, 3),
    )

    ref_fa = MagicMock()
    ref_fa.fetch.side_effect = ['AA', 'C']
    kw = dict(contig=ctg, strand='+', ref_clv=14, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'ACTTAAGC'
    assert ref_fa.fetch.call_count == 2
    ref_fa.fetch.assert_has_calls([call('chr3', 11, 13), call('chr3', 8, 9)])

    # use a new mock, couldn't make ref_fa.reset_mock() work
    ref_fa = MagicMock()
    ref_fa.fetch.return_value = 'AA'
    kw.update(ref_fa=ref_fa)
    assert extract_seq(**kw, window=3) == 'AGC'
    assert ref_fa.fetch.call_count == 1


def test_extract_seq_with_skipped_region_and_insertions_mismatches_for_plus_strand_clv():
    """
        G
        ┬       AAA             <-tail of suffix contig
       A TA--GCG┘||             <-suffix contig with skip
       0 23  456789            <-contig coord
       | |x  | |
       ctg_clv ^  ^init_ctg_idx <-contig coord
    ...A TTCCGCGXXX...          <-genome
       7 8901234567            <-genome coord
           1   |  |
        ref_clv^  ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr3'
    ctg.query_sequence = 'AGTAGCGAAA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CINS, 1),
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CSOFT_CLIP, 3),
    )

    ref_fa = MagicMock()
    ref_fa.fetch.return_value = 'CC'
    kw = dict(contig=ctg, strand='+', ref_clv=14, ref_fa=ref_fa, ctg_clv=6)
    assert extract_seq(**kw) == 'AGTACCGCG'
    ref_fa.fetch.assert_called_once_with('chr3', 10, 12)
    assert extract_seq(**kw, window=8) == 'GTACCGCG'


def test_extract_seq_with_skipped_and_deleted_regions_for_plus_strand_clv():
    """
               AAA             <-tail of suffix contig
       A_TT--GC┘||             <-suffix contig with skip
       0 12  345678            <-contig coord
              |
       ctg_clv^  ^init_ctg_idx <-contig coord
    ...ACTTAAGCGGT...          <-genome
       789012345678            <-genome coord
          1   |  |
       ref_clv^  ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr3'
    ctg.query_sequence = 'ATTGCAAA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CDEL, 1),
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 2),
        (S.BAM_CSOFT_CLIP, 3),
    )

    ref_fa = MagicMock()
    ref_fa.fetch.return_value = 'AA'
    kw = dict(contig=ctg, strand='+', ref_clv=14, ref_fa=ref_fa, ctg_clv=4)
    assert extract_seq(**kw) == 'ATTAAGC'
    ref_fa.fetch.assert_called_once_with('chr3', 11, 13)
    assert extract_seq(**kw, window=5) == 'TAAGC'


def test_extract_seq_with_three_skipped_region_and_mismatches_for_plus_strand_clv():
    """
                     AA             <-tail of suffix contig
       A---CC-GTA--GC┘|             <-suffix contig with skip
       0|||12|345||678              <-contig coord
        |||x | x || |
        |||x ctg_clv^ ^init_ctg_idx <-contig coord
    ...ACTGTCAGAATTGC...            <-genome
       78901234567890               <-genome coord
          1         | |
             ref_clv^ ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr3'
    ctg.query_sequence = 'ACCGTAGCAA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 3),
        (S.BAM_CMATCH, 2),
        (S.BAM_CREF_SKIP, 1),
        (S.BAM_CMATCH, 3),
        (S.BAM_CREF_SKIP, 2),
        (S.BAM_CMATCH, 2),
        (S.BAM_CSOFT_CLIP, 2),
    )

    ref_fa = MagicMock()
    ref_fa.fetch.side_effect = ['TT', 'A', 'CTG']
    kw = dict(contig=ctg, strand='+', ref_clv=20, ref_fa=ref_fa, ctg_clv=7)
    assert extract_seq(**kw) == 'ACTGCCAGTATTGC'
    assert ref_fa.fetch.call_count == 3
    ref_fa.fetch.assert_has_calls([call('chr3', 17, 19), call('chr3', 13, 14), call('chr3', 8, 11)])

    ref_fa = MagicMock()
    ref_fa.fetch.side_effect = ['TT', 'A']
    kw.update(ref_fa=ref_fa)
    assert extract_seq(**kw, window=9) == 'CAGTATTGC'
    assert ref_fa.fetch.call_count == 2


def test_extract_seq_with_indel_and_skipped_regions_and_mismatches_for_plus_strand_clv():
    """
             TC
             ┬       AA             <-tail of suffix contig
       A---CC GTA__GC┘|             <-suffix contig with skip
       0   12 567  8901             <-contig coord
           x   x    |1
           x ctg_clv^ ^init_ctg_idx <-contig coord
    ...ACTGTC GAATTGC...            <-genome
       789012 345678901             <-genome coord
          1         | |
             ref_clv^ ^init_ref_idx
    """
    ctg = MagicMock()
    ctg.reference_name = 'chr3'
    ctg.query_sequence = 'ACCTCGTAGCAA'
    ctg.cigartuples = (
        (S.BAM_CMATCH, 1),
        (S.BAM_CREF_SKIP, 3),
        (S.BAM_CMATCH, 2),
        (S.BAM_CINS, 2),
        (S.BAM_CMATCH, 3),
        (S.BAM_CDEL, 2),
        (S.BAM_CMATCH, 2),
        (S.BAM_CSOFT_CLIP, 2),
    )

    ref_fa = MagicMock()
    ref_fa.fetch.return_value = 'CTG'
    kw = dict(contig=ctg, strand='+', ref_clv=19, ref_fa=ref_fa, ctg_clv=9)
    assert extract_seq(**kw) == 'ACTGCCTCGTAGC'
    ref_fa.fetch.assert_called_once_with('chr3', 8, 11)

    assert extract_seq(**kw, window=10) == 'GCCTCGTAGC'
