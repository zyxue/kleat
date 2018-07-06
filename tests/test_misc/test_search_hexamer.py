import os
import unittest

import pysam

from kleat.misc.search_hexamer import (
    plus_search,
    minus_search,
    search,

    fetch_seq,
    gen_coords,
    search_reference_genome,
)


class TestSearchHexamer(unittest.TestCase):
    def test_plus_search(self):
        self.assertEqual(plus_search('GGGAATAAAG', 9), ('AATAAA', 1, 3))
        self.assertEqual(plus_search('GGGAATAAA', 9), ('AATAAA', 1, 4))
        self.assertEqual(plus_search('GGGAATAAAGG', 9), ('AATAAA', 1, 2))
        self.assertEqual(plus_search('GGGATTAAAGG', 9), ('ATTAAA', 2, 2))
        self.assertEqual(plus_search('GGGAATAA', 9), None)

        self.assertEqual(plus_search('GAATAAAC', 10), ('AATAAA', 1, 4))
        self.assertEqual(plus_search('GGGGCTAC', 20), ('GGGGCT', 16, 13))
        self.assertEqual(plus_search('GTTTATTC', 6), None)

    def test_plus_search_lowercase(self):
        self.assertEqual(plus_search('GAATaaaC', 10), ('AATAAA', 1, 4))

    def test_plus_search_take_right_most_hexamer(self):
        self.assertEqual(plus_search('CAATAAANAATAAAC', 200), ('AATAAA', 1, 194))

    def test_plus_search_take_right_most_hexamer_with_Ns(self):
        self.assertEqual(plus_search('GCATTAAAAATNAAC', 200), ('ATTAAA', 2, 188))

    def test_plus_search_take_the_strongest_hexamer(self):
        self.assertEqual(plus_search('GCAATAAAATTAAAC', 200), ('AATAAA', 1, 188))

    def test_minus_search(self):
        self.assertEqual(minus_search('ATTTATTCCC', 9), ('AATAAA', 1, 15))
        self.assertEqual(minus_search('ATTTAATCCC', 9), ('ATTAAA', 2, 15))
        self.assertEqual(minus_search('GTTTATTC', 1), ('AATAAA', 1, 7))
        self.assertEqual(minus_search('ATCGTATATTGC', 5), ('AATATA', 7, 14))

    def test_minus_search_lowercase(self):
        self.assertEqual(minus_search('GTttattc', 1), ('AATAAA', 1, 7))

    def test_minus_search_take_left_most_hexamer(self):
        self.assertEqual(minus_search('GTTTATTTTTATTCG', 10), ('AATAAA', 1, 16))

    def test_minus_search_take_left_most_hexamer_with_Ns(self):
        self.assertEqual(minus_search('GTTTATTNTTTATTNNNTGTATTCG', 10), ('AATAAA', 1, 16))

    def test_minus_search_take_the_strongest_hexamer(self):
        self.assertEqual(minus_search('GTTTAATNTTTATTNNNTGTATTCG', 20), ('AATAAA', 1, 33))

    def test_minus_search_take_the_strongest_hexamer_in_lower_case(self):
        self.assertEqual(minus_search('gtttaatntttattnnntgtattcg', 20), ('AATAAA', 1, 33))


class TestSearch(unittest.TestCase):
    def test_plus_strand(self):
        """
         CaataaaGT
        0123456789 <-genome coord
          |      |
          PAS    clv
        """
        seq = 'CaataaaGT'
        clv = 9
        self.assertEqual(search('+', clv, seq, 50), ('AATAAA', 1, 2))

    def test_minus_strand(self):
        """
         GGTTTATT
        0123456789 <-genome coord
         |      |
         clv    PAS
        """
        seq = 'GGTTTATT'
        clv = 1
        self.assertEqual(search('-', clv, seq, 50), ('AATAAA', 1, 8))


# below are about search PAS hexamer on the reference genome

# TODO: if fa doesn't exist, download, or mock it
REF_FA = os.path.join(
    os.path.dirname(__file__),
    './Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa')


class TestFetchSeq(unittest.TestCase):
    def setUp(self):
        self.refseq = pysam.FastaFile(REF_FA)

    def test_fetch_seq_bound_negative_coordinate_by_0(self):
        chrm = 'MT'
        clv = 3
        beg, end = gen_coords(clv, '+', window=6)
        calc_beg = -2            # clv - 6 + 1
        calc_end = 4             # clv + 1
        self.assertEqual((beg, end), (calc_beg, calc_end))
        expected = 'GATC'         # the first 4 bases of hg19 in chr MT
        self.assertEqual(fetch_seq(self.refseq, chrm, beg, end), expected)
        # assert beginning is bound by 0
        self.assertEqual(fetch_seq(self.refseq, chrm, 0,   end), expected)

    def test_fetch_seq_bound_over_large_coordinate_by_chr_length(self):
        # better have another test case for max chr length
        pass

    def test_fetch_seq_based_on_KLEAT_clv_plus_strand_chr12_DRAM1(self):
        """based on hg19"""
        chrm = '12'
        # this is a coord reported by KLEAT, it's 1-based, so -1 to make it
        # 0-based, and it points to the end of 3'UTR, where there is a hexamer
        # right to the upstream of it, see included screenshot in the repo
        clv = 102316878 - 1
        # 'a' is the last base on 3'UTR
        self.assertEqual(fetch_seq(self.refseq, chrm, clv, clv + 1), 'a')
        beg, end = gen_coords(clv, '+', window=6)
        self.assertEqual((beg, end), (clv - 6 + 1, clv + 1))
        self.assertEqual(fetch_seq(self.refseq, chrm, beg, end), 'aataaa')

    def test_fetch_seq_based_on_KLEAT_clv_minus_strand_chr9_GNAQ(self):
        chrm = '9'
        # reported by KLEAT, converted to 0-based, corresponding to a T, the
        # last base of its 3'UTR
        clv = 80335189 - 1
        self.assertEqual(fetch_seq(self.refseq, chrm, clv, clv + 1), 'T')

        beg, end = gen_coords(clv, '-', window=6)
        self.assertEqual((beg, end), (clv, clv + 6))
        self.assertEqual(fetch_seq(self.refseq, chrm, beg, end), 'TTTTAT')

    def test_fetch_seq_based_on_KLEAT_clv_minus_strand_chr2_NFE2L2(self):
        chrm = '2'
        # reported by KLEAT, converted to 0-based, corresponding to a T, the
        # last base of its 3'UTR
        clv = 178095802 - 1

        self.assertEqual(fetch_seq(self.refseq, chrm, clv, clv + 1), 'G')
        beg, end = gen_coords(clv, '-', window=10)
        self.assertEqual((beg, end), (clv, clv + 10))
        self.assertEqual(fetch_seq(self.refseq, chrm, beg, end), 'GCCACTTTAT')


class TestSearchReferenceGenome(unittest.TestCase):
    def setUp(self):
        self.refseq = pysam.FastaFile(REF_FA)

    def test_chr12_DRAM_plus_strand(self):
        chrom = '12'
        clv = 102316878 - 1
        # confirm the corresponding seq in hg19
        self.assertEqual(
            fetch_seq(self.refseq, '12', 102316872, 102316872 + 6),
            'aataaa'
        )
        self.assertEqual(
            search_reference_genome(self.refseq, chrom, clv, '+', 50),
            ('AATAAA', 1, 102316872)
        )

    def test_chr19_AKT2_minus_strand(self):
        chrom = '19'
        clv = 40737005 - 1
        # confirm the corresponding seq in hg19
        self.assertEqual(
            fetch_seq(self.refseq, chrom, 40737011 + 1 - 6, 40737011 + 1),
            'TTTATT'
        )
        self.assertEqual(
            search_reference_genome(self.refseq, chrom, clv, '-', 50),
            ('AATAAA', 1, 40737011)
        )


class TestPysamFetch(unittest.TestCase):
    def setUp(self):
        self.refseq = pysam.FastaFile(
            os.path.join(os.path.dirname(__file__), 'mock_seq.fa'))

    def test_fetch_seq(self):
        # confirm the coordinate system in pysam is 0-based
        self.assertEqual(self.refseq.fetch('chr_mock1', 0, 0), '')
        self.assertEqual(self.refseq.fetch('chr_mock1', 0, 1), 'A')
        self.assertEqual(self.refseq.fetch('chr_mock1', 0, 2), 'AA')
        self.assertEqual(self.refseq.fetch('chr_mock1', 0, 3), 'AAT')

        self.assertEqual(self.refseq.fetch('chr_mock2', 1, 3), 'TC')


if __name__ == "__main__":
    unittest.main()
