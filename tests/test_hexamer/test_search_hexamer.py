import unittest

from kleat.hexamer.search import (
    plus_search,
    minus_search,
    search,
    extract_seq
)


class TestSearchHexamer(unittest.TestCase):
    def test_plus_search(self):
        self.assertEqual(plus_search('GGGAATAAAG', 9), ('AATAAA', 16, 3))
        self.assertEqual(plus_search('GGGAATAAA', 9), ('AATAAA', 16, 4))
        self.assertEqual(plus_search('GGGAATAAAGG', 9), ('AATAAA', 16, 2))
        self.assertEqual(plus_search('GGGATTAAAGG', 9), ('ATTAAA', 15, 2))
        self.assertEqual(plus_search('GGGAATAA', 9), None)

        self.assertEqual(plus_search('GAATAAAC', 10), ('AATAAA', 16, 4))
        self.assertEqual(plus_search('GGGGCTAC', 20), ('GGGGCT', 1, 13))
        self.assertEqual(plus_search('GTTTATTC', 6), None)

    def test_plus_search_lowercase(self):
        seq = 'GAATaaaC'
        #       4567890
        #             1
        self.assertEqual(plus_search(seq, 10), ('AATAAA', 16, 4))

    def test_plus_search_take_right_most_hexamer(self):
        self.assertEqual(plus_search('CAATAAANAATAAAC', 200), ('AATAAA', 16, 194))

    def test_plus_search_take_right_most_hexamer_with_Ns(self):
        self.assertEqual(plus_search('GCATTAAAAATNAAC', 200), ('ATTAAA', 15, 188))

    def test_plus_search_take_the_strongest_hexamer(self):
        self.assertEqual(plus_search('GCAATAAAATTAAAC', 200), ('AATAAA', 16, 188))

    def test_minus_search(self):
        seq = 'ATTTATTCCC'
        #      90123456789 <- one coord
        #       1          <- ten coord
        self.assertEqual(minus_search(seq, 9), ('AATAAA', 16, 15))
        seq = 'ATTTAATCCC'
        #      90123456789 <- one coord
        #       1          <- ten coord
        self.assertEqual(minus_search(seq, 9), ('ATTAAA', 15, 15))
        self.assertEqual(minus_search('GTTTATTC', 1), ('AATAAA', 16, 7))
        self.assertEqual(minus_search('ATCGTATATTGC', 5), ('AATATA', 10, 14))

    def test_minus_search_lowercase(self):
        self.assertEqual(minus_search('GTttattc', 1), ('AATAAA', 16, 7))

    def test_minus_search_take_left_most_hexamer(self):
        self.assertEqual(minus_search('GTTTATTTTTATTCG', 10), ('AATAAA', 16, 16))

    def test_minus_search_take_left_most_hexamer_with_Ns(self):
        self.assertEqual(minus_search('GTTTATTNTTTATTNNNTGTATTCG', 10), ('AATAAA', 16, 16))

    def test_minus_search_take_the_strongest_hexamer(self):
        self.assertEqual(minus_search('GTTTAATNTTTATTNNNTGTATTCG', 20), ('AATAAA', 16, 33))

    def test_minus_search_take_the_strongest_hexamer_in_lower_case(self):
        self.assertEqual(minus_search('gtttaatntttattnnntgtattcg', 20), ('AATAAA', 16, 33))


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
        self.assertEqual(search('+', clv, seq, 50), ('AATAAA', 16, 2))

    def test_minus_strand(self):
        """
         GGTTTATT
        0123456789 <-genome coord
         |      |
         clv    PAS
        """
        seq = 'GGTTTATT'
        clv = 1
        self.assertEqual(search('-', clv, seq, 50), ('AATAAA', 16, 8))



# Good drawing example, utilize them later
# def test_extract_seq_where_for_plus_strand_clv_supported_by_suffix():
#     """
#            AATAAA     AA   <-tail of suffix contig
#        ACGG┘||||└CGGCC┘    <-suffix contig
#        0123456789012345    <-contig coord
#               1      |
#     ...7890123456789012... <-genome coord
#           1         2|
#                      ^ref_clv
#     """
#     clv = 11
#     strand = '+'
#     contig = MagicMock()
#     contig.query_sequence = 'ACGGAATAAACGGCCAA'
#     contig.cigartuples = ((S.BAM_CMATCH, 15), (S.BAM_CSOFT_CLIP, 2))
#     ref_fa = MagicMock()
#     assert extract_seq(contig, strand, clv, ref_fa) == 'ACGGAATAAACGGCC'


# def test_extract_seq_where_for_minus_strand_clv_supported_by_suffix():
#     """
#      TTT  TTTATT        <-tail of suffix contig
#        └AC┘||||└CGGC    <-suffix contig
#         012345678901    <-contig coord
#         |         1
#      ...890123456789... <-genome coord
#         | 1
#         ^ref_clv
#     """
#     clv = 11
#     strand = '+'
#     contig = MagicMock()
#     contig.query_sequence = 'TTACTTTATTCGC'
#     contig.cigartuples = ((S.BAM_CMATCH, 15), (S.BAM_CSOFT_CLIP, 2))
#     ref_fa = MagicMock()
#     assert extract_seq(contig, strand, clv, ref_fa) == 'ACTTTATTCGC'

