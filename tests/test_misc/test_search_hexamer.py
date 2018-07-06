import unittest

from kleat.misc.search_hexamer import (
    plus_search,
    minus_search,
    search,
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
