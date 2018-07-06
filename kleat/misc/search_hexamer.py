"""
This script is adopted from
https://raw.githubusercontent.com/bcgsc/tasrkleat/0a6214bd6db60f947789ae3470b1bad5634db23c/hexamer_search/search_hexamer.py
which search for PAS hexamer on reference genome. However, in this adopted
version, it search for hexamer inside a given sequence (e.g contig) independent
of a reference genome
"""


from kleat.misc.settings import CANDIDATE_HEXAMERS, COMPLEMENT_DICT


def reverse_complement(seq):
    return seq.translate(COMPLEMENT_DICT)[::-1]


def gen_coords(clv, strand, window=50):
    """
    generate the coordinates to be used for cutoff of searching range
    """
    if strand == '+':
        beg = clv - window + 1
        # as clv is the last based of 3'UTR and 0-based, so it should be
        # included in search
        end = clv + 1
    elif strand == '-':
        beg = clv
        end = clv + window
    else:
        raise ValueError('unknown strand: {0}'.format(strand))
    return beg, end


def plus_search(seq, right_coord):
    """
    :param right_coord: the coordinate of the rightmost base in `seq`
    """
    seq = seq.upper()
    left_coord = right_coord - len(seq) + 1
    for (hmr, hid) in CANDIDATE_HEXAMERS:
        idx = seq.rfind(hmr)
        if idx > -1:
            return hmr, hid, idx + left_coord


def minus_search(seq, left_coord):
    """
    :param left_coord: the coordinate of the leftmost base in `seq`
    """
    seq = reverse_complement(seq.upper())
    for (hmr, hid) in CANDIDATE_HEXAMERS:
        idx = seq.rfind(hmr)
        if idx > -1:
            return hmr, hid, len(seq) - idx + left_coord - 1


def search_hexamer(region, strand, beg, end):
    """search for PAS hexamer in the region"""
    if strand == '+':
        return plus_search(region, end)
    elif strand == '-':
        return minus_search(region, beg)
    else:
        raise ValueError('unknown strand: {0}'.format(strand))


def search(strand, clv, seq, window=50):
    """
    :param refseq: an object returned by pysam.FastaFile, see TestSearch for
    its usage

    :param clv: 0-based. suppposed to be the 1-based coordinate of last base of
    3' UTR (e.g. output by KLEAT) - 1. converted to 0-based because pysam is
    0-based.

    return: a tuple of (hexamer, hexamer id (indicates strength), 0-based hexamer location)
    """
    beg, end = gen_coords(clv, strand, window)
    # -1 as it's 0-based
    res = search_hexamer(seq, strand, beg, end - 1)
    return res
