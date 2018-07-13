"""
This script is adopted from
https://raw.githubusercontent.com/bcgsc/tasrkleat/0a6214bd6db60f947789ae3470b1bad5634db23c/hexamer_search/search_hexamer.py
which search for PAS hexamer on reference genome. However, in this adopted
version, it search for hexamer inside a given sequence (e.g contig) independent
of a reference genome
"""

from Bio import Seq

import kleat.misc.settings as S
from kleat.misc.apautils import fetch_seq
from kleat.hexamer import xseq_plus, xseq_minus


def reverse_complement(seq):
    return str(Seq.Seq(seq).reverse_complement().upper())
    # TODO: if prefer to drop dependency on biopython, test this function
    # thoroughly
    # return seq.translate(COMPLEMENT_DICT)[::-1]


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
    for (hmr, hid) in S.CANDIDATE_HEXAMERS:
        idx = seq.rfind(hmr)
        if idx > -1:
            return hmr, hid, idx + left_coord


def minus_search(seq, left_coord):
    """
    :param left_coord: the coordinate of the leftmost base in `seq`
    """
    seq = reverse_complement(seq.upper())
    for (hmr, hid) in S.CANDIDATE_HEXAMERS:
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


def search_ref_genome(refseq, chrom, clv, strand, window=50):
    """
    Different from search, this function search hexamer on reference genome

    :param refseq: an object returned by pysam.FastaFile, see TestSearch for
    its usage

    :param clv: 0-based. suppposed to be the 1-based coordinate of last base of
    3' UTR (e.g. output by KLEAT) - 1. converted to 0-based because pysam is
    0-based.

    return: a tuple of (hexamer, hexamer id (indicates strength), 0-based hexamer location)
    """
    beg, end = gen_coords(clv, strand, window)
    seq = fetch_seq(refseq, chrom, beg, end)
    # -1 as it's 0-based
    res = search_hexamer(seq, strand, beg, end - 1)
    return res


# TODO: remove default value for window
def extract_seq(contig, strand, ref_clv, ref_fa, ctg_clv, window=50):
    """
    extract the upstream sequence for PAS hexamer search

    :param ref_clv: would be necessary needs to combine sequence from both
                    contig and reference genome
    :param ctg_clv: the position of clv in contig coordinate.
    """
    seqname = contig.reference_name
    ctg_seq = contig.query_sequence
    cigartuples = contig.cigartuples

    args = cigartuples, ctg_seq, seqname, strand, ctg_clv, ref_clv, ref_fa, window
    if strand == '+':
        return xseq_plus.extract(*args)
    elif strand == '-':
        return xseq_minus.extract(*args)
    else:
        raise ValueError('unknown strand: "{0}"'.format(strand))


def gen_contig_hexamer_tuple(contig, strand, ref_clv, ref_fa, ctg_clv):
    """
    search PAS hexamer in contig, this ONLY works for suffix and link as the
    ref_clv nees to be at an end of the contig (not including the clipped bases)
    """
    # no need to reverse_complement the seq as the hexamer search function is
    # designed to search reference genome sequence, rev_comp is taken care of
    # within the function
    ctg_seq = extract_seq(contig, strand, ref_clv, ref_fa, ctg_clv)
    ctg_hex_tuple = search(strand, ref_clv, ctg_seq)

    if ctg_hex_tuple is not None:
        return ctg_hex_tuple
    else:
        return 'NA', -1, -1     # ctg_hex, ctg_hex_id, ctg_hex_pos


def gen_reference_hexamer_tuple(ref_fa, chrom_name, strand, ref_clv):
    """search PAS hexamer in reference genome"""
    na_tuple = 'NA', -1, -1     # ref_hex, ref_hex_id, ref_hex_pos
    ref_hex_tuple = search_ref_genome(ref_fa, chrom_name, ref_clv, strand)

    if ref_hex_tuple is None:
        return na_tuple
    else:
        return ref_hex_tuple
