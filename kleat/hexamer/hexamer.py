from kleat.hexamer.search import search, search_ref_genome
from kleat.misc import apautils
from kleat.hexamer import xseq_plus, xseq_minus


# TODO: remove default value for window
def extract_seq(contig, strand, ref_clv, ref_fa, ctg_clv, window=50):
    """
    extract the upstream sequence for PAS hexamer search

    :param ref_clv: would be necessary needs to combine sequence from both
                    contig and reference genome
    :param ctg_clv: the position of clv in contig coordinate.
    """
    seqname = contig.reference_name

    # TODO: this can be done in one step, but needs update a lot of tests
    if apautils.is_hardclipped(contig):
        ctg_seq = apautils.infer_query_sequence(contig, always=True)
    else:
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
