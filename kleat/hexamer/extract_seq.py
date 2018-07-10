import kleat.misc.settings as S

from kleat.hexamer import xseq_plus, xseq_minus

"""Extract relevatn sequence, in which PAS hexamer is searched"""


def extract_seq(contig, strand, ref_clv, ref_fa, ctg_clv, window=50):
    """remove clipped ends before searching for hexamer, the clipped ends would
    affect calculation of genomics coordinates of the found PAS hexamer

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
