from Bio import Seq

import kleat.misc.settings as S
from kleat.misc.calc_genome_offset import calc_genome_offset


def reverse_complement(seq):
    return str(Seq.Seq(seq).reverse_complement().upper())
    # TODO: if prefer to drop dependency on biopython, test this function
    # thoroughly
    # return seq.translate(COMPLEMENT_DICT)[::-1]


def gen_clv_key_tuple(seqname, strand, clv):
    return (seqname, strand, clv)


def gen_clv_key_tuple_from_clv_record(clv_record):
    crec = clv_record
    return gen_clv_key_tuple(crec.seqname, crec.strand, crec.clv)


def gen_clv_key_str(seqname, strand, clv):
    return '{seqname}|{strand}|{clv}'.format(**locals())


def fetch_seq(pysam_fa, seqname, beg, end):
    """
    In addition to pysam_fa.fetch, this wrapper handles fetching seq from circular
    DNA (e.g. chrM), too.

    :param pysam_fa: a pysam.libcalignmentfile.AlignmentFile instance
    """
    seq_len = pysam_fa.get_reference_length(seqname)

    circular_contigs = {'chrM', 'MT'}
    if beg >= 0:
        if end <= seq_len:
            res = pysam_fa.fetch(seqname, beg, end)
        else:
            if seqname in circular_contigs:
                res = (pysam_fa.fetch(seqname, beg, seq_len) +
                       pysam_fa.fetch(seqname, 0, end - seq_len))
            else:
                res = pysam_fa.fetch(seqname, beg, seq_len)
    else:
        if end <= seq_len:
            if seqname in circular_contigs:
                res = (pysam_fa.fetch(seqname, seq_len + beg, seq_len) +
                       pysam_fa.fetch(seqname, 0, end))
            else:
                res = pysam_fa.fetch(seqname, 0, end)
        else:
            raise NotImplementedError('How is this possible? Please report'
                                      'seqname: {0}, '
                                      'seq_length: {1}, '
                                      'beg: {2}, '
                                      'end: {3}')
    return res


def is_hardclipped(contig):
    for (key, val) in contig.cigartuples:
        if key == S.BAM_CHARD_CLIP:
            return True
    return False


# TODO: test this
def infer_query_sequence(contig, always=False):
    """
    :param always: mimic pysam api for infer_query_length, by setting always to
    True, it will try to infer sequence with hardclipped region based on XH
    tag. XH tag must exist, otherwise, it would fail loudly

    http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.infer_query_length

    """
    res = contig.query_sequence
    if always:
        num_cigars = len(contig.cigartuples)
        for k, (key, val) in enumerate(contig.cigartuples):
            if k == 0 and key == S.BAM_CHARD_CLIP:
                res = contig.get_tag('XH') + res
            elif k == num_cigars - 1:
                res += contig.get_tag('XH')
    return res


"""
Below are utility functions that apply to both contig and read as long as
they have a tail
"""


def right_tail(segment, tail_base='A'):
    """
    tseg: tail segment
    default tail_base applies main to alignment to genome, where the polyA tail
    strand is known
    """
    seq = segment.query_sequence
    last_cigar = segment.cigartuples[-1]
    # potential test case == "A0.R100820":
    return (
        seq.endswith(tail_base)
        # right clipped
        and last_cigar[0] == S.BAM_CSOFT_CLIP
        # clipped are all As
        and set(seq[-last_cigar[1]:]) == {tail_base}
    )


def left_tail(segment, tail_base='T'):
    seq = segment.query_sequence
    first_cigar = segment.cigartuples[0]
    # potential test case A0.S36384
    return (
        seq.startswith(tail_base)
        # left clipped
        and first_cigar[0] == S.BAM_CSOFT_CLIP
        # clipped are all Ts
        and set(seq[:first_cigar[1]]) == {tail_base}
    )


def has_tail(segment):
    if left_tail(segment):
        return 'left'
    elif right_tail(segment):
        return 'right'
    else:
        return None


def calc_ref_clv(suffix_segment, tail_side=None):
    """
    Calculate cleavage site position wst the reference

    :param tail_sideed: pass to avoid redundant checking of tail
    """
    if tail_side is None:
        tail_side = has_tail(suffix_segment)

    # the coordinates (+1 or not) are verified against visualization on IGV
    if tail_side == 'left':
        return suffix_segment.reference_start
    elif tail_side == 'right':
        return suffix_segment.reference_end - 1
    else:
        raise ValueError('{0} is not a suffix segment'.format(suffix_segment))


def calc_tail_length(suffix_segment, tail_side=None):
    """
    Calculate A/T length of a contig or a read, this information is extracted
    from softclip in the CIGAR
    """
    if tail_side is None:
        tail_side = has_tail(suffix_segment)

    if tail_side == 'left':
        the_cigar = suffix_segment.cigartuples[0]
    elif tail_side == 'right':
        the_cigar = suffix_segment.cigartuples[-1]
    else:
        raise ValueError('{0} is not a suffix segment'.format(suffix_segment))

    if the_cigar[0] != S.BAM_CSOFT_CLIP:
        cigar_idx = 'first' if tail_side == 'left' else 'last'
        raise ValueError('this may not be a {0} tailed segment as its '
                         '{1} CIGAR is not BAM_CSOFT_CLIP ({2})'.format(
                             tail_side, cigar_idx, S.BAM_CSOFT_CLIP))
    return the_cigar[1]


def calc_strand_from_suffix_segment(suffix_segment):
    """
    calculate the strand of clv (hence the corresponding gene) from a suffix
    segment
    """
    tail_side = has_tail(suffix_segment)
    if tail_side is None:
        raise ValueError('{0} is not a suffix segment, hence strand cannot be '
                         'inferred'.format(suffix_segment))
    return calc_strand(tail_side)


def calc_strand(tail_side):
    if tail_side == 'left':
        return '-'
    elif tail_side == 'right':
        return '+'
    else:
        raise ValueError('tail_side must be "left" or "right", '
                         'but {0} passed'.format(tail_side))


def write_row(clv_record, csvwriter):
    csvwriter.writerow([getattr(clv_record, _) for _ in S.HEADER])
