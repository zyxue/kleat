from settings import BAM_CSOFT_CLIP
"""
Utility functions here apply to both contig and read as long as they have a
tail
"""


def has_tail(segment):
    return right_tail(segment) or left_tail(segment)


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
        and last_cigar[0] == BAM_CSOFT_CLIP
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
        and first_cigar[0] == BAM_CSOFT_CLIP
        # clipped are all Ts
        and set(seq[:first_cigar[1]]) == {tail_base}
    )


def calc_ref_clv(suffix_segment):
    """
    Calculate cleavage site position wst the reference

    TODO: parse addition right_or_left to avoid double checking for right or
    left tail
    """
    if right_tail(suffix_segment):
        return suffix_segment.reference_end + 1
    elif left_tail(suffix_segment):
        return suffix_segment.reference_start
    else:
        raise ValueError(f'{suffix_segment} is not a suffix segment')


def calc_tail_length(suffix_segment):
    """
    Calculate A/T length of a contig or a read, this information is extracted
    from softclip in the CIGAR
    """
    if right_tail(suffix_segment):
        return suffix_segment.cigartuples[-1][1]
    elif left_tail(suffix_segment):
        return suffix_segment.cigartuples[0][1]
    else:
        raise ValueError(f'{suffix_segment} is not a suffix segment')


def calc_strand(suffix_segment):
    """
    calculate the strand of clv (hence the corresponding gene) this contig may
    support
    """
    if right_tail(suffix_segment):
        return '+'
    elif left_tail(suffix_segment):
        return '-'
    else:
        raise ValueError(f'{suffix_segment} is not a suffix segment')
