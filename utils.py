from settings import HEADER
from settings import BAM_CSOFT_CLIP, BAM_CMATCH


def gen_clv_key_tuple(seqname, strand, clv):
    return (seqname, strand, clv)


def gen_clv_key_str(seqname, strand, clv):
    return f'{seqname}|{strand}|{clv}'


def write_row(clv_record, csvwriter):
    csvwriter.writerow([getattr(clv_record, _) for _ in HEADER])


def infer_contig_abs_ref_start(contig):
    """
    infer the absolute reference starting position taking into consideration
    the non-M bases (esp. softclipped bases)
    """
    pos = contig.reference_start
    for key, val in contig.cigartuples:
        if key != BAM_CMATCH:
            pos -= val
        break
    return pos


def infer_contig_abs_ref_end(contig):
    """
    infer the absolute reference starting position taking into consideration
    the non-M bases (esp. softclipped bases)
    """
    pos = contig.reference_end
    for key, val in reversed(contig.cigartuples):
        if key != BAM_CMATCH:
            pos += val
        break
    return pos


"""
Below are utility functions here apply to both contig and read as long as
they have a tail
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