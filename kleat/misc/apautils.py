from .settings import BAM_CMATCH, BAM_CSOFT_CLIP, HEADER


def gen_clv_key_tuple(seqname, strand, clv):
    return (seqname, strand, clv)


def gen_clv_key_str(seqname, strand, clv):
    return '{seqname}|{strand}|{clv}'.format(**locals())


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
        return suffix_segment.reference_start + 1
    elif tail_side == 'right':
        return suffix_segment.reference_end
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

    if the_cigar[0] != BAM_CSOFT_CLIP:
        cigar_idx = 'first' if tail_side == 'left' else 'last'
        raise ValueError('this may not be a {0} tailed segment as its '
                         '{1} CIGAR is not BAM_CSOFT_CLIP ({2})'.format(
                             tail_side, cigar_idx, BAM_CSOFT_CLIP))
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
    csvwriter.writerow([getattr(clv_record, _) for _ in HEADER])
